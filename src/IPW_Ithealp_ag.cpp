/* Modification of Therneau T (2015). _A Package for Survival Analysis in S_. version 2.38, <URL: https://CRAN.R-project.org/package=survival>.
 */

//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;


//[[Rcpp::export()]]
NumericMatrix IPW_ithealp_ag(NumericVector start, NumericVector stop, IntegerVector event, NumericVector eta, IntegerVector strata, IntegerVector sort1, IntegerVector sort2, NumericMatrix covar, NumericMatrix dpi,
  NumericVector weights, int nused, int nvar, int nalp) {
  int i, j, k, person;
  int indx1, istrat, p, p1;

  double dtime;
  double deadwt, denom, denom2;
  double risk;
  int nrisk, deaths;
  double temp;

  NumericVector a(nvar), a2(nvar);
  NumericVector aa0(nalp), aa02(nalp), risk1(nalp), deadwt2(nalp);
  NumericMatrix aa1(nvar, nalp), aa12(nvar, nalp);
  NumericMatrix Ithealp(nvar, nalp);
  IntegerVector keep(nused);

  indx1 = 0;
  deaths = 0;
  istrat = 0;

  for (person = 0; person < nused;) {
    if (person == strata[istrat]) {
      /* first subject in a new stratum */
      /* finish the work for the prior stratum */
      for (; indx1 < person; indx1++) {
        p1 = sort1[indx1];
        keep[p1] += deaths;
      }
      deaths = 0;
      istrat++;
    }
    p = sort2[person];
    keep[p] = -deaths;
    if (event[p]) {
      dtime = stop[p];
      for (person = person + 1; person < strata[istrat]; person++) {
        /* walk forward over any tied times */
        p = sort2[person];
        if (stop[p] != dtime) break;
        keep[p] = -deaths;
      }
      for (; indx1 < person; indx1++) {
        p1 = sort1[indx1];
        if (start[p1] < dtime) break;
        keep[p1] += deaths;
      }
      deaths++;
    } else person++;
  }
  for (; indx1 < nused; indx1++) {
    /* finish up the last strata */
    p1 = sort1[indx1];
    keep[p1] += deaths;
  }

  person = 0;
  indx1 = 0;
  istrat = 0;

  /* this next set is rezeroed at the start of each stratum */
  denom = 0;
  nrisk = 0;

  for (i = 0; i < nvar; i++) {
    a[i] = 0;
  }
  for (j = 0; j < nalp; j++) {
    aa0[j] = 0;
    for (i = 0; i < nvar; i++) {
      aa1(i, j) = 0;
    }
  }

  /* end of the per-stratum set */

  while (person < nused) {

    /* find the next deaSth time */
    for (k = person; k < nused; k++) {
      if (k == strata[istrat]) {
        /* hit a new stratum; reset temporary sums */
        istrat++;
        denom = 0;
        nrisk = 0;

        for (i = 0; i < nvar; i++) {
          a[i] = 0;
        }
        for (j = 0; j < nalp; j++) {
          aa0[j] = 0;
          for (i = 0; i < nvar; i++) {
            aa1(i, j) = 0;
          }
        }

        person = k; /* skip to end of stratum */
        indx1 = k;
      }
      p = sort2[k];
      if (event[p] == 1) {
        dtime = stop[p];
        break;
      }
    }
    if (k == nused) person = k; /* no more deaths to be processed */
    else {
       for (; indx1 < strata[istrat]; indx1++) {
        p1 = sort1[indx1];
        if (start[p1] < dtime) break;
        if (keep[p1] == 0) continue; /* skip any never-at-risk rows */
        nrisk--;
        if (nrisk == 0) {

          denom = 0;
          for (i = 0; i < nvar; i++) {
            a[i] = 0;
          }
          for (j = 0; j < nalp; j++) {
            aa0[j] = 0;
            for (i = 0; i < nvar; i++) {
              aa1(i, j) = 0;
            }
          }

        } else {

          risk = exp(eta[p1]) * weights[p1];
          for (j = 0; j < nalp; j++) {
            risk1[j] = pow(weights[p1], 2) * dpi(p1, j) * exp(eta[p1]);
          }

          denom -= risk;
          for (j = 0; j < nalp; j++)
            aa0[j] -= risk1[j];
          for (i = 0; i < nvar; i++) {
            a[i] -= risk * covar(p1, i);
            for (j = 0; j < nalp; j++)
              aa1(i, j) -= risk1[j] * covar(p1, i);
          }

        }

      }

      /*
       ** add any new subjects who are at risk
       ** denom2, a2, cmat2, deadwt and deaths count only the deaths
       */
      denom2 = 0;
      deadwt = 0;
      deaths = 0;
      for (j = 0; j < nalp; j++) {
        deadwt2[j] = 0;
      }
      for (j = 0; j < nalp; j++)
        aa02[j] = 0;
      for (i = 0; i < nvar; i++) {
        a2[i] = 0;
        for (j = 0; j < nalp; j++) {
          aa12(i, j) = 0;
        }
      }

      for (; person < strata[istrat]; person++) {
        p = sort2[person];
        if (stop[p] < dtime) break; /* no more to add */
        risk = exp(eta[p]) * weights[p];
        for (j = 0; j < nalp; j++) {
          risk1[j] = pow(weights[p], 2) * dpi(p, j) * exp(eta[p]);
        }

        if (event[p] == 1) {

          nrisk++;
          deaths++;
          denom2 += risk * event[p];
          deadwt += weights[p];

          for (j = 0; j < nalp; j++) {
            deadwt2[j] += pow(weights[p], 2) * dpi(p, j);
            aa02[j] += risk1[j];
          }

          for (i = 0; i < nvar; i++) {
            a2[i] += risk * covar(p, i);
            for (j = 0; j < nalp; j++) {
              Ithealp(i, j) += pow(weights[p], 2) * dpi(p, j) * covar(p, i);
              aa12(i, j) += risk1[j] * covar(p, i);
            }
          }
        } else if (keep[p] > 0) {
          nrisk++;
          denom += risk;
          for (j = 0; j < nalp; j++)
            aa0[j] += risk1[j];

          for (i = 0; i < nvar; i++) {
            a[i] += risk * covar(p, i);
            for (j = 0; j < nalp; j++)
              aa1(i, j) += risk1[j] * covar(p, i);
          }
        }
      }
      /*
       ** Add results into u and imat for all events at this time point
       */

      denom += denom2;

      for (j = 0; j < nalp; j++) {
        aa0[j] += aa02[j];
      }
      for (i = 0; i < nvar; i++) {
        a[i] += a2[i];
        temp = a[i] / denom; /*mean covariate at this time */

        for (j = 0; j < nalp; j++) {
          aa1(i, j) += aa12(i, j);
          Ithealp(i, j) -= deadwt2[j] * temp;
          Ithealp(i, j) -= deadwt * aa1(i, j) / denom;
          Ithealp(i, j) += deadwt * temp * aa0[j] / denom;
        }
      }

    }
  } /* end  of accumulation loop  */

  return (Ithealp);
}

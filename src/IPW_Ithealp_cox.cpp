//
//  IPW_Ithealp.cpp
//
//
//  Created by LeeJooyoung on 2019-01-17.
//

//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;


/* $Id: coxscore.c 11357 2009-09-04 15:22:46Z therneau $
 **
 ** Compute the score residuals for a Cox model
 **
 ** Input
 **      nx      number of subjects
 **      nvarx   number of variables in the covariance matrix
 **      y       matrix of time and status values
 **      strata  =1 for the last obs of each strata
 **      covar2  the matrix of covariates, rows=variables, columns=subjects
 **                (the S executive stores matrices in the Fortran ordering)
 **      score   the vector of subject scores, i.e., exp(beta*z)
 **      weights case weight
 **      method  ==1 for efron method
 **
 ** Output
 **      resid   a matrix of the same shape as x
 **
 ** Scratch
 **      scratch,  from which a and a2 are carved
 **
 ** Data must be sorted by strata, ascending time within strata, death before
 **                      censor within time.
 */

//[[Rcpp::export()]]
NumericMatrix IPW_ithealp_cox(NumericVector time, IntegerVector status, NumericVector lp, IntegerVector strata, NumericMatrix covar, NumericMatrix dpi,
  NumericVector weights, int nused, int nvar, int nalp) {
  int i, j, k, person;
  double dtime;
  double deadwt = 0, denom, denom2;
  double risk;
  double zbeta;
  int nrisk, ndead;
  double temp2;

  NumericVector a(nvar), a2(nvar);
  NumericVector aa0(nalp), aa02(nalp), risk1(nalp), deadwt2(nalp);
  NumericMatrix aa1(nvar, nalp), aa12(nvar, nalp);
  NumericMatrix Ithealp(nvar, nalp);

  denom = 0;
  nrisk = 0;
  for (i = 0; i < nvar; i++) {
    a[i] = 0;
    for (j = 0; j < nalp; j++)
      Ithealp(i, j) = 0;
  }

  for (person = nused - 1; person >= 0;) {
    if (strata[person] == 1) {
      denom = 0;
      nrisk = 0;
      for (i = 0; i < nvar; i++) {
        a[i] = 0;

      }
      for (k = 0; k < nalp; k++) {
        aa0[k] = 0;
        for (i = 0; i < nvar; i++) {
          aa1(i, k) = 0;
        }
      }
    }

    dtime = time[person];
    deadwt = 0;
    for (k = 0; k < nalp; k++) {
      deadwt2[k] = 0;
    }
    denom2 = 0;
    ndead = 0;

    while (person >= 0 && time[person] == dtime) {
      nrisk++;
      zbeta = lp[person];
      risk = weights[person] * exp(zbeta);
      for (k = 0; k < nalp; k++) {
        risk1[k] = pow(weights[person], 2) * dpi(person, k) * exp(zbeta);
      }
      if (status[person] == 0) {
        denom += risk;

        for (k = 0; k < nalp; k++)
          aa0[k] += risk1[k];
        for (i = 0; i < nvar; i++) {
          a[i] += risk * covar(person, i);
          for (k = 0; k < nalp; k++) {
            aa1(i, k) += risk1[k] * covar(person, i);

          }
        }
      } else {
        ndead++;
        denom2 += risk;
        deadwt += weights[person];
        for (k = 0; k < nalp; k++) {
          deadwt2[k] += pow(weights[person], 2) * dpi(person, k);
          aa02[k] += risk1[k];
        }
        for (i = 0; i < nvar; i++) {
          a2[i] += risk * covar(person, i);

          for (k = 0; k < nalp; k++) {
            Ithealp(i, k) += pow(weights[person], 2) * dpi(person, k) * covar(person, i);
            aa12(i, k) += risk1[k] * covar(person, i);
          }
        }
      }
      person--;
      if (person > 0 && strata[person] == 1) break;
    }

    if (ndead > 0) {
      /* add up terms*/
      denom += denom2;
      /* Rcout<<"person"<<person<<"\n";
       Rcout<<"denom"<<denom<<"\n";*/
      for (k = 0; k < nalp; k++) {
        aa0[k] += aa02[k];
      }

      for (i = 0; i < nvar; i++) {
        a[i] += a2[i];
        temp2 = a[i] / denom; /* mean */

        for (k = 0; k < nalp; k++) {
          aa1(i, k) += aa12(i, k);
          Ithealp(i, k) -= deadwt2[k] * temp2;
          Ithealp(i, k) -= deadwt * aa1(i, k) / denom;
          Ithealp(i, k) += deadwt * temp2 * aa0[k] / denom;
        }
      }

      denom2 = 0;
      for (k = 0; k < nalp; k++)
        aa02[k] = 0;
      for (i = 0; i < nvar; i++) {
        /*in anticipation */
        a2[i] = 0;
        for (k = 0; k < nalp; k++) aa12(i, k) = 0;
      }
    }

  }

  /* end  of accumulation loop  */

  return (Ithealp);
}

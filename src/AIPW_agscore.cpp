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
Rcpp::NumericMatrix AIPW_agscore_cpp(NumericVector start, NumericVector stop, IntegerVector event, NumericMatrix covar,
  IntegerVector eventid, IntegerVector id, NumericVector score, IntegerVector strata,
  IntegerMatrix marker, IntegerVector R, NumericMatrix pR,
  IntegerMatrix total_R, List marker_r, IntegerVector whereX, IntegerVector whereW,
  NumericVector gamma, NumericMatrix comb_y,
  int nvar, int n_marker, int nR, int ngamma, int nalp,
  bool second_cont_bl, bool second_cont_rr) {

  int i, j, k, l, person, pid, r, ty;
  double denom = 0, risk;
  int deaths;
  double meanwt;
  double hazard;
  double time;
  double zgamma, tmp_denom;
  int col1, col2;
  int nused, nX, nW, nstrat, nevent;
  int tmp_r;

  /* get local copies of some input args */
  nused = score.size();
  nstrat = strata.size();
  nX = whereX.size();
  nW = whereW.size();
  nevent = eventid.size();

  int two_y = comb_y.ncol();
  int ny = n_marker;
  int ny_rr = n_marker;
  if (second_cont_bl == TRUE) ny = ny + two_y;
  if (second_cont_rr == TRUE) ny_rr = ny + two_y;

  NumericVector a(nvar);
  NumericVector mean(nvar);
  NumericMatrix resid(nused, nvar);

  NumericVector tmp_y(ny);
  NumericVector tmp_yr(ny_rr);
  NumericVector tmp_num(ngamma);
  NumericVector tmp_w(ngamma);

  NumericMatrix Ecov(nR - 1, nvar);
  NumericVector EZ((nR - 1) * ngamma * nevent);

  for (person = 0; person < nused; person++) {

    if (event[person] == 1 && R[person] == 1) {
      for (i = 0; i < nevent; i++) {
        if (eventid[i] == id[person]) {
          pid = i;
          break;
        }
      }

      for (r = 1; r < nR; r++) {
        tmp_denom = 0;
        for (i = 0; i < ngamma; i++) {
          tmp_num[i] = 0;

        }

        NumericMatrix tmp_marker_r = marker_r[r - 1];
        for (ty = 0; ty < tmp_marker_r.nrow(); ty++) {
          for (k = 0; k < n_marker; k++) {
            tmp_y[k] = total_R(r, k) * marker(person, k) + (1 - total_R(r, k)) * tmp_marker_r(ty, k); /* General R for k=1,..., n_marker, r=1,..., nR*/
          }
          if (second_cont_bl == TRUE) {
            for (j = 0; j < two_y; j++) {
              col1 = comb_y(0, j) - 1;
              col2 = comb_y(1, j) - 1;
              tmp_y[n_marker + j] = tmp_y[col1] * tmp_y[col2];
            }
          }

          for (k = 0; k < n_marker; k++) {
            tmp_yr[k] = tmp_y[k];
          }
          if (second_cont_rr == TRUE) {
            for (j = 0; j < two_y; j++) {
              col1 = comb_y(0, j) - 1;
              col2 = comb_y(1, j) - 1;
              tmp_yr[n_marker + j] = tmp_y[col1] * tmp_y[col2];
            }
          }

          for (k = 0; k < ny; k++) {
            tmp_w[k] = tmp_y[k];
          }

          for (j = 0; j < nX; j++) {
            for (k = 0; k < ny_rr; k++) {
              tmp_w[ny + ny_rr * j + k] = covar(person, whereX[j] - 1) * tmp_yr[k];
            }
          }

          zgamma = 0;
          for (l = 0; l < ngamma; l++) {
            zgamma += gamma[l] * tmp_w[l];
          }
          tmp_denom += exp(zgamma);
          for (i = 0; i < ngamma; i++) {
            tmp_num[i] += tmp_w[i] * exp(zgamma);

          }

        }
        for (i = 0; i < ngamma; i++) {
          EZ[(r - 1) * ngamma * nevent + i * nevent + pid] = tmp_num[i] / tmp_denom;

        }

      }
    } else if (event[person] == 1 && R[person] != 1) {
      tmp_r = R[person] - 1;
      NumericMatrix tmp_marker_r = marker_r[tmp_r - 1];

      for (i = 0; i < nevent; i++) {
        if (eventid[i] == id[person]) {
          pid = i;
          break;
        }
      }

      tmp_denom = 0;
      for (i = 0; i < ngamma; i++) {
        tmp_num[i] = 0;

      }

      for (ty = 0; ty < tmp_marker_r.nrow(); ty++) {
        for (k = 0; k < n_marker; k++) {
          tmp_y[k] = total_R(tmp_r, k) * marker(person, k) + (1 - total_R(tmp_r, k)) * tmp_marker_r(ty, k); /* General R for k=1,..., n_marker, r=1,..., nR*/
        }
        if (second_cont_bl == TRUE) {
          for (j = 0; j < two_y; j++) {
            col1 = comb_y(0, j) - 1;
            col2 = comb_y(1, j) - 1;
            tmp_y[n_marker + j] = tmp_y[col1] * tmp_y[col2];
          }

        }
        for (k = 0; k < n_marker; k++) {
          tmp_yr[k] = tmp_y[k];
        }
        if (second_cont_rr == TRUE) {
          for (j = 0; j < two_y; j++) {
            col1 = comb_y(0, j) - 1;
            col2 = comb_y(1, j) - 1;
            tmp_yr[n_marker + j] = tmp_y[col1] * tmp_y[col2];
          }
        }

        for (k = 0; k < ny; k++) {
          tmp_w[k] = tmp_y[k];
        }

        for (j = 0; j < nX; j++) {
          for (k = 0; k < ny_rr; k++) {
            tmp_w[ny + ny_rr * j + k] = covar(person, whereX[j] - 1) * tmp_yr[k];
          }
        }

        zgamma = 0;
        for (l = 0; l < ngamma; l++) {
          zgamma += gamma[l] * tmp_w[l];
        }
        tmp_denom += exp(zgamma);
        for (i = 0; i < ngamma; i++) {
          tmp_num[i] += tmp_w[i] * exp(zgamma);

        }

        /*   Rcout<<tmp_denom<<"\n";
        for(i=0;i<nvar;i++)
        Rcout<<tmp_num[i]<<"\n"; */
      }
      for (i = 0; i < ngamma; i++) {
        EZ[(tmp_r - 1) * ngamma * nevent + i * nevent + pid] = tmp_num[i] / tmp_denom;

      }

    }
  }

  for (person = 0; person < nused;) {
    if (event[person] == 0) person++;
    else {

      deaths = 0;
      meanwt = 0;
      denom = 0;

      for (i = 0; i < nvar; i++) {
        a[i] = 0;
      }
      time = stop[person];
      for (k = person; k < nused; k++) {
        if (start[k] < time) {
          risk = score[k];
          denom += risk;
          for (i = 0; i < nvar; i++) {
            a[i] = a[i] + risk * covar(k, i);
          }
          if (stop[k] == time && event[k] == 1) {
            deaths++;
            meanwt += 1;
          }
        }
        if (strata[k] == 1) break;
      }

      hazard = meanwt / denom;
      for (i = 0; i < nvar; i++) mean[i] = a[i] / denom;
      for (k = person; k < nused; k++) {

        if (start[k] < time) {
          risk = score[k];
          for (i = 0; i < nvar; i++)
            resid(k, i) -= (covar(k, i) - mean[i]) * risk * hazard;
          if (stop[k] == time) {
            person++;
            if (event[k] == 1) {

              for (i = 0; i < nevent; i++) {
                if (eventid[i] == id[k]) {
                  pid = i;
                  break;
                }
              }

              for (r = 1; r < nR; r++) {
                for (l = 0; l < nX; l++) {
                  Ecov(r - 1, whereX[l] - 1) = covar(k, whereX[l] - 1);
                }
                for (l = 0; l < nW; l++) {
                  Ecov(r - 1, whereW[l] - 1) = covar(k, whereW[l] - 1);
                }
              }
              for (r = 1; r < nR; r++) {
                for (l = 0; l < ngamma; l++) {
                  Ecov(r - 1, nX + nW + l) = EZ[(r - 1) * ngamma * nevent + l * nevent + pid];
                }
              }

              if (R[k] == 1) {
                for (i = 0; i < nvar; i++) {
                  resid(k, i) += 1 / pR(pid, 0) * covar(k, i);
                  for (r = 1; r < nR; r++) {
                    resid(k, i) -= pR(pid, r) / pR(pid, 0) * Ecov(r - 1, i); /* pi_r/pi_1*E[Z|W_r]dNi(t) */
                  }
                }
              } else {
                for (i = 0; i < nvar; i++)
                  resid(k, i) += Ecov(R[k] - 2, i); /* E[Z|W_r]dNi(t) */
              }
              for (i = 0; i < nvar; i++)
                resid(k, i) -= mean[i];
            }
          }
        }
        if (strata[k] == 1) break;
      }
    }
  }

  return (resid);
}

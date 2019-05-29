
//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;


//[[Rcpp::export()]]
Rcpp::NumericMatrix AIPW_coxscore_cpp(NumericVector time, NumericVector status,
  NumericMatrix covar, IntegerVector eventid, IntegerVector id, NumericVector score, IntegerVector strata, IntegerMatrix marker, IntegerVector R, NumericMatrix pR,
  IntegerMatrix total_R, List marker_r, IntegerVector whereX,
  IntegerVector whereW,
  NumericVector gamma, NumericMatrix comb_y,
  int nvar, int n_marker, int nR, int ngamma, int nalp,
  bool second_cont_bl, bool second_cont_rr) {

  int i, j, k, l, person, pid;

  double denom = 0, risk;
  double temp, temp2;
  int deaths;
  double meanwt;
  double hazard;
  double zgamma;
  int col1, col2;
  int nused, nX, nW, nevent;
  NumericVector a(nvar);

  nused = score.size();
  nX = whereX.size();
  nW = whereW.size();
  nevent = eventid.size();

  int two_y = comb_y.ncol();
  int ny, ny_rr;
  ny = ny_rr = n_marker;
  if (second_cont_bl == TRUE) ny = ny + two_y;
  if (second_cont_rr == TRUE) ny_rr = ny + two_y;


  NumericMatrix resid(nused, nvar);
  double tmp_denom;
  NumericVector tmp_y(ny);
  NumericVector tmp_yr(ny_rr);
  NumericVector tmp_num(ngamma);
  NumericVector tmp_w(ngamma);

  NumericMatrix Ecov(nR - 1, nvar);
  NumericVector EZ((nR - 1) * ngamma * nused);

  int r, ty;
  int tmp_r;

  for (person = 0; person < nused; person++) {

    if (status[person] == 1 && R[person] == 1) {

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
          EZ[(r - 1) * ngamma * nused + i * nused + person] = tmp_num[i] / tmp_denom;
        }

      }
    } else if (status[person] == 1 && R[person] != 1) {

      tmp_r = R[person] - 1;
      NumericMatrix tmp_marker_r = marker_r[tmp_r - 1];

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
        EZ[(tmp_r - 1) * ngamma * nused + i * nused + person] = tmp_num[i] / tmp_denom;

      }

    }
  }

  /* Pass 1-- store the risk denominator in 'expect' */

  deaths = 0;
  meanwt = 0;
  denom = 0;
  strata[nused - 1] = 1;
  for (i = 0; i < nvar; i++) {
    a[i] = 0;
  }

  for (i = nused - 1; i >= 0; i--) {

    if (strata[i] == 1) {
      denom = 0;
      for (j = 0; j < nvar; j++) a[j] = 0;
    }

    risk = score[i];
    denom += risk;

    if (status[i] == 1) {
      deaths++;
      meanwt += 1;

    }
    for (j = 0; j < nvar; j++) {
      a[j] += risk * covar(i, j);
      resid(i, j) = 0;
    }

    if (deaths > 0 && (i == 0 || strata[i - 1] == 1 || time[i] != time[i - 1])) {
      /* last obs of a set of tied death times */

      hazard = meanwt / denom;

      for (j = 0; j < nvar; j++) {
        temp = (a[j] / denom); /* xbar */
        for (k = i; k < nused; k++) {

          temp2 = covar(k, j) - temp;

          if (time[k] == time[i] && status[k] == 1) {

            for (l = 0; l < nevent; l++) {
              if (eventid[l] == id[k]) {
                pid = l;
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
                Ecov(r - 1, nX + nW + l) = EZ[(r - 1) * ngamma * nused + l * nused + k];

              }
            }

            if (R[k] == 1) {
              resid(k, j) += 1 / pR(pid, 0) * covar(k, j);

              for (r = 1; r < nR; r++) {
                resid(k, j) -= pR(pid, r) / pR(pid, 0) * Ecov(r - 1, j); /* pi_r/pi_1*E[Z|W_r]dNi(t) */
              }
            } else {
              resid(k, j) += Ecov(R[k] - 2, j); /* E[Z|W_r]dNi(t) */
            }
            resid(k, j) -= temp;
          }
          resid(k, j) -= temp2 * score[k] * hazard; /* term 2 */

          if (strata[k] == 1) break;

        }
      }

      deaths = 0;
      meanwt = 0;

    }
  }

  return (resid);
}

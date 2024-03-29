// Modification of Therneau T (2015). _A Package for Survival Analysis in S_. version 2.38, <URL: https://CRAN.R-project.org/package=survival>.


//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;

//[[Rcpp::export()]]
Rcpp::List AIPW_coxfit_cpp(int maxiter, NumericVector time, IntegerVector status, NumericMatrix covar,
  IntegerVector eventid, IntegerVector id, NumericVector offset, IntegerVector strata, IntegerMatrix marker, IntegerVector R,
  NumericMatrix pR, List dpR,
  IntegerMatrix total_R, List marker_r, IntegerVector whereX,
  IntegerVector whereW,
  NumericVector gamma, NumericMatrix comb_y,
  int nvar, int n_marker, int nR, int ngamma, int nalp,
  double eps, bool first_cont_rr, bool second_cont_bl, bool second_cont_rr, NumericVector init_beta) {

  int i, j, k, l, h, person, pid, r, ty;

  double denom = 0, zbeta, risk;
  double temp, temp2;
  int ndead;
  double newlk = 0;
  double dtime;
  double deadwt;
  double denom2;
  double sctest;

  int halving;
  int nrisk;
  int conv;
  int col1, col2;
  int nused, nX, nW, nevent;
  int tmp_r;
  int doscale = 0;

  nused = offset.size();
  nX = whereX.size();
  nW = whereW.size();
  nevent = eventid.size();

  arma::mat imat(nvar, nvar), cmat(nvar, nvar), cmat2(nvar, nvar);
  NumericVector a(nvar), newbeta(nvar), a2(nvar);
  NumericVector u2(nvar);
  NumericVector score(nused), means(nvar);
  double zgamma, tmp_denom;

  int two_y = comb_y.ncol();
  int ny, ny_rr;
  ny = ny_rr = n_marker;
  if (second_cont_bl == TRUE) ny = ny + two_y;
  if (second_cont_rr == TRUE) ny_rr = ny + two_y;
  int tngamma;
  tngamma = ny + n_marker*nX;
  if (second_cont_rr == TRUE) ny + n_marker*nX + two_y*nX;
  

  NumericVector beta(nvar);
  NumericVector u(nvar);
  NumericMatrix Ithegam(nvar, ngamma);
  NumericMatrix Ithealp(nvar, nalp);
  NumericMatrix resid(nused, nvar);
  NumericVector loglik(2);
  int iter;

  NumericVector tmp_y(ny);
  NumericVector tmp_yr(ny_rr);
  NumericVector tmp_num(tngamma);
  NumericVector tmp_w(tngamma);
  NumericVector tmp_numw(tngamma);
  NumericMatrix tmp_mat(tngamma, tngamma);
  NumericMatrix Ecov(nR - 1, nvar);

  NumericVector EZ((nR - 1) * tngamma * nevent);
  NumericVector dEZ((nR - 1) * tngamma * nevent * tngamma);
  NumericVector dEcov((nR - 1) * nvar * tngamma);

  beta = init_beta;
  covar = clone(covar);


  for (person = 0; person < nused; person++) {

    if (status[person] == 1 && R[person] == 1) {

      for (i = 0; i < nevent; i++) {
        if (eventid[i] == id[person]) {
          pid = i;
          break;
        }
      }

      for (r = 1; r < nR; r++) {
        tmp_denom = 0;
        for (i = 0; i < tngamma; i++) {
          tmp_num[i] = 0;
          for (j = 0; j < tngamma; j++) {
            tmp_mat(i, j) = 0;
          }
        }
        for (j = 0; j < tngamma; j++) {
          tmp_numw[j] = 0;
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

          if(first_cont_rr == TRUE){
            for (j = 0; j < nX; j++) {
              for (k = 0; k < n_marker; k++) {
                tmp_w[ny + n_marker * j + k] = covar(person, whereX[j] - 1) * tmp_yr[k];
              }
            }
          }

          if(second_cont_rr == TRUE){
            for (j = 0; j < nX; j++) {
              for (k = 0; k < two_y; k++) {
                tmp_w[ny + n_marker * nX + two_y*j + k] = covar(person, whereX[j] - 1) * tmp_yr[n_marker + k];
              }
            }
          }

          zgamma = 0;
          for (l = 0; l < ngamma; l++) {
            zgamma += gamma[l] * tmp_w[l];
          }
          tmp_denom += exp(zgamma);
          for (i = 0; i < tngamma; i++) {
            tmp_num[i] += tmp_w[i] * exp(zgamma);
            for (j = 0; j < tngamma; j++) {
              tmp_mat(i, j) += tmp_w[i] * tmp_w[j] * exp(zgamma);
            }
          }

          for (j = 0; j < tngamma; j++) {
            tmp_numw[j] += tmp_w[j] * exp(zgamma);
          }

        }
        for (i = 0; i < tngamma; i++) {
          EZ[(r - 1) * tngamma * nevent + i * nevent + pid] = tmp_num[i] / tmp_denom;
          for (j = 0; j < tngamma; j++) {
            dEZ[(r - 1) * tngamma * nevent * tngamma + i * nevent * tngamma + pid * tngamma + j] = (tmp_mat(i, j) - tmp_num[i] / tmp_denom * tmp_numw[j]) / tmp_denom;
          }
        }

      }
    } else if (status[person] == 1 && R[person] != 1) {

      tmp_r = R[person] - 1;
      NumericMatrix tmp_marker_r = marker_r[tmp_r - 1];
      for (i = 0; i < nevent; i++) {
        if (eventid[i] == id[person]) {
          pid = i;
          break;
        }
      }

      tmp_denom = 0;
      for (i = 0; i < tngamma; i++) {
        tmp_num[i] = 0;
        for (j = 0; j < tngamma; j++) {
          tmp_mat(i, j) = 0;
        }
      }
      for (j = 0; j < tngamma; j++) {
        tmp_numw[j] = 0;
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

        if(first_cont_rr == TRUE){
          for (j = 0; j < nX; j++) {
            for (k = 0; k < n_marker; k++) {
              tmp_w[ny + n_marker * j + k] = covar(person, whereX[j] - 1) * tmp_yr[k];
            }
          }
        }

        if(second_cont_rr == TRUE){
          for (j = 0; j < nX; j++) {
            for (k = 0; k < two_y; k++) {
              tmp_w[ny + n_marker * nX + two_y*j + k] = covar(person, whereX[j] - 1) * tmp_yr[n_marker + k];
            }
          }
        }


        zgamma = 0;
        for (l = 0; l < ngamma; l++) {
          zgamma += gamma[l] * tmp_w[l];
        }
        tmp_denom += exp(zgamma);
        for (i = 0; i < tngamma; i++) {
          tmp_num[i] += tmp_w[i] * exp(zgamma);
          for (j = 0; j < tngamma; j++) {
            tmp_mat(i, j) += tmp_w[i] * tmp_w[j] * exp(zgamma);
          }
        }

        for (j = 0; j < tngamma; j++) {
          tmp_numw[j] += tmp_w[j] * exp(zgamma);
        }

       }
      for (i = 0; i < tngamma; i++) {
        EZ[(tmp_r - 1) * tngamma * nevent + i * nevent + pid] = tmp_num[i] / tmp_denom;
        for (j = 0; j < tngamma; j++) {
          dEZ[(tmp_r - 1) * tngamma * nevent * tngamma + i * nevent * tngamma + pid * tngamma + j] = (tmp_mat(i, j) - tmp_num[i] / tmp_denom * tmp_numw[j]) / tmp_denom;
        }
      }

    }
  }

  /*
   ** Subtract the mean from each covar, as this makes the regression
   **  much more stable.
   */

/* Do not scaling! */

  temp2 = nused;

  for (i=0; i<nvar; i++) {
    temp=0;
    for (person=0; person<nused; person++)
      temp += covar(person, i);
    temp /= temp2;
    means[i] = temp;
    for (person=0; person<nused; person++) covar(person, i) -=temp;

  }



  /*
   ** do the initial iteration step
   */
  strata[nused - 1] = 1;
  loglik[1] = 0;
  for (i = 0; i < nvar; i++) {
    u[i] = 0;
    a2[i] = 0;
    for (j = 0; j < nvar; j++) {
      imat(i, j) = 0;
      cmat2(i, j) = 0;
    }
  }

  for (person = nused - 1; person >= 0;) {
    if (strata[person] == 1) {
      nrisk = 0;
      denom = 0;
      for (i = 0; i < nvar; i++) {
        a[i] = 0;
        for (j = 0; j < nvar; j++) cmat(i, j) = 0;
      }
    }
    dtime = time[person];
    ndead = 0; /*number of deaths at this time point */
    deadwt = 0; /* sum of weights for the deaths */
    denom2 = 0; /* sum of weighted risks for the deaths */
    while (person >= 0 && time[person] == dtime) {
      /* walk through the this set of tied times */
      nrisk++;
      zbeta = offset[person]; /* form the term beta*z + offset */
      for (i = 0; i < nvar; i++)
        zbeta += beta[i] * covar(person, i);
      risk = exp(zbeta);
      if (status[person] == 0) {
        denom += risk;
        /* a contains weighted sums of x, cmat sums of squares */
        for (i = 0; i < nvar; i++) {
          a[i] += risk * covar(person, i);
          for (j = 0; j <= i; j++)
            cmat(i, j) += risk * covar(person, i) * covar(person, j);
        }
      } else {
        for (i = 0; i < nevent; i++) {
          if (eventid[i] == id[person]) {
            pid = i;
            break;
          }
        }

        ndead++;
        deadwt += 1;
        denom2 += risk;

        for (r = 1; r < nR; r++) {
          for (l = 0; l < nX; l++) {
            Ecov(r - 1, whereX[l] - 1) = covar(person, whereX[l] - 1);
          }
          for (l = 0; l < nW; l++) {
            Ecov(r - 1, whereW[l] - 1) = covar(person, whereW[l] - 1);
          }
        }

        for (r = 1; r < nR; r++) {
          for (i = 0; i < tngamma; i++) {
            Ecov(r - 1, nX + nW + i) = EZ[(r - 1) * tngamma * nevent + i * nevent + pid];
          }
        }

        for (i = 0; i < nvar; i++) {
          if (R[person] == 1) {
            u[i] += 1 / pR(pid, 0) * covar(person, i); /* 1/pi_1*Z*dNi(t)*/
            loglik[1] += 1 / pR(pid, 0) * covar(person, i) * beta[i]; /* sum 1/pi_1*Z*beta*dNi(t)*/
            for (r = 1; r < nR; r++) {
              u[i] -= pR(pid, r) / pR(pid, 0) * Ecov(r - 1, i); /* pi_r/pi_1*E[Z|W_r]dNi(t) */
              loglik[1] -= pR(pid, r) / pR(pid, 0) * Ecov(r - 1, i) * beta[i]; /* sum pi_r/pi_1*E[Z|w_r]*beta*dNi(t) */
            }
          } else {
            u[i] += Ecov(R[person] - 2, i); /* E[Z|W_r]dNi(t) */
            loglik[1] += Ecov(R[person] - 2, i) * beta[i]; /* E[Z|W_r]*beta*dNi(t) */
          }
          a2[i] += risk * covar(person, i);
          for (j = 0; j <= i; j++)
            cmat2(i, j) += risk * covar(person, i) * covar(person, j);
        }

      }
      person--;
      if (person >= 0 && strata[person] == 1) break; /*ties don't cross strata */
    }
    if (ndead > 0) {
      /* we need to add to the main terms */

      denom += denom2;
      loglik[1] -= deadwt * log(denom);

      for (i = 0; i < nvar; i++) {
        a[i] += a2[i];
        temp2 = a[i] / denom; /* mean  S(1)/S(0) */
        u[i] -= deadwt * temp2; /*   -S(1)/S(0)*dNi(t) -> for all individual regardless of missing  */
        for (j = 0; j <= i; j++) {
          cmat(i, j) += cmat2(i, j);
          imat(j, i) += deadwt * (cmat(i, j) - temp2 * a[j]) / denom; /* upper triangle */

        }
      }

      for (i = 0; i < nvar; i++) {
        a2[i] = 0;
        for (j = 0; j < nvar; j++) cmat2(i, j) = 0;
      }
    }
  }

  /* end  of accumulation loop */
  loglik[0] = loglik[1]; /* save the loglik for iter 0 */

  /* am I done?
   **   update the betas and test for convergence
   */
  for (i = 0; i < nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
    a[i] = u[i];

  for (i = 0; i < nvar; i++) {
    for (j = 0; j < i; j++) {
      imat(i, j) = imat(j, i);
    }
  }
  imat = arma::inv(imat); /* Inverse matrix */

  for (i = 0; i < nvar; i++) {
    a[i] = a[i] * imat(i, i);
  }

  sctest = 0;
  for (i=0; i<nvar; i++)
    sctest +=  u[i]*a[i]; /* score test */

  /*
   **  Never, never complain about convergence on the first step.  That way,
   **  if someone HAS to they can force one iter at a time.
   ** A non-finite loglik comes from exp overflow and requires almost
   **  malicious initial values.
   */
  for (i = 0; i < nvar; i++) {
    newbeta[i] = beta[i] + a[i];
  }

  if (maxiter == 0 || is_na(loglik)[0] || 0 != is_infinite(loglik)[0]) {
    conv = 0;
    goto finish;
  }
  /*
   ** here is the main loop
   */
  halving = 0; /* =1 when in the midst of "step halving" */
  for (iter = 1; iter <= maxiter;
    (iter) ++) {
    R_CheckUserInterrupt();
    newlk = 0;
    for (i = 0; i < nvar; i++) {
      u[i] = 0;
      for (j = 0; j < nvar; j++)
        imat(i, j) = 0;
    }

    /*
     ** The data is sorted from smallest time to largest
     ** Start at the largest time, accumulating the risk set 1 by 1
     */
    for (person = nused - 1; person >= 0;) {
      if (strata[person] == 1) {
        denom = 0;
        nrisk = 0;
        for (i = 0; i < nvar; i++) {
          a[i] = 0;
          for (j = 0; j < nvar; j++) cmat(i, j) = 0;
        }
      }

      dtime = time[person];
      deadwt = 0;
      ndead = 0;
      denom2 = 0;
      while (person >= 0 && time[person] == dtime) {
        nrisk++;
        zbeta = offset[person];
        for (i = 0; i < nvar; i++)
          zbeta += newbeta[i] * covar(person, i);
        risk = exp(zbeta);

        if (status[person] == 0) {
          denom += risk;

          for (i = 0; i < nvar; i++) {
            a[i] += risk * covar(person, i);
            for (j = 0; j <= i; j++)
              cmat(i, j) += risk * covar(person, i) * covar(person, j);
          }
        } else {
          for (i = 0; i < nevent; i++) {
            if (eventid[i] == id[person]) {
              pid = i;
              break;
            }
          }

          ndead++;
          denom2 += risk;
          deadwt += 1;
          for (r = 1; r < nR; r++) {
            for (l = 0; l < nX; l++) {
              Ecov(r - 1, whereX[l] - 1) = covar(person, whereX[l] - 1);
            }
            for (l = 0; l < nW; l++) {
              Ecov(r - 1, whereW[l] - 1) = covar(person, whereW[l] - 1);

            }
          }
          for (r = 1; r < nR; r++) {
            for (i = 0; i < tngamma; i++) {
              Ecov(r - 1, nX + nW + i) = EZ[(r - 1) * tngamma * nevent + i * nevent + pid] - means[nX + nW + i];
            }
          }

          for (i = 0; i < nvar; i++) {
            if (R[person] == 1) {
              u[i] += 1 / pR(pid, 0) * covar(person, i); /* 1/pi_1*Z*dNi(t)*/
              newlk += 1 / pR(pid, 0) * covar(person, i) * newbeta[i]; /* sum 1/pi_1*Z*dNi(t)*/
              for (r = 1; r < nR; r++) {
                u[i] -= pR(pid, r) / pR(pid, 0) * Ecov(r - 1, i); /* pi_r/pi_1*E[Z|W_r]dNi(t) */
                newlk -= pR(pid, r) / pR(pid, 0) * Ecov(r - 1, i) * newbeta[i];
              }
            } else {
              u[i] += Ecov(R[person] - 2, i); /* E[Z|W_r]dNi(t) */
              newlk += Ecov(R[person] - 2, i) * newbeta[i]; /* E[Z|W_r]*beta*dNi(t) */
            }

            a2[i] += risk * covar(person, i);
            for (j = 0; j <= i; j++)
              cmat2(i, j) += risk * covar(person, i) * covar(person, j);
          }
        }
        person--;
        if (person > 0 && strata[person] == 1) break; /*tied times don't cross strata*/
      }

      if (ndead > 0) {
        /* add up terms*/
        denom += denom2;
        newlk -= deadwt * log(denom);
        for (i = 0; i < nvar; i++) {
          a[i] += a2[i];
          temp2 = a[i] / denom; /* mean */
          u[i] -= deadwt * temp2;

          for (j = 0; j <= i; j++) {
            cmat(i, j) += cmat2(i, j);
            imat(j, i) += deadwt * (cmat(i, j) - temp2 * a[j]) / denom;
          }
        }

        denom2 = 0;
        for (i = 0; i < nvar; i++) {
          /*in anticipation */
          a2[i] = 0;
          for (j = 0; j < nvar; j++) cmat2(i, j) = 0;
        }
      }
    } /* end  of accumulation loop  */

              /* am I done?
     **   update the betas and test for convergence
     */
     if (fabs(1 - (loglik[1] / newlk)) <= eps && halving == 0) {
      /* all done */
      loglik[1] = newlk;

      for (i = 0; i < nvar; i++) {
        for (j = 0; j < i; j++) {
          imat(i, j) = imat(j, i);
        }
      }

      imat = arma::inv(imat); /* Inverse matrix */
      conv = 1;
      goto finish;
    }
    /*
     ** a non-finite loglik is very rare: a step so bad that we get
     ** an overflow of the exp function.
     **  When this happens back up one iteration and quit
     */

    NumericVector nnewlk(1);
    nnewlk[0] = newlk;

    if (is_nan(nnewlk)[0] || 0 != is_infinite(nnewlk)[0]) {
      for (i = 0; i < nvar; i++) newbeta[i] = beta[i];
       maxiter = iter+1;
      continue;
    }

    if (iter == maxiter) break; /*skip the step halving calc*/
    if (newlk < loglik[1]) {
      /*it is not converging ! */
      halving = 1;
      for (i = 0; i < nvar; i++)
        newbeta[i] = (newbeta[i] + beta[i]) / 2; /*half of old increment */
    } else {
      halving = 0;
      loglik[1] = newlk;
      u2 = u;

      for (i = 0; i < nvar; i++) {
        for (j = 0; j < i; j++) {
          imat(i, j) = imat(j, i);
        }
      }
      imat = arma::inv(imat); /* Inverse matrix */

      for (i = 0; i < nvar; i++) {
        u2[i] = u2[i] * imat(i, i);
      }

      j = 0;
      for (i = 0; i < nvar; i++) {
        beta[i] = newbeta[i];
        newbeta[i] = newbeta[i] + u2[i];
           }
    }
  } /* return for another iteration */

  /*
   ** We end up here only if we ran out of iterations
   */
  loglik[1] = newlk;
  for (i = 0; i < nvar; i++) {
    for (j = 0; j < i; j++) {
      imat(i, j) = imat(j, i);
    }
  }
  imat = arma::inv(imat); /* Inverse matrix */
  conv = 2;
  finish:

    /* Calculate Ithealp , Ithegam */
    /*
     ** The data is sorted from smallest time to largest
     ** Start at the largest time, accumulating the risk set 1 by 1
     */
    for (person = nused - 1; person >= 0;) {

      dtime = time[person];
      while (person >= 0 && time[person] == dtime) {

        if (status[person] == 1) {

          for (i = 0; i < nevent; i++) {
            if (eventid[i] == id[person]) {
              pid = i;
              break;
            }
          }

          NumericMatrix dpR_tmp = dpR[pid];

          for (r = 1; r < nR; r++) {
            for (l = 0; l < nX; l++) {
              Ecov(r - 1, whereX[l] - 1) = covar(person, whereX[l] - 1);
            }
            for (l = 0; l < nW; l++) {
              Ecov(r - 1, whereW[l] - 1) = covar(person, whereW[l] - 1);
            }
          }

          for (r = 1; r < nR; r++) {
            for (i = 0; i < tngamma; i++) {
              Ecov(r - 1, nX + nW + i) = EZ[(r - 1) * tngamma * nevent + i * nevent + pid] - means[nX + nW + i];
              for (k = 0; k < tngamma; k++) {
                dEcov[(r - 1) * nvar * tngamma + (nX + nW + i) * tngamma + k] = dEZ[(r - 1) * tngamma * nevent * tngamma + i * nevent * tngamma + pid * tngamma + k];
              }
            }
          }

          for (i = 0; i < nvar; i++) {
            if (R[person] == 1) {
              for (h = 0; h < nalp; h++) {
                Ithealp(i, h) += 1 / pow(pR(pid, 0), 2) * dpR_tmp(0, h) * covar(person, i);
                for (r = 1; r < nR; r++) {
                  Ithealp(i, h) += dpR_tmp(r, h) / pR(pid, 0) * Ecov(r - 1, i);
                  Ithealp(i, h) -= pR(pid, r) * dpR_tmp(0, h) / pow(pR(pid, 0), 2) * Ecov(r - 1, i);
                }
              }
              for (k = 0; k < ngamma; k++)
                for (r = 1; r < nR; r++) {
                  Ithegam(i, k) += pR(pid, r) / pR(pid, 0) * dEcov[(r - 1) * nvar * tngamma + i * tngamma + k];

                }

            } else {
              for (k = 0; k < ngamma; k++) {
                Ithegam(i, k) -= dEcov[(R[person] - 2) * nvar * tngamma + i * tngamma + k];
              }

            }
          }
        }
        person--;
        if (person > 0 && strata[person] == 1) break; /*tied times don't cross strata*/
      }

    } /* end  of accumulation loop  */

  /*
   ** create the output list
   */
  List rlist = List::create(Named("coef") = beta,
    Named("means") = means,
    Named("u") = u,
    Named("imat") = imat,
    Named("Ithegam") = Ithegam,
    Named("Ithealp") = Ithealp,
    Named("loglik") = loglik,
    Named("sctest") = sctest,
    Named("iter") = iter,
    Named("conv") = conv);

  return (rlist);
}


#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix dpR_multinom (NumericVector X, NumericVector score, int nR, bool twostage, int nrp, int nalp, int nalp0){

  int i, j, k;
  int ncp = (nalp - nalp0)/nrp;

  int fnR;

  if(twostage == TRUE) fnR = nR-1;
  else fnR = nR;

  NumericMatrix rep(nR, nalp);


  double sscore;
  double tmp;

  sscore = sum(score);
  for(i=0; i<nrp; i++){
    for(j=0; j<ncp; j++)
      rep(0, j*nrp + i) = - score[i]*X[j]/pow(1+sscore,2);
  }

  for(k=1; k<fnR; k++){
    for(i=0; i<nrp; i++){
      for(j=0; j<ncp; j++)
        if(i == (k-1)) rep(k, j*nrp + i) =  score[i]*X[j]*(1+sscore - score[i])/pow(1+sscore,2);
        else rep(k, j*nrp + i) = - score[k-1]*score[i]*X[j]/pow(1+sscore,2);
    }
  }



  return(rep);
}



#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector findR (NumericVector R, NumericVector event, int fonR, NumericMatrix Rmat, NumericMatrix ototal_R){

  int i, j, k;
  int n = event.size();
  LogicalVector v(Rmat.ncol());

  for (i=0; i<n; i++) {
    if (R[i] == 1 && event[i]==1) {
      for (j=0; j<fonR; j++) {
         v = (Rmat(i,_) == ototal_R(j, _));
         if(is_true(all(v))) R[i] = j+1;
        }

      }
    }

 return(R);
}
// [[Rcpp::export]]
Rcpp::NumericVector findcause (NumericVector R, NumericVector cause, NumericVector event, NumericMatrix marker, int on_subtype,  NumericMatrix ototal_subtype){

  int i, j, k;
  int n = event.size();


  LogicalVector v(marker.ncol());

  for (i=0; i<n; i++) {
    if(event[i]==1){
      if (R[i] != 1) {
        cause[i] = NA_REAL;
      } else {
        for (j =0; j<on_subtype; j++) {
          v = (marker(i,_) == ototal_subtype(j, _));
          if(is_true(all(v))) cause[i] = j+1;

         }
      }
    }
  }

  return(cause);

}


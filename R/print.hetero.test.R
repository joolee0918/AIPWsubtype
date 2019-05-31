
#' @export
#' @method print hetero.test
#' @keywords internal

print.hetero.test <- function(object){

  n <- length(object)
  cat("Heterogeneity Tests (Likelihood ratio test) ")
  cat("\n")
  for(i in 1:n){
    cat("---------------------------------------------")
    cat("\n")
    print(object[[i]], row.names = FALSE)

  }
}

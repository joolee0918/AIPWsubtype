
#' @method print hetero.test
#' @keywords internal

#' @export
print.loglik.test <- function(object){

  n <- length(object)
  pdig <- max(1, getOption("digits")-4)

  cat("\n")
  cat("Heterogeneity Tests (Likelihood ratio test) ")
  cat("\n")
  cat("--------------------------------------------")
  cat("\n")
  object <- rbindlist(object)
  object$pvalue <- format.pval(object$pvalue, digits=pdig)
  print(object, row.names = F)

  cat("\n")
}

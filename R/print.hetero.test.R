
#' @method print hetero.test
#' @keywords internal

#' @export
print.hetero.test <- function(object){

  pdig <- max(1, getOption("digits")-4)

  cat("\n")
  cat("Heterogeneity Tests (Likelihood ratio test) ")
  cat("\n")
  cat("--------------------------------------------")
  cat("\n")

  object$pvalue <- format.pval(object$pvalue, digits=pdig)

  print(object)
  invisible()
}

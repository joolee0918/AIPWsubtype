
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

  print(object[[1]], row.names = F)
  cat("\n")
  invisible()
}

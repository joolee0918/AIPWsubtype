# =========================================================
# Control List
# =========================================================

#' Control parameters

#' @param eps Tolerance for convergence in the parameters. Default is 1e-12.
#' @param maxiter The maxiumum number of interations. Default is 2000.

controlList <- function(control) {

  control.default <- list(eps = 1e-12, iter.max = 2000)

  control <- c(control)
  namec <- names(control)

  if (length(uname <- namec[!namec %in% names(control.default)]) > 0) {
    warning("\n unknown names in 'control': ", paste(uname, collapse = ", "))
  }

  control.default[namec] <- control

  return(control.default)

}


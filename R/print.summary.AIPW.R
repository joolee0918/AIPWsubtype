# =========================================================
# Print class summary.AIPWcprisk
# =========================================================

#' @export
#' @method print summary.AIPWcprisk
#' @keywords internal

print.summary.AIPWcprisk <- function(x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"),
    ...) {

    cat("\n")
    cat("#### AIPW semiparametric cause-specific competing risks model  ####")
    cat("\n\n")

    if (!is.null(x$call)) {

        cat("Call:\n")
        dput(x$call)
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" Coxreg failed.", x$fail, "\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    cat("  n=", x$n)
    if (!is.null(x$nevent))
        cat(", number of events=", x$nevent, "\n") else cat("\n")

    if (nrow(x$coef) == 0) {
        # Null model
        cat("   Null model\n")
        return()
    }


    if (!is.null(x$coefficients)) {
        cat("\n")
        printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
    }
    if (!is.null(x$conf.int)) {
        cat("\n")
        print(x$conf.int)
    }
    cat("\n")


    invisible()
}

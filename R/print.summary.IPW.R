# =========================================================
# Print class summary.IPWcprisk
# =========================================================

#' @export
#' @method print summary.IPWcprisk
#' @keywords internal

print.summary.IPWcprisk <- function(x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"),
    ...) {

    cat("\n")
    cat("#### IPW semiparametric cause-specific competing risks model  ####")
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
    cat("  complete-case n=", x$nc)
    if (!is.null(x$ncevent))
        cat(", number of complete-case events=", x$ncevent, "\n") else cat("\n")

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

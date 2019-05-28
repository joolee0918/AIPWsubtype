# =========================================================
# Print class IPWcprisk
# =========================================================

#' @export
#' @method print IPWcprisk
#' @keywords internal
print.IPWcprisk <- function(object, ...) {

    cat("\n")
    cat("#### IPW semiparametric cause-specific competing risks model  ####")
    cat("\n\n")

    if (!is.null(cl <- object$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }


    if (!is.null(object$fail)) {
        cat("Fitting failed.", object$fail, "\n")
        return()
    }

    digits = max(1L, getOption("digits") - 3L)
    signif.stars = FALSE
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    J <- object$neventtype

    coef <- object$coef
    se <- sqrt(diag(object$var))
    if (is.null(est) | is.null(se))
        stop("Input is not valid")

    tmp <- cbind(coef, exp(coef), se, coef/se, pchisq((coef/se)^2, 1, lower.tail = FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "z", "p"))


    printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE, signif.stars = signif.stars, ...)

    if (!is.null(object$nevent))
        cat(", number of events=", object$nevent, "\n")


    invisible(object)

}

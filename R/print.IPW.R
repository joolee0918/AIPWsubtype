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
    if (is.null(coef) | is.null(se))
        stop("Input is not valid")

    tmp <- cbind(coef, exp(coef), se, coef/se, pchisq((coef/se)^2, 1, lower.tail = FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "z", "p"))


    printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE, signif.stars = signif.stars, ...)
    logtest <- -2 * (object$loglik[1] - object$loglik[2])
    df <- sum(!is.na(coef))
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
        df, " df,", " p=",
        format.pval(pchisq(logtest, df, lower.tail=FALSE), digits=digits),
        "\n",  sep="")


    cat("  n=", object$n)
    if (!is.null(object$nevent))
      cat(", number of events=", object$nevent, "\n") else cat("\n")
    cat("  complete-case n=", object$nc)
    if (!is.null(object$ncevent))
      cat(", number of complete-case events=", object$ncevent, "\n") else cat("\n")


    invisible(object)

}

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

    pdig <- max(1, getOption("digits")-4)  # default it too high IMO
    cat("Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
        x$logtest["df"], " df,", "   p=",
        format.pval(x$logtest["pvalue"], digits=pdig),
        "\n", sep = "")
    cat("Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
        x$waldtest["df"], " df,", "   p=",
        format.pval(x$waldtest["pvalue"], digits=pdig),
        "\n", sep = "")
    cat("Score (logrank) test = ", format(round(x$sctest["test"], 2)), "  on ",
        x$sctest["df"]," df,", "   p=",
        format.pval(x$sctest["pvalue"], digits=pdig), sep ="")
    if (is.null(x$robscore))
      cat("\n\n")
    else cat(",   Robust = ", format(round(x$robscore["test"], 2)),
             "  p=",
             format.pval(x$robscore["pvalue"], digits=pdig), "\n\n", sep="")

    if (x$used.robust)
      cat("  (Note: the likelihood ratio and score tests",
          "assume independence of\n     observations within a cluster,",
          "the Wald and robust score tests do not).\n")

    invisible()
}

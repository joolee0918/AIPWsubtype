# =========================================================
# class summary.IPWcprisk
# =========================================================

#' @export
#' @method summary IPWcprisk

summary.IPWcprisk <- function(object, conf.int = 0.95, scale = 1, ...) {

    if (!is.null(object$fail)) {
        class(object) <- "summary.IPWcprisk"
        return(object)
    }


    beta <- object$coefficients
    if (is.null(object$coefficients)) {
        # Null model
        return(object)  #The summary method is the same as print in this case
    }
    nabeta <- !(is.na(beta))  #non-missing coefs
    beta2 <- beta[nabeta]

    if (is.null(beta) | is.null(object$var))
        stop("Input is not valid")
    se <- sqrt(diag(object$var))

    rval <- list(call = object$call, fail = object$fail, loglik = object$loglik, n = object$n, nc = object$nc,
        nevent = object$nevent, ncevent = object$ncevent)

    tmp <- cbind(beta, exp(beta), se, beta/se, pchisq((beta/se)^2, 1, lower.tail = FALSE))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
    rval$coefficients <- tmp
    if (conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se), exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)", paste("lower .", round(100 * conf.int,
            2), sep = ""), paste("upper .", round(100 * conf.int, 2), sep = "")))
        rval$conf.int <- tmp
    }

    df <- length(beta2)
    logtest <- -2 * (object$loglik[1] - object$loglik[2])
    rval$logtest <- c(test=logtest,
                      df=df,
                      pvalue= pchisq(logtest, df, lower.tail=FALSE))
    rval$sctest <- c(test=object$score,
                     df=df,
                     pvalue= pchisq(object$score, df, lower.tail=FALSE))
    rval$waldtest<-c(test=as.vector(round(object$wald.test, 2)),
                     df=df,
                     pvalue= pchisq(as.vector(object$wald.test), df,
                                    lower.tail=FALSE))
    if (!is.null(object$rscore))
      rval$robscore<-c(test=object$rscore,
                       df=df,
                       pvalue= pchisq(object$rscore, df, lower.tail=FALSE))
    rval$used.robust<-!is.null(object$naive.var)

    class(rval) <- "summary.IPWcprisk"
    rval
}

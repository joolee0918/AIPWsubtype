#' The heterogeneity test for the subtype-specific effet

#' Testing whether the exposures are not associated with the specific subtypes of disease.

#' @param fit an object of class \code{subtype}, \code{IPWsubtype}, or \code{AIPWsubtype}. If a single object \code{fit} was specified, it tests whether all unconstrained variables, so called, exposures are not associated with disease subtypes.
#' @param fit2 an object of the same class as \code{fit} to compare the heterogeneity effect.
#' @param data a data frame which was used in \code{fit}. It must be specified if a single object \code{fit} is specified.
#'
#' @details Log-likelihood tests are conducted to test whether the heterogeneity effects are significant. To test the heterogeneity effects, \code{marker_rr} can be used to drop out the marker-specific effects of exposures.

#' @return a data frame: names of variables which are compared, test statistics, the degree of freedom, and p-value.
#' @importFrom data.table data.table
#' @export
hetero.test <- function(fit, fit2 = NULL, data = NULL){

  mf <- list()

  Mloglik1 <- fit$loglik[2]
  beta <- fit$coefficients
  nabeta <- !(is.na(beta))  #non-missing coefs
  beta2 <- beta[nabeta]
  Mdf1 <- length(beta2)

  if(!is.null(fit2)){
    Mloglik2 <- fit2$loglik[2]
    beta <- fit2$coefficients
    nabeta <- !(is.na(beta))  #non-missing coefs
    beta2 <- beta[nabeta]
    Mdf2 <- length(beta2)

    if(Mdf1 > Mdf2) {
      logtest <- -2 * (Mloglik2 - Mloglik1)
      df <- Mdf1 - Mdf2
    }else {
      logtest <- -2 * (Mloglik1 - Mloglik2)
      df <- Mdf2 - Mdf1
    }
    pvalue = pchisq(logtest, df, lower.tail=FALSE)
    if(is.null(fit$subtype$marker_rr)){
      who1 <- fit$subtype$marker_name
    }else{
      who1 <- fit$subtype$marker_name[fit$subtype$marker_rr]
    }
    if(is.null(fit2$subtype$marker_rr)){
      who2 <- fit2$subtype$marker_name
    }else{
      who2 <- fit2$subtype$marker_name[fit2$subtype$marker_rr]
    }

    label <- paste(setdiff(who1, who2), collapse = " ")
    mf[[1]] <-  data.frame(label = label, test = logtest, df = df, pvalue = pvalue)
  } else{
    if(is.null(data)) stop("data should be included")
    Call <- fit$call
    tcall <- Call[c(1, match(c('formula', 'id', 'missing_model', 'missing_indep', 'two_stage', 'tstage_name', 'marker_name', 'second_cont_bl', 'constvar',
                             'init', 'control'), names(Call), nomatch=0))]

    tcall$data <- data
    tcall$first_cont_rr <- FALSE
    tcall$second_cont_rr <- FALSE
    tcall[[1L]] <- quote(AIPWsubtype)
    mm <- eval(tcall)
    if(is.null(fit$subtype$marker_rr)){
      who <- fit$subtype$marker_name
    }else{
      who <- fit$subtype$marker_name[fit$subtype$marker_rr]
    }
    logtest <- -2 * (mm$loglik[2] - Mloglik1)
    beta <- mm$coefficients
    nabeta <- !(is.na(beta))  #non-missing coefs
    beta2 <- beta[nabeta]
    df <-  Mdf1 - length(beta2)
    pvalue = pchisq(logtest, df, lower.tail=FALSE)
    label <- paste(who, collapse = " ")
    mf[[1]] <-  data.frame(label = label, test = logtest, df = df, pvalue = pvalue)

  }
  class(mf) <- c("hetero.test")

  return(mf)
}


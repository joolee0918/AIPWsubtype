

#' @importFrom data.table data.table
#' @export
loglik.test <- function(fit1, fit2 = NULL, data = NULL){

  Mloglik1 <- fit1$loglik[2]
  beta <- fit1$coefficients
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
    who1 <- fit1$subtype$marker_name[fit1$subtype$marker_rr]
    who2 <- fit2$subtype$marker_name[fit2$subtype$marker_rr]
    label <- paste(setdiff(who1, who2), collapse = " ")
    mf <-  data.frame(label = label, test = logtest, df = df, pvalue = pvalue)
  } else{
    if(!is.null(data)) stop("data should be included")

    tcall <- Call[c(1, match(c('formula', 'id', 'missing_model', 'missing_indep', 'two_stage', 'tstage_name', 'marker_name', 'second_cont_bl', 'constvar',
                             'init', 'control'), names(Call), nomatch=0))]

    tcall$data <- data
    tcall$first_cont_rr <- FALSE
    tcall$second_cont_rr <- FALSE
    tcall[[1L]] <- quote(AIPWsubtype)
    mm <- eval(tcall)
    who <-  paste(marker_name, collapse = "+")
    logtest <- -2 * (mm$loglik[2] - Mloglik1)
    beta <- mm$coefficients
    nabeta <- !(is.na(beta))  #non-missing coefs
    beta2 <- beta[nabeta]
    df <-  Mdf1 - length(beta2)
    pvalue = pchisq(logtest, df, lower.tail=FALSE)

    mf <-  data.frame(label = who, test = logtest, df = df, pvalue = pvalue)

  }
  class(mf) <- "hetero.test"

  return(mf)
}


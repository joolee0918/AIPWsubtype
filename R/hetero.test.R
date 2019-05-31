

#' @importFrom data.table data.table
#' @export
hetero.test <- function(object, data){

  Call <- object$call

  Mloglik <- object$loglik[1]
  beta <- object$coefficients
  nabeta <- !(is.na(beta))  #non-missing coefs
  beta2 <- beta[nabeta]
  Mdf <- length(beta2)

  nullmodel <- coxph(object$formula, data, method="breslow")

  marker_name <- object$subtype$marker_name
  missing_model <- eval(Call[c(match(c('missing_model'), names(Call), nomatch=0))][[1]])
  n_marker <- length(marker_name)
  marker_rr <- object$subtype$marker_rr
  if(is.null(marker_rr)) marker_rr <- rep(TRUE, n_marker)

  mf <- vector('list', n_marker)
  for(k in 1:(n_marker-1)){
    comb <- combn(n_marker, k)
    J <- ncol(comb)
    who <- logtest <- df <- pvalue <- rep(0, J)
      for(i in 1:J){
        tcall <- Call[c(1, match(c('formula', 'id', 'missing_model', 'missing_indep', 'two_stage', 'tstage_name', 'marker_name', 'second_cont_bl', 'first_cont_rr', 'second_cont_rr', 'constvar',
                             'init', 'control'), names(Call), nomatch=0))]

        tmp <- marker_rr
        tmp[c(comb[, i])] <- FALSE
        tcall$marker_rr <- tmp
        #print(tcall$marker_name)
        tcall$data <- data
        tcall[[1L]] <- quote(AIPWsubtype)
        mm <- eval(tcall)
        print(summary(mm))
        who[i] <- paste(marker_name[c(comb[, i])], collapse=" ")
        logtest[i] <- -2 * (Mloglik - mm$loglik[2])
        beta <- mm$coefficients
        nabeta <- !(is.na(beta))  #non-missing coefs
        beta2 <- beta[nabeta]
        df[i] <-  Mdf - length(beta2)
        pvalue[i] = pchisq(logtest[i], df[i], lower.tail=FALSE)
      }
    mf[[k]] <- data.frame(label = who, test = logtest, df = df, pvalue = pvalue)
  }


  tcall <- Call[c(1, match(c('formula', 'id', 'missing_model', 'missing_indep', 'two_stage', 'tstage_name', 'marker_name', 'second_cont_bl', 'second_cont_rr', 'constvar',
                             'init', 'control'), names(Call), nomatch=0))]

  tcall$data <- data
  tcall$first_cont_rr <- FALSE
  #print(tcall$marker_name)
  tcall[[1L]] <- quote(AIPWsubtype)
  mm <- eval(tcall)
  who <-  paste(marker_name, collapse = "+")
  logtest <- -2 * (Mloglik - mm$loglik[2])
  beta <- mm$coefficients
  nabeta <- !(is.na(beta))  #non-missing coefs
  beta2 <- beta[nabeta]
  df <-  Mdf - length(beta2)
  pvalue = pchisq(logtest, df, lower.tail=FALSE)

  mf[[n_marker]] <-  data.frame(label = who, test = logtest, df = df, pvalue = pvalue)

  class(mf) <- "hetero.test"

  return(mf)
}


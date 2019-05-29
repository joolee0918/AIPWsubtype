## Modify getHaz in eha packages

cumHaz <- function(Y, strats, score, weights = NULL){
  if (NCOL(Y) == 2) {
    Y <- cbind(rep(0, NROW(Y)), Y)
   }else{
    if (NCOL(Y) != 3) stop("'Y' is of wrong type.")
    }

  if (missing(score)) score <- rep(1, NROW(Y))
  if (missing(strats)) strats <- rep(1, NROW(Y))
  strats <- as.factor(strats)
  Strata <- levels(strats)
  ns <- length(Strata)

  out <- vector("list", ns)
  for (j in seq_along(Strata)){
    enter <- Y[strats == Strata[j], 1]
    exit <- Y[strats == Strata[j], 2]
    event <- Y[strats == Strata[j], 3] != 0
    sco <- score[strats == Strata[j]]
    time <- sort(unique(exit[event]))
    haz <- matrix(0, ncol = 2, nrow = length(time))
    haz[, 1] <- time
    for (i in seq_along(time)){
      rs <- (enter < time[i] & exit >= time[i])
      if(is.null(weights)) {
        haz[i, 2] <- sum(event[exit == time[i]]) /sum(sco[rs])
      } else{
        haz[i, 2] <- sum(weights[event==time[i]])*event[exit == time[i]]) /sum(sco[rs])
      }
    }
    out[[j]] <- haz
  }
  names(out) <- Strata
  class(out) <- "hazdata"
  out
}

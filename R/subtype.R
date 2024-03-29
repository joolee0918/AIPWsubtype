#' Cox proportional hazard model for competing risks data
#'
#' Fitting a Cox proportional hazard model for competing risks data, in which marker variables define the subyptes of outcome.
#'
#' @param formula a formula object with an obect of the type \code{Surv} on the left side and the terms on the right side. It should not include marker variables.
#' @param data a data.frame which has variables in formula and markers.
#' @param id a charhacter string specifying subject IDs.
#' @param marker_name a vector of charhacter strings specifying markers defining cause of failures.
#' @param marker_rr a vector of logical value. Dafault is NULL, in which a model includes all markers named in \code{marker_name} in modeling heterogeneity effects. Otherwise, a vector of logical value can spcify whether each marker's heterogeneity effect will be examined or not. A length of this should be equal to that of \code{marker_name}.
#' @param first_cont_rr a logical value: if \code{TRUE}, the first order contrasts are included in modeling cause-specific relative risks based on log-linear representation. Otherwise the first contrasts are only included.
#' @param second_cont_bl a logical value: if \code{TRUE}, the second order contrasts are included in modeling cause-specific baseline functions based on log-linear representation. Otherwise the first contrasts are only included.
#' @param second_cont_rr a logical value: if \code{TRUE}, the second order contrasts are included in modeling cause-specific relative risks based on log-linear representation. Otherwise the first contrasts are only included.
#' @param constvar a vector of character strings specifying constrained varaibles of which the effects on the outcome are to be the same across subtypes of outcome. The variables which are not specified in \code{constvar} are considered as unconstrained variabales of which the associations with the outcome may be different across the outcome subtypes.
#' @param init a vector of initial values of the iteration. Default value is zero for all variables.
#' @param control an object of class \code{\link[survival]{coxph.control}} in \code{survival} packages.
#'
#' @details The Cox proportional hazard model is used to model cause-specific hazard functions. To examine the association between exposures and the specific subtypes of disease, the log-linear is used for reparameterization.
#' This is a wrapper function for \code{\link[survival]{coxph}} after data duplication. The returned value \code{x}, \code{y} are duplicated the number of subtypes of outcome. \code{+cluster()} should be avoided since we automatically include it in the formula.
#'
#' For marker variables, 0 indicates censored events.
#'
#' @return A returned object is an object of class \code{\link[survival]{coxph}}. See \code{\link[survival]{coxph.object}} for details. The additional returned values include the following:
#' \item{basehaz}{estimated baseline cause-specific hazard functions the reference disease subtype corresponding to marker variables equal to 1.}
#' \item{subtype}{a list of values related to subtypes including the number of subtypes, character strings of marker names, etc.}
#'

#' @examples
#' m1 <- subtype(Surv(start, time, status)~ X + W,  data = data, id = "id", marker_name = c("y1", "y2"), second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = "W")
#'

#' @export
subtype <- function(formula, data, id,  marker_name, marker_rr = NULL,
                    first_cont_rr = TRUE, second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = NULL, init, control,  x = FALSE, y = TRUE, model = FALSE) {


  Call <- match.call()
  rx <- x
  rmodel <- model

  if (missing(control)) {
    control = coxph.control()
     }


  data <- data[order(data[, id]), ]
  data$rowid <- seq(1, nrow(data))

  n <- nrow(data)
  n_marker <- length(marker_name)
  marker <- data.frame(data[, marker_name])
  names(marker) <- marker_name
  marker[marker == 0] <- NA

  n_subtype = 1
  for (i in 1:n_marker) {
    n_subtype = n_subtype * (nlevels(factor(marker[, i])))
  }


  # possible marker combination

  level_y = list()
  for (k in 1:n_marker) {
    level_y[[k]] <- seq(nlevels(factor(marker[, k])), 1)
  }
  tmpy <- list()
  for (k in 1:n_marker) {
    tmpy[[k]] <- as.vector(level_y[[k]])
  }
  total_subtype <- as.data.frame(expand.grid(tmpy))
  names(total_subtype) <- marker_name

  ## Data frame for cause

  umarker <- unique(na.omit(marker))
  tmpmar <- rep(0, nrow(umarker))
  for(i in 1:nrow(umarker)){
    for (j in 1:n_subtype) {
      if (all(umarker[i, ] == total_subtype[j, 1:n_marker]))
        tmpmar[i] <- j
    }
  }

  ototal_subtype <- total_subtype[sort(tmpmar),]
  on_subtype <- nrow(ototal_subtype)

  event <- tail(terms.inner(formula[1:2]), 1)

  cause <- rep(NA, n)
  R <- rep(1, n)
  cause <- findcause(R, cause, data[, event], as.matrix(marker), on_subtype, as.matrix(ototal_subtype))

  # observed cause
  ocause <- sort(unique(cause))

  for (i in 1:n_marker) {
    marker[, i] <- factor(marker[, i])
  }
  marker <- as.data.frame(marker)
  data[, marker_name] <- marker


  lf <- function(x) {
    if (!is.na(x)) {
      res <- c(x, ocause[!(ocause %in% x)])
    } else {
      res <- ocause
    }
    res
  }

  newcause <- unlist(lapply(1:n, function(i) lf(cause[i])))
  newdata <- data[rep(1:n, each = on_subtype), ]
  newdata[, event] <- rep(0, n*on_subtype)
  newdata[seq(1, n*on_subtype, by=on_subtype), event] <- data[, event]
  newdata[, marker_name] <- data.frame(ototal_subtype[newcause, ])
  term_marker <- rep(0, n_marker)

  for (i in 1:n_marker) term_marker[i] <- paste("factor(", marker_name[i], ")", sep = "")
  if(is.null(marker_rr)) {
    term_marker_rr <- term_marker
  }else{
    term_marker_rr <- term_marker[marker_rr]
  }


  special <- c("strata", "cluster")
  Tf <- terms(formula, special)
  strats <- attr(Tf, "specials")$strata
  if (length(strats)) {
    strterm <- survival:::untangle.specials(Tf, "strata", 1)$terms
    Xattr <- attr(Tf, "term.labels")[-strterm]
  } else {
    Xattr <- attr(Tf, "term.labels")
  }

  if(is.null(constvar)){
    unconstvar <- Xattr
  }else{
    unconstvar <- Xattr[-c(grep(paste(constvar, collapse="|"), Xattr))] # Xattr[!(Xattr %in% constvar)]
    constvar <- Xattr[!(Xattr %in% unconstvar)]
  }


  nuvar <- length(unconstvar)
  order_rr <- NULL
  order_bl <- NULL
  if(first_cont_rr == TRUE){
    if (nuvar == 0) {
      order_rr <- paste(term_marker_rr, collapse = "+")
    } else {
      for (i in 1:nuvar) {
        tmp_rr <- paste(term_marker_rr, "*", unconstvar[i], collapse = "+")
        order_rr <- paste(order_rr, tmp_rr, sep = "+")
      }
    }
  }

  order_bl <- paste(term_marker, collapse= "+")



  if(n_marker > 1) pairm <- combn(n_marker, 2)

  if (second_cont_rr == TRUE) {
    order_rr <- NULL
    for (i in 1:nuvar) {
      tmp_rr <- paste(attr(terms(as.formula(paste("~", unconstvar[i], "*", "(", paste(term_marker_rr, collapse = "+"),
                                                  ")", "^", 2))), "term.labels"), collapse = "+")
      order_rr <- paste(order_rr, tmp_rr, sep = "+")
    }
  }


  if (second_cont_bl == TRUE & second_cont_rr == FALSE) {
    for (i in 1:ncol(pairm)) {
      tmp_bl <- paste(term_marker[pairm[1, i]], term_marker[pairm[2, i]], sep = ":")
      order_bl <- paste(order_bl, tmp_bl, sep = "+")
    }
  }

  drop_bl <- NULL
  if (second_cont_bl == FALSE & second_cont_rr == TRUE) {
    drop_bl <- rep(0, ncol(pairm))
    for (i in 1:ncol(pairm)) {
      drop_bl[i] <- paste(term_marker[pairm[1, i]], term_marker[pairm[2, i]], sep = ":")
    }
    order_bl <- paste("-", paste(drop_bl, collapse = "-"))
  }


  newformula <- update.formula(formula, paste("~.+", order_bl, order_rr, "+", "cluster", "(", id, ")"))
  fit <- coxph(formula = newformula, data = newdata, control = control, robust = T, method = "breslow", model = rmodel, x = TRUE)

  if(is.null(fit$strata)) {
    stratum <- rep(1, nrow(fit$y))
  }else {
    stratum <- fit$strata
  }
  lp <- fit$linear.predictors
  s0 <- exp(lp)
  basehaz <- baseHaz(fit$y, stratum, s0)

  fit$n <- n
  fit$subtype = list(n_subtype = on_subtype, marker_name = marker_name, total_subtype = ototal_subtype, marker_rr = marker_rr)
  fit$basehaz <- basehaz
  fit$call <- Call

  if(!rx) fit$x <- NULL

  return(fit)

}



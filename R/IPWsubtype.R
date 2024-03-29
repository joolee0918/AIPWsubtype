#' Inverse probability weighted Cox proportional hazard model for competing risks data
#'
#' Fitting an inverse probability weighted Cox proportional hazard model for competing risks data with partially missing markers, in which marker variables define the subyptes of outcome.
#'
#' @param formula a formula object with an obect of the type \code{Surv} on the left side and the terms on the right side. It should not include marker variables.
#' @param data a data.frame which has variables in formula and markers.
#' @param id a charhacter string specifying subject IDs.
#' @param missing_model a character string specifying the approach of the missingness model.\code{missing_model = "condi"} represents the conditional approach and \code{missing_model = "multinom"} uses a multinomial model.
#' @param missing_formula a list of the missingness model formula for each marker. A right side of formula object including a ~ operator. If \code{two_stage} is \code{TRUE}, a model for the first-stage missingness model should be included at the last element of a list.
#' @param missing_indep a logical value: if \code{TRUE}, markers are assumed to be independent.
#' @param two_stage  a logical value: if \code{TRUE}, a two-stage missingness model is used. Default is \code{FALSE}.
#' @param tstage_name a charhacter string specifying the first-stage missingness variable if \code{two_stage} is \code{TRUE}.
#' @param marker_name a vector of charhacter strings specifying markers defining cause of failures.
#' @param marker_rr a vector of logical value. Dafault is NULL, in which a model includes all markers named in \code{marker_name} in modeling heterogeneity effects. Otherwise, a vector of logical value can spcify whether each marker's heterogeneity effect will be examined or not. A length of this should be equal to that of \code{marker_name}.
#' @param first_cont_rr a logical value: if \code{TRUE}, the first order contrasts are included in modeling cause-specific relative risks based on log-linear representation. Otherwise the first contrasts are only included.
#' @param second_cont_bl a logical value: if \code{TRUE}, the second order contrasts are included in modeling cause-specific baseline functions based on log-linear representation. Otherwise the first contrasts are only included.
#' @param second_cont_rr a logical value: if \code{TRUE}, the second order contrasts are included in modeling cause-specific relative risks based on log-linear representation. Otherwise the first contrasts are only included.
#' @param constvar a vector of character strings specifying constrained varaibles of which the effects on the outcome are to be the same across subtypes of outcome. The variables which are not specified in \code{constvar} are considered as unconstrained variabales of which the associations with the outcome may be different across the outcome subtypes.
#' @param init a vector of initial values of the iteration. Default value is zero for all variables.
#' @param control an object of class \code{\link[survival]{coxph.control}} in \code{survival} packages. The default value of \code{iter.max} is 2000 and that of \code{eps} is 1e-12. See \code{\link[survival]{coxph.control}} for other values.

#' @details The Cox proportional hazard model with weights is used to model cause-specific hazard functions. To examine the association between exposures and the specific subtypes of disease, the log-linear is used for reparameterization. The weights for the complete cases are obtained by fitting logistic regression models \code{\link[stats]{glm}}.
#' The data duplication method is used so that the returned value \code{x}, \code{y}, and \code{weights} are duplicated the number of subtypes of outcome. Special terms including \code{+strata()} and \code{+offset()} can be used. \code{+cluster()} should be avoided since we automatically include it in the formula. Breslow method is used for handling tied event times.
#'
#' For marker variables, 0 indicates censored events or missing values.
#'
#' @return An object of class \code{IPWsubtype} representing the fit.
#' \item{coefficients}{the estimated regressoin coefficients.}
#' \item{naive.var}{the inverse of estimated Information matrix.}
#' \item{var}{the robust sandwich variance estimate.}
#' \item{linear.predictors}{a vector of linear predictors. This is not centered.}
#' \item{loglik}{a vector of length 2 containing the log-likelihood with the initial values and with the estimated coefficients.}
#' \item{score}{a value of the score test at the initial value of the coefficients.}
#' \item{rscore}{a value of the robust log-rank statistic.}
#' \item{score.residual}{the score residuals for each subject. For incomplece cases, the value is 0.}
#' \item{iter}{a number of iterations used.}
#' \item{weights}{a vector of weights used, which is the inverse of the probability of complete case given the event occurs.}
#' \item{basehaz}{estimated baseline cause-specific hazard functions the reference disease subtype corresponding to marker variables equal to 1.}
#' \item{Ithealp}{a matrix of the partial derivative of the score functions with respect to the parameters from the missingness models.}
#' \item{model.missing}{a list of an object of class \code{glm} fitting the missingness models.}
#' \item{n}{the number of observations.}
#' \item{nc}{the number of complete-case observations.}
#' \item{nevent}{the number of events.}
#' \item{ncevent}{the number of complete-case events.}
#' \item{subtype}{a list of values related to subtypes including the number of subtypes, character strings of marker names, etc.}

#' The object will also contain the following: strata, formula, call, terms, y, offset, xlevels, optionally x, and model.

#' @importFrom stats model.frame model.extract model.matrix model.offset aggregate update.formula printCoefmat rexp qlnorm rnorm pchisq pnorm predict glm
#' @import survival
#' @importFrom Matrix bdiag
#' @importFrom data.table rbindlist
#' @importFrom sandwich estfun
#'
#' @examples
#' m1 <- IPWsubtype(Surv(start, time, status)~ X + W,  data = data, id = "id", missing_formula = list(~time + X, ~time + X),
#'  two_stage = FALSE, marker_name = c("y1", "y2"), second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = "W")
#'
#' # Two-stage missingness model
#' m2 <- IPWsubtype(Surv(start, time, status)~ X + W,  data = data, id = "id", missing_formula = list(~time + X, ~time + X, ~time + X + W),
#'  two_stage = TRUE, tstage_name = c("R0"), marker_name = c("y1", "y2"), second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = "W")

#' @export
IPWsubtype <- function(formula, data, id, missing_model = c("condi", "multinom"), missing_formula, missing_indep = FALSE, two_stage = FALSE, tstage_name = NULL,  marker_name,
                       marker_rr = NULL, first_cont_rr = TRUE, second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = NULL, init, control,  x = FALSE, y = TRUE, model = FALSE) {


    Call <- match.call()

    if(missing(id)) stop("id must be specified")
    if(missing(marker_name)) stop("marker_name must be specified")
    if(missing(missing_model)) stop("missing_model must be specified")
    if(missing(missing_formula)) stop("missing_formula must be specified")

    if (missing(control)) {
      control = coxph.control()
      control$eps = 1e-09
      control$iter.max = 2000
      control$toler.inf = sqrt(control$eps)
    }

    rx <- x
    ry <- y
    rmodel <- model

    data <- data[order(data[, id]), ]
    data$rowid <- seq(1, nrow(data))

    n <- nrow(data)
    uniqid <- unique(data[, id])
    n_marker <- length(marker_name)
    marker <- data.frame(data[, marker_name])
    names(marker) <- marker_name
    marker[marker == 0] <- NA

    n_subtype = 1
    for (i in 1:n_marker) {
      n_subtype = n_subtype * (nlevels(factor(marker[, i])))
    }


    ## One-stage / Two-stage -> number of possible case of R

    if (two_stage == T) {
      nR = 2^(n_marker) + 1
      nvecR = n_marker + 1
    } else {
      nR = 2^(n_marker)
      nvecR = n_marker
    }

    if (two_stage == T) {
      tmpr = rep(list(c(1, 0)), n_marker)
      total_R <- as.matrix(expand.grid(tmpr))
      total_R <- cbind(total_R, rep(1, 2^n_marker))
      total_R <- rbind(total_R, rep(0, n_marker + 1))
      total_R <- cbind(total_R, seq(1, 2^(n_marker) + 1))
      colnames(total_R) <- c(paste("R", c(1:n_marker), sep = ""), "R0", "R")
    } else {
      tmpr <- rep(list(c(1, 0)), n_marker)
      total_R <- as.matrix(expand.grid(tmpr))
      total_R <- cbind(total_R, seq(1, 2^(n_marker)))
      colnames(total_R) <- c(paste("R", c(1:n_marker), sep = ""), "R")
    }
    total_R <- as.matrix(total_R)


    ### Marker

    event <- tail(terms.inner(formula[1:2]), 1)
    Rmat <- matrix(1, nrow = n, ncol = n_marker)
    for (i in 1:n_marker) {
      Rmat[, i] <- ifelse(data[, event] == 1, 1 - as.integer(is.na(marker[, i])), 1)
    }

    uRmat <- unique(Rmat)
    tmpR <- rep(0, nrow(uRmat))
    for(i in 1:nrow(uRmat)){
      for (j in 1:2^n_marker) {
        if (all(uRmat[i, ] == total_R[j, 1:n_marker]))
          tmpR[i] <- j
      }
    }
    ototal_R <- total_R[sort(tmpR),]
    if(two_stage == TRUE) ototal_R <- rbind(ototal_R, total_R[nR, ])


    onR <- nrow(ototal_R)
    fonR <- onR-1*(two_stage==T)

    R = rep(1, n)
    if (two_stage == T)
      R[data[, event]==1 & data[, tstage_name] == 0] <- onR
    R <- findR(R, data[, event], fonR, Rmat, ototal_R[, 1:n_marker])

    # observed R
    oR <- sort(unique(R))

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

    ccdata <- data[R == 1, ]
    ndata <- nrow(ccdata)
    missing_prob <- rep(1, ndata)
    newedata <- ccdata[ccdata[, event] == 1, ]
    nnevent <- nrow(newedata)
    nuniqid <- unique(ccdata[, id])
    newmarker <- as.data.frame(marker[R == 1, ])
    names(newmarker) <- marker_name

    cause <- rep(NA, ndata)
    cause <- findcause(rep(1,ndata), cause, ccdata[, event], as.matrix(newmarker), on_subtype, as.matrix(ototal_subtype))

    # observed cause
    ocause <- sort(unique(cause))


    for (i in 1:n_marker) {
      marker[, i] <- factor(marker[, i])
    }
    marker <- as.data.frame(marker)
    data[, marker_name] <- marker



    ## Data frame for missing indicator

    # model missing

    missing_data <- as.data.frame(Rmat)
    names(missing_data) <- c(paste("R", c(1:n_marker), sep = ""))
    if (two_stage == TRUE) {
      missing_data$R0 <- data[, tstage_name]
    }
    missing_data$R <- R
    missing_name <- names(missing_data)

    data$R <- R
    edata <- data[data[, event] == 1, ]
    eventid <- edata[, id]
    nevent <- nrow(edata)

    model_missing <- list()

    if(missing_model == "multinom"){
      if(two_stage == FALSE){
        tmpedata <- mlogit.data(edata, shape="wide", choice="R")
        mult_formula <- as.formula(paste("R", "~1|", paste(missing_formula[[1]][2], collapse = "")))
        model_missing[[1]] <- mlogit(mult_formula, data=tmpedata, Hess=TRUE)

      } else{

        tmpedata <- mlogit.data(edata[edata[, tstage_name] == 1, ], shape="wide", choice="R")
        mult_formula <- as.formula(paste("R", "~1|", paste(missing_formula[[1]][2], collapse = "")))
        model_missing[[1]] <- mlogit(mult_formula, data=tmpedata, Hess=TRUE)

        tmp_model <- as.formula(paste(tstage_name, paste(missing_formula[[2]], collapse = "")))
        model_missing[[2]] <- glm(tmp_model, data = edata, family = binomial)
      }
    } else{

      edata <- cbind(edata, missing_data[data[, event]==1,])

      model_missing <- list()
      tmp_model <- as.formula(paste(missing_name[1], paste(missing_formula[[1]], collapse = "")))

      if (two_stage == FALSE) {

        for (i in 1:n_marker) {
          model_missing[[i]] <- glm(tmp_model, data = edata, family = binomial)
          if (i == n_marker)
            break
          tmp_mname <- missing_name[1:i]
          if (missing_indep == FALSE)
            missing_formula[[i + 1]] <- paste(paste(missing_formula[[i + 1]], collapse = ""), paste(tmp_mname,
                                                                                                    collapse = "+"), sep = "+")
          tmp_model <- as.formula(paste(missing_name[i + 1], paste(missing_formula[[i + 1]], collapse = "")))
        }
      } else {
        for (i in 1:n_marker) {
          model_missing[[i]] <- glm(tmp_model, data = edata[edata[, tstage_name] == 1, ], family = binomial)
          if (i == n_marker)
            break
          tmp_mname <- c(missing_name[1:i])
          if (missing_indep == FALSE)
            missing_formula[[i + 1]] <- paste(paste(missing_formula[[i + 1]], collapse = ""), paste(tmp_mname,
                                                                                                    collapse = "+"), sep = "+")
          tmp_model <- as.formula(paste(missing_name[i + 1], paste(missing_formula[[i + 1]], collapse = "")))
        }

        tmp_model <- as.formula(paste(tstage_name, paste(missing_formula[[nvecR]], collapse = "")))

        model_missing[[nvecR]] <- glm(tmp_model, data = edata, family = binomial)

      }

    }




    if(missing_model == "multinom"){

      
      nalp0 = 0
      if (two_stage == T)
        nalp0 = length(model_missing[[2]]$coefficients)

      p1 <- rep(1, ndata)
      if(two_stage == FALSE){

        tmpdata <- mlogit.data(edata, shape="wide", choice="R")

        whereR1 <- (edata$R == 1)
        tmpp1 <- model_missing[[1]]$probabilities[whereR1 ,1]

        p1[ccdata[, event] == 1] <- tmpp1
        nalp <- length(coef(model_missing[[1]]))
      }else {

        edata$RR <- ifelse(edata$R == max(oR), 1, edata$R)
        tmpdata <- mlogit.data(edata[edata[, tstage_name]==1,], shape="wide", choice="RR")

        whereR1 <- (edata[edata[, tstage_name] == 1,]$R == 1)
        tmpp1 <- model_missing[[1]]$probabilities[whereR1 ,1]
        p0 <-  as.vector(predict(model_missing[[2]], newdata = newedata, type = "response"))

        p1[ccdata[, event] == 1] <- tmpp1*p0
        nalp <- length(coef(model_missing[[1]])) + length(model_missing[[2]]$coefficients)

      }

      dpR <- as.data.frame(matrix(0, ncol = nalp, nrow = ndata))

      tmpf <- mFormula(mult_formula)
      newX <- model.matrix(tmpf, tmpdata)
      tmpX <-  model.matrix(missing_formula[[1]], newedata)
      scoreX <- matrix(exp(as.vector(newX%*%coef(model_missing[[1]]))), byrow=T, ncol=fonR)[whereR1,-1]

      dpR[ccdata[, event] == 1,] <- dpR1_multinom(tmpX, scoreX, onR, two_stage, ncol(scoreX), nalp, nalp0)

      if (two_stage == T) {
        tmp <- model.matrix(model_missing[[2]]$formula, newedata) *p0/(1 + exp(predict(model_missing[[2]], newdata = newedata)))
        dpR[ccdata[, event] == 1, 1:(nalp - nalp0)] <- dpR[ccdata[, event] == 1, 1:(nalp - nalp0)]*p0
        dpR[ccdata[, event] == 1, (nalp - nalp0 + 1):(nalp)] <- tmp*tmpp1
      }


    } else{
      newedata <- cbind(newedata, missing_data[data[, event]==1 & R==1,])

      nalp = 0
      for (k in 1:nvecR) {
        nalp = nalp + length(model_missing[[k]]$coefficients)
      }
      nalp0 = 0
      if (two_stage == T)
        nalp0 = length(model_missing[[nvecR]]$coefficients)


      PI <- matrix(1, nrow = ndata, ncol = nvecR)
      est_pi <- matrix(0, nrow = nnevent, ncol = nvecR)

      for (k in 1:nvecR) {
        tmp <- predict(model_missing[[k]], newdata = newedata, type = "response")
        est_pi[, k] <- tmp
        PI[ccdata[, event] == 1, k] <- tmp
      }
      p1 <- apply(PI, 1, prod)

      dpR <- as.data.frame(matrix(0, ncol = nalp, nrow = ndata))

      a = 1
      for (k in 1:nvecR) {
        b <- length(model_missing[[k]]$coefficients)
        tmp <- model.matrix(model_missing[[k]]$formula, newedata) * predict(model_missing[[k]], newdata = newedata,
                                                                            type = "response")/(1 + exp(predict(model_missing[[k]], newdata = newedata))) * apply(as.matrix(est_pi[,-k]), 1, prod)
        dpR[ccdata[, event] == 1, a:(a + b - 1)] <- tmp
        a = a + b
      }


    }

    ## Data duplication
    lf <- function(x) {
      if (!is.na(x)) {
        res <- c(x, ocause[!(ocause %in% x)])
      } else {
        res <- ocause
      }
      res
    }


    newcause <- unlist(lapply(1:ndata, function(i) lf(cause[i])))
    ccdata$pi1 <- p1
    newdata <- ccdata[rep(1:ndata, each = on_subtype), ]
    newdata[, event] <- rep(0, ndata*on_subtype)
    newdata[seq(1, ndata*on_subtype, by=on_subtype), event] <- ccdata[, event]
    newdata[, marker_name] <- data.frame(ototal_subtype[newcause, ])
    dpR <- dpR[rep(seq_len(nrow(dpR)), each = on_subtype), ]

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


    #whichX <- which(dimnames(attr(Tf, "factors"))[[2]] %in% c(unconstvar))
    #whereX <- which(Xattr$assign == whichX)
    # where is the main effect for uonstrained variables in X matrix
    #whichW <- which(dimnames(attr(Tf, "factors"))[[2]] %in% c(constvar))
    #whereW <- which(Xattr$assign == whichW)

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


    fit <- coxph(formula = newformula, data = newdata, weights = 1/pi1, control = control, robust = T, model = TRUE, x = TRUE, method = "breslow", init=init)

    if (!is.null(fit$fail)) {
        fit <- list(fail = fit)
        class(fit) <- "IPWcprisk"
    } else {
        mf <- model.frame(newformula, data = newdata)

        eid <- newdata[newdata[, event] == 1, id][1]
        ss <- newdata[, id] %in% eid
        y <- fit$y
        ny <- ncol(y)
        if (ny == 2) {
            time <- as.double(y[, 1])
            status <- as.integer(y[, 2])
        } else {
            start <- as.double(y[, 1])
            stop <- as.double(y[, 2])
            status <- as.integer(y[, 3])
        }

        x <- fit$x
        weights <- fit$weights
        strata <- fit$strata
        lp <- fit$linear.predictors

        nused <- nrow(y)
        nvar <- ncol(x)


        if (ny == 2) {
            if (length(strata) == 0) {
                ord <- order(time)
                strata <- NULL
                newstrat <- as.integer(rep(0, nused))
            } else {
                ord <- order(strata, time)
                newstrat <- as.integer(c(1 * (diff(as.numeric(strata[ord])) != 0), 1))
            }
        } else {
            if (length(strata) == 0) {
                sort.end <- as.integer(order(-y[, 2])) - 1L  #indices start at 0 for C code
                sort.start <- as.integer(order(-y[, 1])) - 1L
                newstrat <- nused
            } else {
                sort.end <- as.integer(order(strata, -y[, 2])) - 1L
                sort.start <- as.integer(order(strata, -y[, 1])) - 1L
                newstrat <- cumsum(table(strata))
            }
        }


        if (ny == 2) {
            time2 <- time[ord]
            status2 <- status[ord]
            x2 <- x[ord, ]
            status2 <- status[ord]
            lp2 <- lp[ord]
            weights2 <- weights[ord]
            dpR2 <- as.matrix(dpR[ord, ])

            Ithealp <- IPW_ithealp_cox(time2, status2, lp2, newstrat, x2, dpR2, weights2, nused, nvar, nalp)
        } else {

            Ithealp <- IPW_ithealp_ag(start, stop, status, lp, newstrat, sort.start, sort.end, x, as.matrix(dpR),
                weights, nused, nvar, nalp)
        }


        score_theta <- as.data.frame(residuals(fit, type = "score", collapse = newdata[, id], weighted = T))
        Stheta <- as.data.frame(uniqid)
        names(Stheta) <- id
        score_theta <- cbind(nuniqid, score_theta)
        colnames(score_theta)[1] <- id
        Stheta <- merge(Stheta, score_theta, by = id, all = T)
        Stheta[is.na(Stheta)] <- 0
        Stheta <- as.matrix(Stheta[, -1])

        fit2 <- fit
        fit2$linear.predictors <- 0*fit$linear.predictors
        score_theta0 <- as.data.frame(residuals(fit2, type = "score", collapse = newdata[, id], weighted = T))
        Stheta0 <- as.data.frame(uniqid)
        names(Stheta0) <- id
        score_theta0 <- cbind(nuniqid, score_theta0)
        colnames(score_theta0)[1] <- id
        Stheta0 <- merge(Stheta0, score_theta0, by = id, all = T)
        Stheta0[is.na(Stheta0)] <- 0
        Stheta0 <- as.matrix(Stheta0[, -1])


        # P I_thealp = t(tmp_score)%*%as.matrix(dpR)

        ## Only fully observed data

        Itheta <- fit$naive.var
        RVar_theta <- fit$var

        # For all data

        Salp = as.data.frame(uniqid)
        names(Salp) <- id

        if(missing_model == "multinom"){
          Ialp <- as.matrix(vcov(model_missing[[1]]))
          class(Ialp) <- "matrix"
          Ualp <- estfun(model_missing[[1]])

          if(two_stage == TRUE){
            Ialp <- bdiag(Ialp, vcov(model_missing[[2]]))
            Ialp <- as.matrix(Ialp)
            Ualp_ts <- estfun(model_missing[[2]])
          }
        } else{

          Ialp <- vcov(model_missing[[1]])
          if (nvecR > 1) {
            for (k in 2:nvecR) {
              Ialp <- bdiag(Ialp, vcov(model_missing[[k]]))
            }
          }
          Ialp <- as.matrix(Ialp)

          Ualp <- list()

          for (k in 1:n_marker) {
            Ualp[[k]] <- estfun(model_missing[[k]])
          }
          if (two_stage == TRUE) {
            Ualp_ts <- estfun(model_missing[[nvecR]])
          }
          Ualp <- do.call(cbind, Ualp)

        }

        if (two_stage == FALSE) {
          Ualp <- cbind(eventid, Ualp)
          colnames(Ualp)[1] <- id
        } else {
          Ualp <- cbind(edata[edata[, tstage_name] == 1, id], Ualp)
          colnames(Ualp)[1] <- id
          Ualp_ts <- cbind(eventid, Ualp_ts)
          colnames(Ualp_ts)[1] <- id
          Ualp <- suppressWarnings(merge(Ualp, Ualp_ts, by = id, all = T))
          Ualp[is.na(Ualp)] <- 0
        }

        Salp <- suppressWarnings(merge(Salp, Ualp, by = id, all = T))
        Salp[is.na(Salp)] <- 0
        Salp <- as.matrix(Salp[, -1])

        resid <- as.matrix(Stheta) - as.matrix(Salp) %*% Ialp %*% t(Ithealp)
        var <- fit$naive.var %*% t(resid) %*% resid %*% fit$naive.var

        ## robust log-rank statistic
        temp0 <-  as.matrix(Stheta0) - as.matrix(Salp) %*% Ialp %*% t(Ithealp)
        u <- apply(as.matrix(temp0), 2, sum)
        rscore <- survival::coxph.wtest(t(temp0)%*%temp0, u, control$toler.chol)$test

        #Wald test
        if (length(fit$coefficients)) {
          #not for intercept only models, or if test is already done
          nabeta <- !is.na(fit$coefficients)
          # The init vector might be longer than the betas, for a sparse term
          if (missing(init)) temp <- fit$coefficients[nabeta]
          else temp <- (fit$coefficients -
                          init[1:length(fit$coefficients)])[nabeta]
          wald.test <-  survival::coxph.wtest(var[nabeta,nabeta], temp,
                                              control$toler.chol)$test
        }

        if(is.null(fit$strata)) {
          stratum <- rep(1, nrow(fit$y))
        }else {
          stratum <- fit$strata
        }
        s0 <- exp(lp)*weights
        basehaz <- baseHaz(fit$y, stratum, s0, weights = weights)

        afit <- list(coefficients = fit$coefficients, naive.var = fit$naive.var, var = var,  linear.predictors = lp, loglik = fit$loglik, score = fit$score,
                     rscore = rscore, wald.test = wald.test, score.residual =  as.matrix(Stheta), iter = fit$iter,
                     weights = fit$weights, basehaz = basehaz, Ithealp = Ithealp, model.missing = model_missing, n = n, nevent = nevent,
                     nc = ndata, ncevent = nnevent, subtype = list(n_subtype = n_subtype, marker_name = marker_name, total_subtype = ototal_subtype, marker_rr = marker_rr), formula = formula, call = Call, terms = fit$terms, assign = fit$assign, method = "IPW")

        if(rmodel){afit$model <- mf}
        if (rx)  {
          afit$x <- fit$x
        }
        if (length(strata)) {
          afit$strata <- strata
        }
        if (ry)  afit$y <- fit$y
        afit$xlevels <- fit$xlevels
        afit$offset <- fit$offset

        class(afit) <- "IPWcprisk"
    }
    return(afit)

}



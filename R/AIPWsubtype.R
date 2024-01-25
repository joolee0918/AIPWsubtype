## Modification of Therneau T (2015). _A Package for Survival Analysis in S_. version 2.38, <URL: https://CRAN.R-project.org/package=survival>.*/

#' Augmented inverse probability weighted Cox proportional hazard model for competing risks data
#'
#' Fitting an augmented inverse probability weighted Cox proportional hazard model for competing risks data with partially missing markers, in which marker variables define the subyptes of outcome.
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
#' @param control an object of class \code{\link[survival]{coxph.control}} in \code{survival} packages.The default value of \code{iter.max} is 2000 and that of \code{eps} is 1e-12. See \code{\link[survival]{coxph.control}} for other values.

#' @details The Cox proportional hazard model is used to model cause-specific hazard functions. To examine the association between exposures and the specific subtypes of disease, the log-linear is used for reparameterization. Logistic regression models \code{\link[stats]{glm}} are used for the missiness models and a conditional logistic regression model \code{\link[survival]{clogit}} is used for the marker model.
#' The data duplication method is used so that the returned value \code{x} and \code{y} are duplicated the number of subtypes of outcome. Special terms including \code{+strata()} and \code{+offset()} can be used. \code{+cluster()} should be avoided since we automatically include it in the formula. Breslow method is used for handling tied event times. The first order contrasts are included as default in modeling cause-specific baseline functions based on log-linear representation.
#'
#' For marker variables, 0 indicates censored events or missing values.
#'formula = formula, call = Call, terms = Terms, assign = assign, method = "AIPW"

#' @return An object of class \code{AIPWsubtype} representing the fit.
#' \item{coefficients}{the estimated regressoin coefficients.}
#' \item{naive.var}{the inverse of estimated Information matrix.}
#' \item{var}{the robust sandwich variance estimate.}
#' \item{linear.predictors}{a vector of linear predictors. This is not centered.}
#' \item{loglik}{a vector of length 2 containing the log-likelihood with the initial values and with the estimated coefficients.}
#' \item{score}{a value of the score test at the initial value of the coefficients.}
#' \item{rscore}{a value of the robust log-rank statistic.}
#' \item{wald.test}{a value of the Wald test statistics for whether the estimated coefficients are different from the initial value of the coefficients.}
#' \item{score.residual}{the score residuals for each subject.}
#' \item{iter}{a number of iterations used.}
#' \item{conv}{an integer code for the convergence. 0 indicates successful convergence, 1 indicates a failure to convergence, and 2 indicates it reaches the maximum number of iterations.}
#' \item{basehaz}{estimated baseline cause-specific hazard functions the reference disease subtype corresponding to marker variables equal to 1.}
#' \item{Ithealp}{a matrix of the partial derivative of the score functions with respect to the parameters from the missingness models.}
#' \item{Ithegam}{a matrix of the partial derivative of the score functions with respect to the parameters from the marker model.}
#' \item{model.missing}{a list of an object of class \code{glm} fitting the missingness models.}
#' \item{model.subtype}{an object of class \code{clogit} fitting the marker model.}
#' \item{n}{the number of observations.}
#' \item{nevent}{the number of events.}
#' \item{subtype}{a list of values related to subtypes including the number of subtypes, character strings of marker names, etc.}
#'
#'
#' The object will also contain the following: strata, formula, call, terms, y, xlevels, offset, optionally x, and model.

#' @importFrom stats model.frame model.extract model.matrix model.offset aggregate update.formula printCoefmat rexp qlnorm rnorm pchisq pnorm predict glm
#' @import survival
#' @importFrom Matrix bdiag
#' @importFrom data.table rbindlist
#' @importFrom sandwich estfun
#' @importFrom mlogit mlogit mlogit.data mFormula
#'
#' @examples
#' m1 <- AIPWsubtype(Surv(start, time, status)~ X + W,  data = data, id = "id", missing_model = "multinom",  missing_formula = list(~time + X, ~time + X),
#'  two_stage = FALSE, marker_name = c("y1", "y2"), second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = "W")
#'
#' # Two-stage missingness model
#' m2 <- AIPWsubtype(Surv(start, time, status)~ X + W,  data = data, id = "id", missing_model = "multinom", missing_formula = list(~time + X, ~time + X, ~time + X + W),
#'  two_stage = TRUE, tstage_name = c("R0"), marker_name = c("y1", "y2"), second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = "W")
#'
#'
#' @export
AIPWsubtype <- function(formula, data, id, missing_model = c("condi", "multinom"), missing_formula, marker_formula, missing_indep = FALSE, two_stage = FALSE, tstage_name = NULL,
     marker_name, marker_rr = NULL, first_cont_rr = TRUE, second_cont_bl = FALSE,  second_cont_rr = FALSE, constvar = NULL, init, control, x = FALSE, y = TRUE, model = FALSE) {


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


    if (missing(init))
      init <- NULL

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

    cause <- rep(NA, n)
    cause <- findcause(R, cause, data[, event], as.matrix(marker), on_subtype, as.matrix(ototal_subtype))

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
    eventrid <- edata[, "rowid"]
    nevent <- nrow(edata)
    obseid <- data[R == 1 & data[, event] == 1, id]

    ## Drop id for incomplete data
    dropid <- data[R != 1 & data[, event] == 1, id]

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


    ## pi

    pR <- matrix(rep(c(1, rep(0, onR - 1)), nevent), byrow = T, nrow = nevent, ncol = onR)

    if(missing_model == "multinom"){
        tmpf <- mFormula(mult_formula)

        if(two_stage == FALSE){

          tmpdata <- mlogit.data(edata, shape="wide", choice="R")
          pR <- predict(model_missing[[1]], type="probs", newdata=tmpdata)

        }else {
          edata$RR <- ifelse(edata$R == max(oR), 1, edata$R)
          tmpdata <- mlogit.data(edata, shape="wide", choice="RR")
          pR <- predict(model_missing[[1]], type="probs", newdata=tmpdata)

          pR <-  cbind(pR, as.vector(1 - predict(model_missing[[2]], newdata = edata, type = "response")))
          pR[, 1:(onR-1)] <- pR[, 1:(onR-1)]*(1-pR[, onR])
        }

        nalp <- length(coef(model_missing[[1]])) + length(model_missing[[2]]$coefficients)

        nalp0 = 0
        if (two_stage == T)
          nalp0 = length(model_missing[[2]]$coefficients)

        newX <- model.matrix(tmpf, tmpdata)
        tmpX <-  model.matrix(missing_formula[[1]], edata)
        scoreX <- matrix(exp(as.vector(newX%*%coef(model_missing[[1]]))), byrow=T, ncol=fonR)[,-1]

        dpR <- lapply(1:nevent, function(i) dpR_multinom(tmpX[i,], scoreX[i,], onR, two_stage, ncol(scoreX), nalp, nalp0))

        if (two_stage == T) {
           tmp <- -model.matrix(model_missing[[2]]$formula, edata) *(1-pR[,onR])/(1 + exp(predict(model_missing[[2]], newdata = edata)))

           for(i in 1:nevent) {
            dpR[[i]][onR, (nalp - nalp0 + 1):(nalp)] <- tmp[i,]
            dpR[[i]][1:(onR-1),] <- dpR[[i]][1:(onR-1),]*(1-pR[i, onR])
            dpR[[i]][1:(onR-1), (nalp - nalp0 + 1):(nalp)] <- t(-tmp[i,] %o% pR[i, -onR]/(1-pR[i, onR]))
           }
        }

    } else{

    est_pi = list()

    for (r in 1:fonR) {
       l <- ifelse(r == 1, 1, 0)
        PI <- matrix(l, nrow = nevent, ncol = nvecR)
        est_pi[[r]] <- matrix(0, nrow = nevent, ncol = nvecR)
        for (k in 1:nvecR) {
            tmp <- ifelse(rep(ototal_R[r, k], nevent) == 1, predict(model_missing[[k]], newdata = edata,
                type = "response"), 1 - predict(model_missing[[k]], newdata = edata, type = "response"))
            est_pi[[r]][, k] <- tmp
            PI[, k] <- tmp
        }
        pR[, r] <- apply(PI, 1, prod)
    }
    if (two_stage == T) {
        est_pi[[nR]] <- 1 - predict(model_missing[[nvecR]], newdata = edata, type = "response")
        pR[, nR] <- est_pi[[nR]]
    }


    nalp = 0
    for (k in 1:nvecR) {
        nalp = nalp + length(model_missing[[k]]$coefficients)
    }
    nalp0 = 0
    if (two_stage == T)
        nalp0 = length(model_missing[[nvecR]]$coefficients)

    dpR <- replicate(nevent, matrix(0, nrow = nR, ncol = nalp), simplify = F)

    for (r in 1:fonR) {
        a = 1
        tt <- matrix(0, nrow = nevent, ncol = nalp)

           for (k in 1:nvecR) {
            b <- length(model_missing[[k]]$coefficients)
            tmp <- model.matrix(model_missing[[k]]$formula, edata) * predict(model_missing[[k]], newdata = edata,
                type = "response")/(1 + exp(predict(model_missing[[k]], newdata = edata))) * apply(as.matrix(est_pi[[r]][,
                -k]), 1, prod)
            if (ototal_R[r, k] == 0)
                tmp <- -tmp

            tt[, a:(a + b - 1)] <- tmp

            for (i in 1:nevent) dpR[[i]][r, a:(a + b - 1)] <- tt[i, a:(a + b - 1)]
            a = a + b
        }
    }
    if (two_stage == T) {
        tt <- matrix(0, nrow = nevent, ncol = length(model_missing[[nvecR]]$coefficients))
        tmp <- -model.matrix(model_missing[[nvecR]]$formula, edata) * predict(model_missing[[nvecR]], newdata = edata,
            type = "response")/(1 + exp(predict(model_missing[[nvecR]], newdata = edata)))
        tt <- tmp

        for (i in 1:nevent) dpR[[i]][nR, (nalp - nalp0 + 1):(nalp)] <- tt[i, ]
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
    newcause <- unlist(lapply(1:n, function(i) lf(cause[i])))
    newdata <- data[rep(1:n, each = on_subtype), ]
    newdata[, event] <- rep(0, n*on_subtype)
    newdata[seq(1, n*on_subtype, by=on_subtype), event] <- data[, event]
    newdata[, marker_name] <- data.frame(ototal_subtype[newcause, ])
    newid <- newdata[, id]


    ## Do not duplicate: pr, dpr,

    marker <- as.data.frame(marker[rep(seq_len(nrow(marker)), each = on_subtype), ])
    names(marker) <- marker_name
    R <- R[rep(seq_len(length(R)), each = on_subtype)]

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


   if(n_marker >1) pairm <- combn(n_marker, 2)

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



    
    if(marker_formula==FALSE) {
         Xformula <- as.formula(paste("~", order_bl))

         Xterms <- terms(Xformula)
         eid <- newdata[newdata[, event]==1, id]
         subset_data <- newdata[ newdata[, id] %in% eid & newdata[, id] %in% obseid, ]
         s_X <- model.matrix(Xformula, data = subset_data)

         drop_bl <- NULL
         if (second_cont_bl == FALSE & second_cont_rr == TRUE) {
           drop_bl <- rep(0, ncol(pairm))
           for (i in 1:ncol(pairm)) {
             drop_bl[i] <- paste(term_marker[pairm[1, i]], term_marker[pairm[2, i]], sep = ":")
           }
           order_bl <- paste("-", paste(drop_bl, collapse = "-"))
         }

         whichdrop <- which(dimnames(attr(Xterms, "factors"))[[2]] %in% c(unconstvar, drop_bl))
         xdrop <- attributes(s_X)$assign %in% c(0, whichdrop)
         s_X <- s_X[, !xdrop, drop = FALSE]
         s_y <- subset_data[, event]
         s_uid <- subset_data[, id]
         model_subtype <- clogit(s_y ~ s_X + strata(s_uid))  
         }
     else {
          Xformula <- as.formula(paste("~", order_bl, order_rr))

         Xterms <- terms(Xformula)
         eid <- newdata[newdata[, event]==1, id]
         subset_data <- newdata[ newdata[, id] %in% eid & newdata[, id] %in% obseid, ]
         s_X <- model.matrix(Xformula, data = subset_data)

         drop_bl <- NULL
         if (second_cont_bl == FALSE & second_cont_rr == TRUE) {
           drop_bl <- rep(0, ncol(pairm))
           for (i in 1:ncol(pairm)) {
             drop_bl[i] <- paste(term_marker[pairm[1, i]], term_marker[pairm[2, i]], sep = ":")
           }
           order_bl <- paste("-", paste(drop_bl, collapse = "-"))
         }

         whichdrop <- which(dimnames(attr(Xterms, "factors"))[[2]] %in% c(unconstvar, drop_bl))
         xdrop <- attributes(s_X)$assign %in% c(0, whichdrop)
         s_X <- s_X[, !xdrop, drop = FALSE]
         s_y <- subset_data[, event]
         s_uid <- subset_data[, id]

          model_subtype <- clogit(s_y ~ s_X + strata(s_uid))
     }
    subset_data$lp <- model_subtype$linear.predictors

    gamma <- model_subtype$coefficients
    ngamma <- length(gamma)

    newformula <- update.formula(formula, paste("~.+", order_bl, order_rr, "+", "cluster", "(", id, ")"))
    mf <- model.frame(newformula, data = newdata)
    special <- c("strata", "cluster")
    Terms <- terms(newformula, special)


    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    type <- attr(Y, "type")
    if (type != "right" && type != "counting")
        stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))

    data.n <- nrow(Y)  #remember this before any time transforms

    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        stemp <- survival:::untangle.specials(Terms, "strata", 1)
        if (length(stemp$vars) == 1)
            strata.keep <- mf[[stemp$vars]] else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }

    cluster <- attr(Terms, "specials")$cluster

    if (length(cluster)) {
        robust <- TRUE  #flag to later compute a robust variance
        tempc <- survival:::untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1))
            stop("Cluster can not be used in an interaction")
        cluster <- strata(mf[, tempc$vars], shortlabel = TRUE)  #allow multiples
        dropterms <- tempc$terms  #we won't want this in the X matrix
        # Save away xlevels after removing cluster (we don't want to save upteen levels of that variable, which we
        # will never need).
        xlevels <- .getXlevels(Terms[-tempc$terms], mf)
    } else {
        dropterms <- NULL
        if (missing(robust))
            robust <- FALSE
        xlevels <- .getXlevels(Terms, mf)
    }

    contrast.arg <- NULL  #due to shared code with model.matrix.coxph
    attr(Terms, "intercept") <- 1
    stemp <- untangle.specials(Terms, "strata", 1)
    hasinteractions <- FALSE
    if (length(stemp$vars) > 0) {
        # if there is a strata statement multiple strata terms are allowed The factors attr has one row for each
        # variable in the frame, one col for each term in the model.  Pick rows for each strata var, and find if it
        # participates in any interactions.
        for (i in stemp$vars) {
            if (any(attr(Terms, "order")[attr(Terms, "factors")[i, ] > 0] > 1))
                hasinteractions <- TRUE
        }
        if (!hasinteractions)
            dropterms <- c(dropterms, stemp$terms)
    }

    if (length(dropterms)) {
        Terms2 <- Terms[-dropterms]
        X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
        # we want to number the terms wrt the original model matrix
        temp <- attr(X, "assign")
        shift <- sort(dropterms)
        temp <- temp + 1 * (shift[1] <= temp)
        if (length(shift) == 2)
            temp + 1 * (shift[2] <= temp)
        attr(X, "assign") <- temp
    } else {
        Terms2 <- Terms
        X <- model.matrix(Terms2, mf, contrasts = contrast.arg)
    }


    # drop intercept and drop baseline part if necessary
    Xattr <- attributes(X)
    whichdrop <- which(dimnames(attr(Terms2, "factors"))[[2]] %in% c(drop_bl))
    if (hasinteractions) {
        xdrop <- Xattr$assign %in% c(0, whichdrop, untangle.specials(Terms2, "strata")$terms)
    } else {
        xdrop <- Xattr$assign %in% c(0, whichdrop)

    }
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xattr$assign[!xdrop]


    # where is the main effect for unconstrained variables in X matrix
    whichX <- which(dimnames(attr(Terms2, "factors"))[[2]] %in% c(unconstvar))
    whereX <- which( attr(X, "assign")  %in% whichX)
    # where is the main effect for uonstrained variables in X matrix
    whichW <- which(dimnames(attr(Terms2, "factors"))[[2]] %in% c(constvar))
    whereW <- which( attr(X, "assign")  %in% whichW)

    offset <- model.offset(mf)
    if (is.null(offset) | all(offset == 0)) {
        offset <- rep(0, nrow(mf))
    } else {
        if (any(!is.finite(offset)))
            stop("offsets must be finite")
    }


    # new marker|R by model matrix

    newmarker <- as.matrix(model.matrix.lm(as.formula(paste("~", paste(term_marker, collapse = "+"))), data = marker,
        na.action = na.pass)[, -1])
    nc_marker <- ncol(newmarker)

    ntotal_subtype <- as.matrix(model.matrix.lm(as.formula(paste("~", paste(term_marker, collapse = "+"))), data = ototal_subtype,
                                 na.action = na.pass)[, -1])


    tmp <- list()
    for (i in 1:n_marker) {
        tmp[[i]] <- ototal_R[, rep(i, length(level_y[[i]]) - 1)]
    }

    ntotal_R <- as.matrix(do.call(cbind, tmp))


    if (two_stage == T) {
        nmarker_r <- replicate(onR - 1, ntotal_subtype, simplify = F)
        for (i in 2:onR) {
            # exclude R=1; complete case
            for (k in 1:on_subtype) {
                nmarker_r[[i - 1]][k, ] <- as.vector(ntotal_subtype[k, ]) * as.vector(1 - ntotal_R[i, ])
            }
        }
        for (i in 2:onR) {
            # exclude R=1; complete case
            nmarker_r[[i - 1]] <- as.matrix(unique(nmarker_r[[i - 1]]))
        }
    } else {
        nmarker_r <- replicate(onR - 1, ntotal_subtype, simplify = F)
        for (i in 2:onR) {
            # exclude R=1; complete case
            for (k in 1:on_subtype) {
                nmarker_r[[i - 1]][k, ] <- as.vector(ntotal_subtype[k, ]) * as.vector(1 - ntotal_R[i, ])
            }
        }
        for (i in 2:onR) {
            # exclude R=1; complete case
            nmarker_r[[i - 1]] <- as.matrix(unique(nmarker_r[[i - 1]]))
        }

    }

    if (sum(Y[, ncol(Y)]) == 0) {
        stop("No events in the data!")
    }

    if (type == "right") {
        fitter <- get("coxph.fit")
    } else {
        fitter <- get("agreg.fit")
    }

    if (missing(init)) init <- NULL

    ## Fitting model

    if (type == "right") {
        fitter <- get("AIPW_coxph.fit")
    } else {
        fitter <- get("AIPW_agreg.fit")
    }


    fit <- fitter(x = X, y = Y, eventid = eventid, id = newid, strata = strats, offset = offset, whereX = whereX,
        whereW = whereW, init = init, control = control, marker = newmarker, gamma = gamma, pR = pR, R = R,
        dpR = dpR, nR = onR, total_R = ntotal_R, marker_r = nmarker_r, two_stage = two_stage, n_marker = nc_marker, first_cont_rr = first_cont_rr, second_cont_rr = second_cont_rr,
        second_cont_bl = second_cont_bl, rownames = row.names(mf), collapse = cluster)


    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "AIPWcprisk"
    } else {
        if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
            vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
            msg <- paste("X matrix deemed to be singular; variable", paste(vars, collapse = " "))
            if (singular.ok)
                warning(msg) else stop(msg)
        }



        Igam <- vcov(model_subtype)


        Salp = as.data.frame(uniqid)
        names(Salp) <- id

        tmp_id = edata[, id]

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
          Ualp <- cbind(tmp_id, Ualp)
          colnames(Ualp)[1] <- id
        } else {
          Ualp <- cbind(edata[edata[, tstage_name] == 1, id], Ualp)
          colnames(Ualp)[1] <- id
          Ualp_ts <- cbind(tmp_id, Ualp_ts)
          colnames(Ualp_ts)[1] <- id
          Ualp <- suppressWarnings(merge(Ualp, Ualp_ts, by = id, all = T))
          Ualp[is.na(Ualp)] <- 0
        }

        Salp <- suppressWarnings(merge(Salp, Ualp, by = id, all = T))
        Salp[is.na(Salp)] <- 0
        Salp <- as.matrix(Salp[, -1])

        Sgam <- as.data.frame(uniqid)
        names(Sgam) <- id

        subset_data <- cbind(subset_data, s_X)
        s0 <- by(subset_data, subset_data[, id], function(x) sum(exp(x$lp)))
        s1 <- by(subset_data, subset_data[, id], function(x) t(x[, colnames(s_X)]) %*% exp(x$lp))
        term1 <- by(subset_data, subset_data[, id], function(x) x[, colnames(s_X)][x[, event] == 1, ])
        Ugam <- lapply(1:length(term1), function(i) term1[[i]] - s1[[i]]/s0[[i]])
        Ugam <- as.matrix(rbindlist(Ugam))

        Ugam <- cbind(as.numeric(names(term1)), Ugam)
        colnames(Ugam)[1] <- id

        Sgam <- suppressWarnings(merge(Sgam, Ugam, by = id, all = T))
        Sgam[is.na(Sgam)] <- 0
        Sgam <- as.matrix(Sgam[, -1])

        resid <- fit$resid - as.matrix(Sgam) %*% Igam %*% t(fit$Ithegam) - as.matrix(Salp) %*% Ialp %*% t(fit$Ithealp)
        var <- fit$naive.var %*% t(resid) %*% resid %*% fit$naive.var

        ## robust log-rank statistic
        temp0 <- fit$score0 - as.matrix(Sgam) %*% Igam %*% t(fit$Ithegam) - as.matrix(Salp) %*% Ialp %*% t(fit$Ithealp)
        u <- apply(as.matrix(temp0), 2, sum)
        rscore <- survival::coxph.wtest(t(temp0)%*%temp0, u, control$toler.chol)$test


        #Wald test
        if (length(fit$coefficients)) {
          #not for intercept only models, or if test is already done
          nabeta <- !is.na(fit$coefficients)
          # The init vector might be longer than the betas, for a sparse term
          if (is.null(init)) temp <- fit$coefficients[nabeta]
          else temp <- (fit$coefficients -
                          init[1:length(fit$coefficients)])[nabeta]
          wald.test <-  survival::coxph.wtest(var[nabeta,nabeta], temp,
                                                   control$toler.chol)$test
        }

        if(is.null(strats)) {
          stratum <- rep(1, nrow(Y))
        }else {
          stratum <- strats
        }
        s0 <- exp(fit$linear.predictors)
        basehaz <- baseHaz(Y, stratum, s0)

        afit <- list(coefficients = fit$coef, naive.var = fit$naive.var, var = var, linear.predictors = fit$linear.predictors,
                     score = fit$sctest, loglik = fit$loglik, rscore = rscore, wald.test = wald.test, score.residual = fit$resid, iter = fit$iter, conv = fit$conv, basehaz = basehaz,
                     Ithealp = fit$Ithealp, Ithegam = fit$Ithegam, model.missing = model_missing, model.subtype = model_subtype,
                     n = n, nevent = nevent, subtype = list(n_subtype = on_subtype, marker_name = marker_name, total_subtype = ototal_subtype, marker_rr = marker_rr),
                     formula = formula, call = Call, terms = Terms, assign = assign, method = "AIPW")


        if(rmodel){afit$model <- mf}
        if (rx)  {
          afit$x <- X
        }
        if (length(strats)) {
          afit$strata <- strata.keep
        }
        if (ry)  afit$y <- Y

        if (length(xlevels) >0) afit$xlevels <- xlevels
        if (any(offset !=0)) afit$offset <- offset
        class(afit) <- "AIPWcprisk"
    }

    return(afit)

}


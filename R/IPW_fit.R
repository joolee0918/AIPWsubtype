
#' @export
IPW_fit <- function(formula, data, id, missing_model, missing_indep = FALSE, two_stage = FALSE, tstage_name = NULL, n_marker, marker_name,
    second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = NULL, init, control, ...) {


    Call <- match.call()
    data <- data[order(data[, id]), ]
    n <- nrow(data)
    marker <- data[, marker_name]
    marker[marker == 0] <- NA


    for (i in 1:n_marker) {
        marker[, i] <- factor(marker[, i])
    }
    data[, marker_name] <- marker

    n_subtype = 1
    for (i in 1:n_marker) {
        n_subtype = n_subtype * (nlevels(factor(data[, marker_name[i]])))
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
        tmpr <- rep(list(c(1, 0)), n_marker)
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

    # possible marker combination

    level_y = list()
    for (k in 1:n_marker) {
        level_y[[k]] <- seq(nlevels(factor(data[, marker_name[k]])), 1)
    }
    tmpy = list()
    for (k in 1:n_marker) {
        tmpy[[k]] <- as.vector(level_y[[k]])
    }
    total_subtype <- as.data.frame(expand.grid(tmpy))


    # possible marker combination

    event <- tail(survival:::terms.inner(formula[1:2]), 1)
    Rmat <- matrix(1, nrow = n, ncol = n_marker)
    for (i in 1:n_marker) {
        Rmat[, i] <- ifelse(data[, event] == 1, 1 - as.integer(is.na(marker[, i])), 1)
    }

    R = rep(1, n)
    if (two_stage == T)
        R[data[, tstage_name] == 0] <- nR
    for (i in 1:n) {
        if (R[i] == 1) {
            for (j in 1:2^n_marker) {
                if (all(Rmat[i, ] == total_R[j, 1:n_marker]))
                  R[i] <- j
            }
        }
    }


    ## Data frame for missing indicator

    # model missing
    missing_data <- as.data.frame(Rmat)
    names(missing_data) <- c(paste("R", c(1:n_marker), sep = ""))
    if (two_stage == TRUE) {
        missing_data$R0 <- data[, tstage_name]
    }

    missing_name <- names(missing_data)
    data <- cbind(data, missing_data)
    edata <- data[data[, event] == 1, ]
    eventid <- edata[, id]
    nevent <- nrow(edata)
    uniqid <- unique(data[, id])


    model_missing <- list()
    tmp_model <- as.formula(paste(missing_name[1], paste(missing_model[[1]], collapse = "")))

    if (two_stage == FALSE) {

        for (i in 1:n_marker) {
            model_missing[[i]] <- glm(tmp_model, data = edata, family = binomial)
            if (i == n_marker)
                break
            tmp_mname <- missing_name[1:i]
            if(missing_indep == FALSE) missing_model[[i + 1]] <- paste(paste(missing_model[[i + 1]], collapse = ""), paste(tmp_mname, collapse = "+"),
                sep = "+")
            tmp_model <- as.formula(paste(missing_name[i + 1], paste(missing_model[[i + 1]], collapse = "")))
        }
    } else {
        for (i in 1:n_marker) {
            model_missing[[i]] <- glm(tmp_model, data = edata[edata[, tstage_name] == 1, ], family = binomial)
            if (i == n_marker)
                break
            tmp_mname <- c(missing_name[1:i])
            if(missing_indep == FALSE) missing_model[[i+1]] <- paste(paste( missing_model[[i+1]], collapse=''), paste(tmp_mname, collapse ='+'), sep='+')
            tmp_model <- as.formula(paste(missing_name[i + 1], paste(missing_model[[i + 1]], collapse = "")))
        }

        tmp_model <- as.formula(paste(tstage_name, paste(missing_model[[nvecR]], collapse = "")))

        model_missing[[nvecR]] <- glm(tmp_model, data = edata, family = binomial)

    }

    newdata <- data[R == 1, ]
    ndata <- nrow(newdata)
    missing_prob <- rep(1, ndata)
    newedata <- newdata[newdata[, event] == 1, ]
    nnevent <- nrow(newedata)
    nuniqid <- unique(newdata[, id])
    newmarker <- marker[R == 1, ]


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
        PI[newdata[, event] == 1, k] <- tmp
    }
    missing_prob <- apply(PI, 1, prod)

    dpR <- as.data.frame(matrix(0, ncol = nalp, nrow = ndata))

    a = 1
    for (k in 1:nvecR) {
        b <- length(model_missing[[k]]$coefficients)
        tmp <- model.matrix(model_missing[[k]]$formula, newedata) * predict(model_missing[[k]], newdata = newedata,
            type = "response")/(1 + exp(predict(model_missing[[k]], newdata = newedata))) * apply(as.matrix(est_pi[,
            -k]), 1, prod)
        dpR[newdata[, event] == 1, a:(a + b - 1)] <- tmp
        a = a + b
    }

    dpR <- dpR[rep(seq_len(nrow(dpR)), each = n_subtype), ]

    ## Data frame for cause

    cause <- rep(0, ndata)
    for (i in 1:ndata) {
        if (anyNA(newmarker[i, ])) {
            cause[i] <- NA
        } else {
            for (j in 1:n_subtype) {
                if (all(newmarker[i, ] == total_subtype[j, 1:n_marker]))
                  cause[i] <- j
            }
        }
    }

    newdata$cause <- cause
    newdata$pi1 <- missing_prob
    newdata <- lapply(1:ndata, function(i) newdata[rep(i, each = n_subtype), ])
    lf <- function(x) {
        if (!is.na(x$cause[1])) {
            x$cause <- c(x$cause[1], seq(1, n_subtype)[!(seq(1, n_subtype) %in% x$cause[1])])
        } else {
            x$cause <- seq(1, n_subtype)
        }
        x[, event] <- c(x[, event][1], rep(0, (n_subtype - 1)))
        x
    }

    newdata <- lapply(newdata, lf)
    newdata <- as.data.frame(rbindlist(newdata))

    newdata[, marker_name] <- data.frame(total_subtype[newdata$cause, ])
    term_marker <- rep(0, n_marker)

    for (i in 1:n_marker) term_marker[i] <- paste("factor(", marker_name[i], ")", sep = "")

    special <- c("strata", "cluster")
    Tf <- terms(formula, special)
    strats <- attr(Tf, "specials")$strata
    if (length(strats)) {
        strterm <- survival:::untangle.specials(Tf, "strata", 1)$terms
        Xattr <- attr(Tf, "term.labels")[-strterm]
    } else {
        Xattr <- attr(Tf, "term.labels")
    }

    unconstvar <- Xattr[!(Xattr %in% constvar)]

    nuvar <- length(unconstvar)
    order_rr <- NULL

    if (nuvar == 0) {
        order_rr <- paste(term_marker, collapse = "+")
    } else {
        for (i in 1:nuvar) {
            tmp_rr <- paste(term_marker, "*", unconstvar[i], collapse = "+")
            order_rr <- paste(order_rr, tmp_rr, sep = "+")
        }
    }


    pairm <- combn(n_marker, 2)

    if (second_cont_rr == TRUE) {
        order_rr <- NULL
        for (i in 1:nuvar) {
            tmp_rr <- paste(attr(terms(as.formula(paste("~", unconstvar[i], "*", "(", paste(term_marker, collapse = "+"),
                ")", "^", 2))), "term.labels"), collapse = "+")
            order_rr <- paste(order_rr, tmp_rr, sep = "+")
        }
    }

    order_bl <- NULL
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
    m_c <- coxph(formula = newformula, data = newdata, robust = T, method = "breslow")


    fit <- coxph(formula = newformula, data = newdata, weights = 1/pi1, robust = T, model = TRUE, x = TRUE,
        method = "breslow")
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
        vv <- fit$naive.var
        weights <- fit$weights
        strata <- fit$strata
        lp <- fit$linear.predictors + sum(fit$coefficients * fit$means)

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

        # P I_thealp = t(tmp_score)%*%as.matrix(dpR)

        ## Only fully observed data

        Itheta <- fit$naive.var
        RVar_theta <- fit$var

        # For all data

        Ialp <- vcov(model_missing[[1]])
        if (nvecR > 1) {
            for (k in 2:nvecR) {
                Ialp <- bdiag(Ialp, vcov(model_missing[[k]]))
            }
        }
        Ialp <- as.matrix(Ialp)

        Salp <- as.data.frame(uniqid)
        names(Salp) <- id

        Ualp <- list()
        tmp_id <- edata[, id]
        for (k in 1:n_marker) {
            Ualp[[k]] = estfun(model_missing[[k]])
        }
        if (two_stage == TRUE) {
            Ualp_ts = estfun(model_missing[[nvecR]])
        }
        Ualp <- do.call(cbind, Ualp)
        if (two_stage == FALSE) {
            Ualp <- cbind(tmp_id, Ualp)
            colnames(Ualp)[1] <- id
        } else {
            Ualp <- cbind(edata[edata[, tstage_name] == 1, id], Ualp)
            colnames(Ualp)[1] <- id
            Ualp_ts <- cbind(tmp_id, Ualp_ts)
            colnames(Ualp_ts)[1] <- id
            Ualp <- merge(Ualp, Ualp_ts, by = id, all = T)
            Ualp[is.na(Ualp)] <- 0
        }

        Salp <- merge(Salp, Ualp, by = id, all = T)
        Salp[is.na(Salp)] <- 0
        Salp <- as.matrix(Salp[, -1])


        resid <- as.matrix(Stheta)
        resid <- as.matrix(Stheta) - as.matrix(Salp) %*% Ialp %*% t(Ithealp)
        var_IPW <- fit$naive.var %*% t(resid) %*% resid %*% fit$naive.var

        afit <- list(coefficients = fit$coefficients, var = var_IPW, loglik = fit$loglik, iter = fit$iter,
            weights = fit$weights, Ithealp = Ithealp, model_missing = model_missing, res = resid, n = n, nevent = nevent,
            ndata = ndata, nnevent = nnevent, call = Call, terms = fit$terms, assign = fit$assign, method = "IPW")

        class(afit) <- "IPWcprisk"
    }
    return(afit)

}



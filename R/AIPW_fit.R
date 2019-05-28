
#' @importFrom stats model.frame model.extract model.matrix model.offset aggregate update.formula printCoefmat rexp qlnorm rnorm pchisq pnorm predict glm
#' @import survival
#' @importFrom Matrix bdiag
#' @importFrom data.table rbindlist
#' @importFrom sandwich estfun
#'
#' @export
AIPWsubtype <- function(formula, data, id, missing_model, missing_indep = FALSE, two_stage, tstage_name = NULL,
    n_marker, marker_name, second_cont_bl = FALSE, second_cont_rr = FALSE, constvar = NULL, init, control, ...) {


    Call <- match.call()
    data <- data[order(data[, id]), ]

    n <- nrow(data)

    uniqid <- unique(data[, id])
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

    # possible marker combination

    level_y = list()
    for (k in 1:n_marker) {
        level_y[[k]] <- seq(nlevels(factor(data[, marker_name[k]])), 1)
    }
    tmpy <- list()
    for (k in 1:n_marker) {
        tmpy[[k]] <- as.vector(level_y[[k]])
    }
    total_subtype <- as.data.frame(expand.grid(tmpy))

    ## marker|R
    if (two_stage == T) {
        marker_r <- replicate(nR - 1, total_subtype, simplify = F)
        for (i in 2:nR) {
            # exclude R=1; complete case
            for (k in 1:n_subtype) {
                marker_r[[i - 1]][k, ] <- as.vector(total_subtype[k, ]) * as.vector(1 - total_R[i, 1:n_marker])
            }
        }
        for (i in 2:nR) {
            # exclude R=1; complete case
            marker_r[[i - 1]] <- as.matrix(unique(marker_r[[i - 1]]))
        }
    } else {
        marker_r <- replicate(nR - 1, total_subtype, simplify = F)
        for (i in 2:nR) {
            # exclude R=1; complete case
            for (k in 1:n_subtype) {
                marker_r[[i - 1]][k, ] <- as.vector(total_subtype[k, ]) * as.vector(1 - total_R[i, 1:n_marker])
            }
        }
        for (i in 2:nR) {
            # exclude R=1; complete case
            marker_r[[i - 1]] <- as.matrix(unique(marker_r[[i - 1]]))
        }

    }

    ### Marker

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
    obseid <- data[R == 1 & data[, event] == 1, id]

    ## Drop id for incomplete data
    dropid <- data[R != 1 & data[, event] == 1, id]


    model_missing <- list()
    tmp_model <- as.formula(paste(missing_name[1], paste(missing_model[[1]], collapse = "")))

    if (two_stage == FALSE) {

        for (i in 1:n_marker) {
            model_missing[[i]] <- glm(tmp_model, data = edata, family = binomial)
            if (i == n_marker)
                break
            tmp_mname <- missing_name[1:i]
            if (missing_indep == FALSE)
                missing_model[[i + 1]] <- paste(paste(missing_model[[i + 1]], collapse = ""), paste(tmp_mname,
                  collapse = "+"), sep = "+")
            tmp_model <- as.formula(paste(missing_name[i + 1], paste(missing_model[[i + 1]], collapse = "")))
        }
    } else {
        for (i in 1:n_marker) {
            model_missing[[i]] <- glm(tmp_model, data = edata[edata[, tstage_name] == 1, ], family = binomial)
            if (i == n_marker)
                break
            tmp_mname <- c(missing_name[1:i])
            if (missing_indep == FALSE)
                missing_model[[i + 1]] <- paste(paste(missing_model[[i + 1]], collapse = ""), paste(tmp_mname,
                  collapse = "+"), sep = "+")
            tmp_model <- as.formula(paste(missing_name[i + 1], paste(missing_model[[i + 1]], collapse = "")))
        }

        tmp_model <- as.formula(paste(tstage_name, paste(missing_model[[nvecR]], collapse = "")))

        model_missing[[nvecR]] <- glm(tmp_model, data = edata, family = binomial)

    }



    ## pi

    NR_data <- data[data[, event] == 1, !(names(data) %in% missing_name)]
    pR <- matrix(rep(c(1, rep(0, nR - 1)), nevent), byrow = T, nrow = nevent, ncol = nR)


    est_pi = list()
    for (r in 1:2^(n_marker)) {
        RR <- matrix(total_R[r, ], byrow = T, nrow = 1)
        colnames(RR) <- colnames(total_R)
        NR_data_tmp <- cbind(NR_data, RR)
        l <- ifelse(r == 1, 1, 0)
        PI <- matrix(l, nrow = nevent, ncol = nvecR)
        est_pi[[r]] <- matrix(0, nrow = nevent, ncol = nvecR)
        for (k in 1:nvecR) {
            tmp <- ifelse(rep(total_R[r, k], nevent) == 1, predict(model_missing[[k]], newdata = NR_data_tmp,
                type = "response"), 1 - predict(model_missing[[k]], newdata = NR_data_tmp, type = "response"))
            est_pi[[r]][, k] <- tmp
            PI[, k] <- tmp
        }
        pR[, r] <- apply(PI, 1, prod)
    }
    if (two_stage == T) {
        RR <- matrix(total_R[nR, nvecR], byrow = T, nrow = 1)
        colnames(RR) <- colnames(total_R)[nvecR]
        NR_data_tmp <- cbind(NR_data, RR)
        est_pi[[nR]] <- matrix(0, nrow = nevent, ncol = 1)
        est_pi[[nR]] <- ifelse(rep(total_R[nR, nvecR], nevent) == 1, predict(model_missing[[nvecR]], newdata = NR_data_tmp,
            type = "response"), 1 - predict(model_missing[[nvecR]], newdata = NR_data_tmp, type = "response"))
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

    for (r in 1:2^(n_marker)) {
        a = 1
        tt <- matrix(0, nrow = nevent, ncol = nalp)
        RR <- matrix(total_R[r, ], byrow = T, nrow = 1)
        colnames(RR) <- colnames(total_R)
        NR_data_tmp <- cbind(NR_data, RR)

        for (k in 1:nvecR) {
            b <- length(model_missing[[k]]$coefficients)
            tmp <- model.matrix(model_missing[[k]]$formula, NR_data_tmp) * predict(model_missing[[k]], newdata = NR_data_tmp,
                type = "response")/(1 + exp(predict(model_missing[[k]], newdata = NR_data_tmp))) * apply(as.matrix(est_pi[[r]][,
                -k]), 1, prod)
            if (total_R[r, k] == 0)
                tmp <- -tmp

            tt[, a:(a + b - 1)] <- tmp

            for (i in 1:nevent) dpR[[i]][r, a:(a + b - 1)] <- tt[i, a:(a + b - 1)]
            a = a + b
        }
    }
    if (two_stage == T) {
        tt <- matrix(0, nrow = nevent, ncol = length(model_missing[[nvecR]]$coefficients))
        RR <- matrix(total_R[nR, nvecR], byrow = T, nrow = 1)
        colnames(RR) <- colnames(total_R)[nvecR]
        NR_data_tmp <- cbind(NR_data, RR)

        tmp <- -model.matrix(model_missing[[nvecR]]$formula, NR_data_tmp) * predict(model_missing[[nvecR]], newdata = NR_data_tmp,
            type = "response")/(1 + exp(predict(model_missing[[nvecR]], newdata = NR_data_tmp)))
        tt <- tmp

        for (i in 1:nevent) dpR[[i]][nR, (nalp - nalp0 + 1):(nalp)] <- tt[i, ]
    }




    ## Data frame for cause

    cause <- rep(0, n)
    for (i in 1:n) {
        if (anyNA(marker[i, ])) {
            cause[i] <- NA
        } else {
            for (j in 1:n_subtype) {
                if (all(marker[i, ] == total_subtype[j, 1:n_marker]))
                  cause[i] <- j
            }
        }
    }


    data <- cbind(data, cause = cause)
    ## Data duplication
    lf <- function(x) {
        if (!is.na(x$cause[1])) {
            x$cause <- c(x$cause[1], seq(1, n_subtype)[!(seq(1, n_subtype) %in% x$cause[1])])
        } else {
            x$cause <- seq(1, n_subtype)
        }
        x[, event] <- c(x[, event][1], rep(0, (n_subtype - 1)))
        x
    }


    newdata <- lapply(1:n, function(i) data[rep(i, each = n_subtype), ])
    newdata <- lapply(newdata, lf)
    newdata <- as.data.frame(rbindlist(newdata))
    newdata[, marker_name] <- data.frame(total_subtype[newdata$cause, ])
    newid <- newdata[, id]


    ## Do not duplicate: pr, dpr,

    marker <- marker[rep(seq_len(nrow(marker)), each = n_subtype), ]
    R <- R[rep(seq_len(length(R)), each = n_subtype)]

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

    Xformula <- as.formula(paste("~", order_bl, order_rr))

    Xterms <- terms(Xformula)
    subset_data <- newdata[newdata[, id] %in% obseid, ]
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
    subset_data$lp <- predict(model_subtype, type = "lp")

    gamma <- model_subtype$coefficients
    ngamma <- length(gamma)


    Igam <- vcov(model_subtype)
    Ialp <- vcov(model_missing[[1]])
    if (nvecR > 1) {
        for (k in 2:nvecR) {
            Ialp <- bdiag(Ialp, vcov(model_missing[[k]]))
        }
    }
    Ialp <- as.matrix(Ialp)

    Salp = as.data.frame(data[, id])
    names(Salp) <- id

    Salp = as.data.frame(uniqid)
    names(Salp) <- id
    Ualp <- list()
    tmp_id = edata[, id]
    for (k in 1:n_marker) {
        Ualp[[k]] <- estfun(model_missing[[k]])
    }
    if (two_stage == TRUE) {
        Ualp_ts <- estfun(model_missing[[nvecR]])
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

    Sgam <- merge(Sgam, Ugam, by = id, all = T)
    Sgam[is.na(Sgam)] <- 0
    Sgam <- as.matrix(Sgam[, -1])


    Salp = as.data.frame(uniqid)
    names(Salp) <- id
    Ualp <- list()
    tmp_id = edata[, id]
    for (k in 1:n_marker) {
        Ualp[[k]] <- estfun(model_missing[[k]])
    }
    if (two_stage == TRUE) {
        Ualp_ts <- estfun(model_missing[[nvecR]])
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

    Sgam <- merge(Sgam, Ugam, by = id, all = T)
    Sgam[is.na(Sgam)] <- 0
    Sgam <- as.matrix(Sgam[, -1])


    newformula <- update.formula(formula, paste("~.+", order_bl, order_rr, "+", "cluster", "(", id, ")"))
    mf <- model.frame(newformula, data = newdata)
    special <- c("strata", "cluster")
    Terms <- terms(newformula, special)

    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(coxph.control))  #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
        if (any(indx == 0L))
            stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 0L]), domain = NA)
    }

    if (missing(control)) {
        control = coxph.control()
        control$eps = 1e-12
        control$iter.max = 1000
        control$toler.inf = sqrt(control$eps)
    }

    if (missing(init))
        init <- NULL

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
        temppred <- attr(Terms, "predvars")
        Terms2 <- Terms[-dropterms]
        if (!is.null(temppred)) {
            # subscripting a Terms object currently drops predvars, in error
            attr(Terms2, "predvars") <- temppred[-(1 + dropterms)]  # 'Call' object
        }
        X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
        # we want to number the terms wrt the original model matrix Do not forget the intercept, which will be a
        # zero
        renumber <- match(colnames(attr(Terms2, "factors")), colnames(attr(Terms, "factors")))
        attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
    } else {
        X <- model.matrix(Terms, mf, contrasts = contrast.arg)
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
        X <- model.matrix(Terms, mf, contrasts = contrast.arg)
    }


    # drop intercept and drop baseline part if necessary
    Xattr <- attributes(X)
    whichdrop <- which(dimnames(attr(Terms, "factors"))[[2]] %in% c(drop_bl))
    if (hasinteractions) {
        xdrop <- Xattr$assign %in% c(0, whichdrop, untangle.specials(Terms, "strata")$terms)
    } else {
        xdrop <- Xattr$assign %in% c(0, whichdrop)

    }
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xattr$assign[!xdrop]

    # where is the main effect for unconstrained variables in X matrix
    whichX <- which(dimnames(attr(Terms, "factors"))[[2]] %in% c(unconstvar))

    # where is the main effect for uonstrained variables in X matrix
    whichW <- which(dimnames(attr(Terms, "factors"))[[2]] %in% c(constvar))

    offset <- model.offset(mf)
    if (is.null(offset) | all(offset == 0)) {
        offset <- rep(0, nrow(mf))
    } else {
        if (any(!is.finite(offset)))
            stop("offsets must be finite")
    }

    newmarker <- model.matrix.lm(as.formula(paste("~", paste(term_marker, collapse = "+"))), data = marker,
        na.action = na.pass)[, -1]


    # new marker|R by model matrix
    nc_marker <- ncol(newmarker)
    tmpy <- list()
    for (k in 1:nc_marker) {
        tmpy[[k]] <- as.vector(c(1, 0))
    }

    ntotal_subtype = as.data.frame(expand.grid(tmpy))


    tmp <- list()
    for (i in 1:n_marker) {
        tmp[[i]] <- total_R[, rep(i, length(level_y[[i]]) - 1)]
    }

    ntotal_R <- do.call(cbind, tmp)


    if (two_stage == T) {
        nmarker_r <- replicate(nR - 1, ntotal_subtype, simplify = F)
        for (i in 2:nR) {
            # exclude R=1; complete case
            for (k in 1:n_subtype) {
                nmarker_r[[i - 1]][k, ] <- as.vector(ntotal_subtype[k, ]) * as.vector(1 - ntotal_R[i, ])
            }
        }
        for (i in 2:nR) {
            # exclude R=1; complete case
            nmarker_r[[i - 1]] <- as.matrix(unique(nmarker_r[[i - 1]]))
        }
    } else {
        nmarker_r <- replicate(nR - 1, ntotal_subtype, simplify = F)
        for (i in 2:nR) {
            # exclude R=1; complete case
            for (k in 1:n_subtype) {
                nmarker_r[[i - 1]][k, ] <- as.vector(ntotal_subtype[k, ]) * as.vector(1 - ntotal_R[i, ])
            }
        }
        for (i in 2:nR) {
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


    miss <- !newdata[, id] %in% dropid
    dim(Y[miss, ])
    m_c <- fitter(X[miss, ], Y[miss, ], strats[miss], offset[miss], init = init, control, weights = NULL, method = "breslow",
        row.names(mf)[miss])
    init = m_c$coefficients

    ## Fitting model

    if (type == "right") {
        fitter <- get("AIPW_coxph.fit")
    } else {
        fitter <- get("AIPW_agreg.fit")
    }


    fit <- fitter(x = X, y = Y, eventid = eventid, id = newid, strata = strats, offset = offset, whereX = whichX,
        whereW = whichW, init = init, control = control, marker = newmarker, gamma = gamma, pR = pR, R = R,
        dpR = dpR, total_R = ntotal_R, marker_r = nmarker_r, two_stage = two_stage, n_marker = nc_marker, second_cont_rr = second_cont_rr,
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

        print(dim(Sgam))
        print(dim(Igam))
        print(dim(fit$Ithegam))
        print(dim(Salp))
        print(dim(Ialp))
        print(dim(fit$Ithealp))
        print(dim(resid))
        resid = resid - as.matrix(Sgam) %*% Igam %*% t(fit$Ithegam) - as.matrix(Salp) %*% Ialp %*% t(fit$Ithealp)

        fit$var <- fit$naive.var %*% t(resid) %*% resid %*% fit$naive.var

        afit <- list(coefficients = fit$coef, var = fit$var, loglik = fit$loglik, iter = fit$iter, conv = fit$conv,
            Ithealp = fit$Ithealp, Ithegam = fit$Ithegam, model_missing = model_missing, model_subtype = model_subtype,
            res = resid, n = n, nevent = nevent, call = Call, terms = Terms, assign = assign, method = "AIPW")

        class(afit) <- "AIPWcprisk"
    }

    return(afit)

}


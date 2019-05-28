
AIPW_coxph.fit <- function(x, y, eventid, id, strata, offset, whereX, whereW, init, control, marker, gamma,  pR, R, dpR, total_R, marker_r, two_stage, n_marker, second_cont_bl, second_cont_rr,  rownames, collapse)
{

    n <-  nrow(y)
    if (is.matrix(x)) {
        nvar <- ncol(x)
    } else {
        if (length(x)==0) {nvar <-0
        } else { nvar <-1
        }
    }
    if(two_stage==T) {
        nR = 2^n_marker + 1
    } else {
        nR = 2^n_marker
    }
    ny <- ncol(y)
    ngamma = as.integer(length(gamma))
    nalp = as.integer(dim(dpR[[1]])[2])

    time <- y[,1]
    status <- y[,2]

    pR <- as.matrix(pR)

    # Sort the data
    if(length(strata) == 0){
        sorted <- order(time)
        strata <- NULL
        newstrat <- rep(0, n)
    }else{
        sorted <- order(strata, time)
        newstrat <- c(1*(diff(as.numeric(strata[sorted]))!=0), 1)
    }

    stime <- time[sorted]
    sstat <- status[sorted]
    sx <- as.matrix(x[sorted,])
    sid <- id[sorted]
    soffset <- offset[sorted]
    sR <- R[sorted]
    smarker <- as.matrix(marker[sorted,])


    if (nvar==0) {
        # A special case: Null model.
        #  (This is why I need the rownames arg- can't use x' names)
        # Set things up for 0 iterations on a dummy variable
        x <- as.matrix(rep(1.0, n))
        nullmodel <- TRUE
        nvar <- 1
        init <- 0
        maxiter <- 0
    }else {
        nullmodel <- FALSE
        maxiter <- control$iter.max
        if (!missing(init) && length(init)>0) {
            if (length(init) != nvar) stop("Wrong length for inital values")
        }else {
            init <- rep(0,nvar)
        }
    }

    storage.mode(init) <- "double"

    comb_y = combn(n_marker,2)
    smarker[is.na(smarker)] <- -999


AIPW_fit <- AIPW_coxfit_cpp(control$iter.max,
        stime,
        sstat,
        sx,
	      eventid,
	      sid,
        soffset,
        newstrat,
        smarker,
	      sR,
	      pR,
        dpR,
	      as.matrix(total_R),
        marker_r,
        whereX,
        whereW,
        gamma,
        comb_y,
	      nvar,
	      n_marker,
        nR,
        ngamma,
        nalp,
        control$eps,
        as.double(control$toler.chol),
        second_cont_bl,
        second_cont_rr,
        init)

    var <- matrix(AIPW_fit$imat,nvar,nvar)
    coef <- AIPW_fit$coef

    infs <- abs(AIPW_fit$u %*% var)
    if (maxiter >1) {
        if (AIPW_fit$conv == 2)
        warning("Ran out of iterations and did not converge")
    }

    names(coef) <- dimnames(x)[[2]]
    lp <- c(x %*% coef)

      # Sort the data
	    if(length(strata) == 0){
	        ord <- order(y[,ny-1], -status)
	        newstrat <- as.integer(rep(0, n))
	    }else{
	        ord <- order(strata, y[,ny-1], -status)
	        newstrat <- c(diff(as.numeric(strata[ord]))!=0 ,1)
	    }
	    newstrat[n] <- 1

	    stime <-  y[ord,1]
	    sstatus <- y[ord,2]
	    sx <- x[ord,]
	    sid <- id[ord]
	    soffset <- offset[ord]
	    score <- exp( sx %*% coef + soffset)
	    sR <- R[ord]
	    smarker <- as.matrix(marker[ord,])
	    smarker[is.na(smarker)] <- -999


	resid <- AIPW_coxscore_cpp(
	        stime,
	        sstatus,
	        sx,
	        eventid,
	        sid,
	        score,
	        newstrat,
	        smarker,
		      sR,
	        pR,
	        as.matrix(total_R),
	        marker_r,
	        whereX,
	        whereW,
	        gamma,
	        comb_y,
		      nvar,
		      n_marker,
	        nR,
	        ngamma,
	        nalp,
	        second_cont_bl,
	        second_cont_rr)

	       	rr <- resid
	        rr[ord,] <- rr
	        dimnames(rr) =list(names(rownames), names(coef))

		rr <- drop(rowsum(rr, collapse))

    afit <- list(coefficients  = coef,
    var    = var,
    Ithegam = AIPW_fit$Ithegam,
    Ithealp = AIPW_fit$Ithealp,
    loglik = AIPW_fit$loglik,
    iter   = AIPW_fit$iter,
    conv   = AIPW_fit$conv,
    linear.predictors = as.vector(lp),
    resid = rr,
    method='AIPW')

    return(afit)

}


## Modification of Therneau T (2015). _A Package for Survival Analysis in S_. version 2.38, <URL: https://CRAN.R-project.org/package=survival>.*/

# =========================================================
# Fitting AIPW with (start, stop, status)
# =========================================================


AIPW_agreg.fit <- function(x, y, eventid, id, strata, offset, whereX, whereW, init, control, marker, gamma,pR, R, dpR, nR, total_R, marker_r, two_stage, n_marker, second_cont_bl, first_cont_rr, second_cont_rr,  rownames, collapse)
{
    n <-  nrow(y)
    event <- y[,3]
    if (all(event==0)) stop("Can't fit a Cox model with 0 failures")

    if (is.matrix(x)) {
        nvar <- ncol(x)
    } else {
        if (length(x)==0) {nvar <-0
        } else { nvar <-1
        }
    }

    ngamma = as.integer(length(gamma))
    nalp = as.integer(dim(dpR[[1]])[2])

    # Sort the data (or rather, get a list of sorted indices)
    #  For both stop and start times, the indices go from last to first
    if (length(strata)==0) {
        sort.end  <- order(-y[,2]) -1L #indices start at 0 for C code
        sort.start <- order(-y[,1]) -1L
        newstrat  <- n
    }
    else {
        sort.end  <- order(strata, -y[,2]) -1L
        sort.start <- order(strata, -y[,1]) -1L
        newstrat  <- cumsum(table(strata))
    }

    ny <- ncol(y)
    start <- y[,1,drop=TRUE]
    stop <- y[,2,drop=TRUE]
    status <- y[,ny,drop=TRUE]


    x <- as.matrix(x)
    pR <- as.matrix(pR)
    marker <- as.matrix(marker)
    if (nvar==0) {
        # A special case: Null model.
        #  (This is why I need the rownames arg- can't use x' names)
        # Set things up for 0 iterations on a dummy variable
        x <- matrix(as.double(1:n), ncol=1)
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
            init <- rep(0, nvar)
        }
    }


    if(n_marker >1) comb_y <- combn(n_marker,2)
    else comb_y <- as.matrix(0)
    marker[is.na(marker)] <- -999

AIPW_fit <- AIPW_agreg_cpp(maxiter,
        start,
	      stop,
	      status,
        x,
	      eventid,
	      id,
        offset,
        newstrat,
        sort.start,
        sort.end,
        marker,
        R,
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
        first_cont_rr,
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
    lp <- c(x %*% coef) + offset - sum(coef*AIPW_fit$means)

    names(coef) <- dimnames(x)[[2]]

	  # Sort the data
	    if(length(strata) == 0){
	        ord <- order(y[,ny-1], -status)
	        newstrat <- as.integer(rep(0, n))
	    }else{
	        ord <- order(strata, y[,ny-1], -status)
	        newstrat <- c(diff(as.numeric(strata[ord]))!=0 ,1)
	    }
	    newstrat[n] <- 1

	    sy <- as.matrix(y[ord,])
	    sx <- as.matrix(x[ord,])
      sid <- id[ord]
	    soffset <- offset[ord]
	    score <- exp( lp[ord])
	    sR <- R[ord]
	    smarker <- as.matrix(marker[ord,])
	    smarker[is.na(smarker)] <- -999


	    resid <- AIPW_agscore_cpp(as.double(sy[,1]),
	        as.double(sy[,2]),
	        as.double(sy[,3]),
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
		      first_cont_rr,
	        second_cont_bl,
	        second_cont_rr)

	    rr <- resid
	    rr[ord,] <- rr
	    dimnames(rr) =list(names(rownames), names(coef))
	    rr <- drop(rowsum(rr, collapse))

	    temp <- 0*coef
	    score <- exp( sx %*% temp)
	    resid0 <- AIPW_agscore_cpp(as.double(sy[,1]),
	                              as.double(sy[,2]),
	                              as.double(sy[,3]),
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
	                              first_cont_rr,
	                              second_cont_bl,
	                              second_cont_rr)

	    rr0 <- resid0
	    rr0[ord,] <- rr0
	    dimnames(rr0) =list(names(rownames), names(coef))
	    rr0 <- drop(rowsum(rr0, collapse))


    afit <- list(coefficients  = coef,
    naive.var    = var,
    Ithegam = AIPW_fit$Ithegam,
    Ithealp = AIPW_fit$Ithealp,
    loglik = AIPW_fit$loglik,
    sctest = AIPW_fit$sctest,
    iter   = AIPW_fit$iter,
    conv   = AIPW_fit$conv,
    linear.predictors = as.vector(lp),
    resid = rr,
    score0 = rr0,
	  method='AIPW')

    return(afit)

}


## Modification of Therneau T (2015). _A Package for Survival Analysis in S_. version 2.38, <URL: https://CRAN.R-project.org/package=survival>.*/

#' Predicting the cumulative incidence function
#'
#' Predicting the cumulative incidence function for \code{subtype}, \code{IPWsubtype}, and \code{AIPWsubtype} models.
#'
#' @param fit \code{subtype}, \code{IPWsubtype} or \code{AIPWsubtype} object.
#' @param newdata a data.frame which has variables in formula and markers. This value must be specified.
#' @param individual a logical value: if \code{FALSE}, each row of newdata is considered as a different subject. If \code{TRUE}, \code{id} must be specified indicating variable name of subjects and time interval must be included in \code{newdata} (start, stop, status). Multiple rows represent different time points for one subject.
#' @param id a variable name identifying subject id.
#' @param na.action na.action to be passed
#'
#' @details The cumulative incidence function with a subtype defined by marker values in \code{newdata}. If a fitted object has strata, strata variable must be included in \code{newdata}.
#' @return An object of class \code{cif} representing the fit. A list of time and cumulative incidence functions for distinct subjects.
#'

#'
#' @export
cif <- function(fit, newdata, individual = FALSE, id, na.action = na.pass){

  Call <- match.call()

  if (is.null(fit)){
    cat("object must be either AIPWsubtype, IPWsubtype or subtype.")
    return(NULL)
  }

  if(is.null(newdata)){
    stop("newdata must be specified.")
  }

  missid <- missing(id) # I need this later, and setting id below makes
  # "missing(id)" always false
  if (!missid) individual <- TRUE
  else if (missid && individual) id <- rep(0,nrow(newdata))  #dummy value
  else id <- NULL


  n_marker <- length(fit$subtype$marker_name)
  cumhaz <- lapply(1:length(fit$basehaz), function(i) cbind(fit$basehaz[[i]][,1], cumsum(fit$basehaz[[i]][,2])))

  event <- tail(survival:::terms.inner(fit$formula[1:2]), 1)
  nvar <- length(fit$coefficients)
  n_subtype <- fit$subtype$n_subtype
  total_subtype <- fit$subtype$total_subtype

  marker <- newdata[, fit$subtype$marker_name]
  if(any(marker==0)) stop("Newdata should not have misisng markers")


  n <- nrow(newdata)
  cause <- rep(0, n)
  for (i in 1:n) {
    for (j in 1:fit$subtype$n_subtype) {
      if (all(marker[i, ] == fit$subtype$total_subtype[j, 1:n_marker]))
        cause[i] <- j
      }
    }


  newdata <- cbind(newdata, cause = cause)
  ## Data duplication
  lf <- function(x) {
    if (!is.na(x$cause[1])) {
      x$cause <- c(x$cause[1], seq(1, fit$subtype$n_subtype)[!(seq(1, fit$subtype$n_subtype) %in% x$cause[1])])
    }
    x
  }

  if(!individual) newdata$id <- seq(1, n)
  dnewdata <- lapply(1:n, function(i) newdata[rep(i, each = fit$subtype$n_subtype), ])
  dnewdata <- lapply(dnewdata, lf)
  dnewdata <- as.data.frame(rbindlist(dnewdata))
  dnewdata[, fit$subtype$marker_name] <- data.frame(fit$subtype$total_subtype[dnewdata$cause, ])

  Terms <- fit$terms
  has.strata <- !is.null(fit$strata)
  if (has.strata){
    stangle <- untangle.specials(Terms, 'strata')
    strata <- fit$strata
  }

  subterms <- function(tt, i) {
    dataClasses <- attr(tt, "dataClasses")
    predvars <- attr(tt, "predvars")
    oldnames <-  dimnames(attr(tt, 'factors'))[[1]]
    tt <- tt[i]
    index <- match(dimnames(attr(tt, 'factors'))[[1]], oldnames)
    if (length(index) >0) {
      if (!is.null(predvars))
        attr(tt, "predvars") <- predvars[c(1, index+1)]
      if (!is.null(dataClasses))
        attr(tt, "dataClasses") <- dataClasses[index]
    }
    tt
  }

  temp <- untangle.specials(Terms, 'cluster')
  if (length(temp$terms))
    Terms <- subterms(Terms, -temp$terms)

  Terms2 <- Terms
  if (!individual)  Terms2 <- delete.response(Terms)

  if (is.null(names(newdata))) {
    stop("Newdata argument must be a data frame")
  }

  if(missid){
  if(has.strata){
    found.strata <- TRUE
    tempenv <- new.env(, parent=emptyenv())
    assign("strata", function(..., na.group, shortlabel, sep)
      list(...), envir=tempenv)
    assign("list", list, envir=tempenv)
    for (svar in stangle$vars) {
      temp <- try(eval(parse(text=svar), dnewdata, tempenv),
                silent=TRUE)
      if (!is.list(temp) ||
        any(unlist(lapply(temp, class))== "function"))
      found.strata <- FALSE
    }
    if (!found.strata) stop("Newdata should not be missing with strata")
    if (found.strata) mf2 <- stats::model.frame(Terms2, data=dnewdata,
                                            na.action=na.action, xlev=fit$xlevels)
  }else{
    mf2 <- stats::model.frame(Terms2, data=dnewdata, na.action=na.action,
                              xlev=fit$xlevels)
    found.strata <- has.strata
  }
  }else{
    tcall <- Call[c(1, match(c('id', "na.action"),
                             names(Call), nomatch=0))]

    tcall$data <- dnewdata
    tcall$formula <- Terms2
    tcall$xlev <- fit$xlevels
    tcall[[1L]] <- quote(stats::model.frame)
    mf2 <- eval(tcall)
    found.strata <- has.strata # would have failed otherwise
  }


  if (has.strata && found.strata) { #pull them off
    temp <- untangle.specials(Terms2, 'strata')
    strata2 <- strata(mf2[temp$vars], shortlabel=TRUE)
    strata2 <- factor(strata2, levels=levels(strata))
    if (any(is.na(strata2)))
      stop("New data set has strata levels not found in the original")
    # An expression like age:strata(sex) will have temp$vars= "strata(sex)"
    #  and temp$terms = integer(0).  This does not work as a subscript
    if (length(temp$terms) >0) Terms2 <- Terms2[-temp$terms]
  }else{
    strata2 <- factor(rep(1, n))
  }
  Strata <- levels(strata2)
  ns <- length(Strata)


  if (individual) {
   if (!missid) {  #grab the id variable
      id <- model.extract(mf2, "id")
      if (is.null(id)) stop("id=NULL is an invalid argument")
    }
    else id <- rep(1, nrow(mf2))

    x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
    if (length(x2)==0) stop("Individual survival but no variables")
    offset2 <- model.offset(mf2)
    if (length(offset2) == 0) offset2 <- 0

    y2 <- model.extract(mf2, 'response')
     if (attr(y2, "type") != "counting")
      stop("Individual=TRUE is only valid for counting process data")
    y2 <- y2[,1:2, drop=F]  #throw away status, it's never used

    newrisk <- exp(c(x2 %*% fit$coefficients) + offset2)

    onecurve <- function(basehaz, cumhaz, x2, y2, strata2,  Strata, newrisk, n_subtype) {
    ntarget <- nrow(x2)/n_subtype
    surv <- newcumhaz <- newhaz <- time <- vector('list', ntarget)
    whichstr <- match(strata2, Strata)
    timeforward <- 0
     i=1; j=1
    while(i <nrow(x2)){
      tt <- cumhaz[[whichstr[i]]][,1]
      indx <- which(tt > y2[i,1] & tt <= y2[i,2])
      if (length(indx)==0) {
        timeforward <- timeforward + y2[i,2] - y2[i,1]
        # No deaths or censors in user interval.  Possible
        # user error, but not uncommon at the tail of the curve.
      }
      else {
        time[[j]] <- diff(c(y2[i,1], cumhaz[[whichstr[i]]][indx,1])) #time increments
        time[[j]][1] <- time[[j]][1] + timeforward
        timeforward <- y2[i,2] - max(cumhaz[[whichstr[i]]][indx,1])
        newcumhaz <- outer(cumhaz[[whichstr[i]]][indx,2], newrisk[i:(i+n_subtype-1)], '*')
        newhaz[[j]] <- fit$basehaz[[whichstr[i]]][indx,2]*newrisk[i]
        surv[[j]] <- exp(-drop(rowSums(newcumhaz)))
      }
      i = i + fit$subtype$n_subtype
      j = j +1
    }
    res <- list(time = cumsum(unlist(time)), cif = cumsum(unlist(newhaz)*unlist(surv)))
    return(res)
    }
    if (all(id ==id[1])) {
      result <- list(onecurve(basehaz, cumhaz, x2, y2, strata2, Strata, newrisk, fit$subtype$n_subtype))
    }
    else {
      uid <- unique(id)
      result <- vector('list', length=length(uid))
      for (i in 1:length(uid)) {
        indx <- which(id==uid[i])
        result[[i]] <- onecurve(basehaz, cumhaz,x2[indx,,drop=FALSE],
                                y2[indx,,drop=FALSE],
                                strata2[indx],  Strata, newrisk[indx], fit$subtype$n_subtype)
      }
      names(result) <- uid
    }
  }else{
  whichstr <- match(strata2, Strata)
  offset2 <- model.offset(mf2)
  if (length(offset2) == 0)  offset2 <- 0
  x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
  risk <- exp(c(x2 %*% fit$coefficients) + offset2)
  ntarget <- nrow(x2)/n_subtype
  newrisk <- split(risk, newdata$id)

  result <- list()
  tmp <- list()
  for(i in 1:ntarget){
    temp1 <- cumhaz[[whichstr[i]]]
    temp2 <- fit$basehaz[[whichstr[i]]]
    newcumhaz <- outer(temp1[,2], newrisk[[i]], '*')
    newhaz <- temp2[,2]*newrisk[[i]][1]
    surv <- exp(-drop(rowSums(newcumhaz)))
    tmp$time <- temp1[,1]
    tmp$cif <- cumsum(newhaz*surv)
    result[[i]] <- tmp
  }

  }


  class(result) <- "cif"
  return(result)

  }







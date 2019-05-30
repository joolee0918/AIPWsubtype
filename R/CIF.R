
#' @export
CIF <- function(fit, newdata, na.action = na.pass){

  if (is.null(fit)){
    cat("object must be either AIPWsubtype, IPWsubtype or subtype.")
    return(NULL)
  }

  n_marker <- length(fit$subtype$marker_name)
  cumhaz <- lapply(1:length(fit$basehaz), function(i) cbind(fit$basehaz[[i]][,1], cumsum(fit$basehaz[[i]][,2])))

  event <- tail(survival:::terms.inner(fit$formula[1:2]), 1)
  nvar <- length(fit$coefficients)
  n_subtype <- fit$subtype$n_subtype
  total_subtype <- fit$subtype$total_subtype

  if(missing(newdata)){
    newdata <- data.frame(t(rep(0, fit$subtype$nX + fit$subtype$nW))) #create a dummy newdata
    names(newdata) <- names(fit$coefficients[1:(fit$subtype$nX + fit$subtype$nW)])
    found.strata <- FALSE
    marker <- data.frame(t(rep(1, n_marker)))
    names(marker) <- fit$subtype$marker_name
    newdata <- cbind(newdata, marker)
  } else{
    marker <- newdata[, fit$subtype$marker_name]
    if(any(marker==0)) stop("Newdata should not have misisng markers")
  }

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

  newdata$id <- seq(1, n)
  dnewdata <- lapply(1:n, function(i) newdata[rep(i, each = fit$subtype$n_subtype), ])
  dnewdata <- lapply(dnewdata, lf)
  dnewdata <- as.data.frame(rbindlist(dnewdata))
  dnewdata[, marker_name] <- data.frame(fit$subtype$total_subtype[dnewdata$cause, ])

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
  Terms2 <- delete.response(Terms)

  if (is.null(names(newdata))) {
    stop("Newdata argument must be a data frame")
  }

  if(has.strata){
    found.strata <- TRUE
    tempenv <- new.env(, parent=emptyenv())
    assign("strata", function(..., na.group, shortlabel, sep)
      list(...), envir=tempenv)
    assign("list", list, envir=tempenv)
    for (svar in stangle$vars) {
      temp <- try(eval(parse(text=svar), newdata, tempenv),
                silent=TRUE)
      if (!is.list(temp) ||
        any(unlist(lapply(temp, class))== "function"))
      found.strata <- FALSE
    }
    if (!found.strata) stop("Newdata should not be missing with strata")
    if (found.strata) mf2 <- stats::model.frame(Terms2, data=newdata,
                                            na.action=na.action, xlev=fit$xlevels)
  }else{
    found.strata <- has.strata
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
  ntarget <- length(strata2)
  whichstr <- match(strata2, Strata)

  mf2 <- stats::model.frame(Terms2, data=dnewdata, na.action=na.action, xlev=fit$xlevels)

  offset2 <- model.offset(mf2)
  if (length(offset2) == 0)  offset2 <- 0
  x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept

  risk <- exp(c(x2 %*% fit$coefficients) + offset2)
  newrisk <- split(risk, newdata$id)

  pt <- list()
  for(i in 1:ntarget){
    temp1 <- cumhaz[[whichstr[i]]]
    temp2 <- fit$basehaz[[whichstr[i]]]
    newcumhaz <- outer(temp1[,2], newrisk[[i]], '*')
    newhaz <- temp2[,2]*newrisk[[i]][1]
    surv <- exp(-drop(rowSums(newcumhaz)))
    pt[[i]] <- cbind(temp1[,1], cumsum(newhaz*surv))
    names(pt[[i]]) <- c("time", "CIF")
  }


  class(pt) <- "CIF"
  return(pt)

  }







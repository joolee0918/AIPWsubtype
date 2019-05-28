
### Generate Competing Risk

# Weibull distribution -> generate time to event


sim_cpt <- function(ntype, cshaz, nid) {
    
    sim.ind <- data.frame()
    pro <- vector()
    
    X = rbinom(1, 1, 0.5)
    C = rnorm(1, 75, 5)
    W = rbinom(1, 1, 0.3)
    A <- function(t, y) {
        res <- 0
        for (k in 1:ntype) {
            res <- res + integrate(cshaz[[k]], lower = 0.001, upper = t, r = k, x = X, w = W, subdivisions = 1000)$value
        }
        res <- res + y
        return(res)
    }
    u <- runif(1)
    iters <- 0
    while (A(0.001, log(1 - u)) * A(200, log(1 - u)) > 0 & iters < 1000) {
        u <- runif(1)
        iters <- iters + 1
    }
    if (iters >= 1000) 
        stop("Error: Values at endpoints not of opposite sign. \n")
    
    tb <- uniroot(A, c(0, 200), tol = 1e-04, y = log(1 - u))$root
    
    sumprob <- 0
    for (k in 1:length(cshaz)) {
        sumprob <- sumprob + cshaz[[k]](t = tb, r = k, x = X, w = W)
    }
    
    for (k in 1:length(cshaz)) {
        pro[k] <- cshaz[[k]](t = tb, r = k, x = X, w = W)/sumprob
    }
    
    cause1 <- rmultinom(1, 1, prob = pro)
    for (k in 1:length(cshaz)) {
        if (cause1[k] == 1) 
            cause <- k
    }
    
    y1 <- ifelse(cause == 1 | cause == 3, 1, 0)
    y2 <- ifelse(cause == 1 | cause == 2, 1, 0)
    
    time <- ifelse(tb < C, tb, C)
    del <- ifelse(tb < C, 1, 0)
    
    sim.ind <- data.frame(nid = nid, cause = cause, time = time, status = del, X = X, W = W)
    return(sim.ind)
}


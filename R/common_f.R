
exposit <- function(x) {
    exp(x)/(1 + exp(x))
}
dexposit <- function(x) {
    exp(x)/(1 + exp(x))^2
}

pi = function(w, alp, r) {
    walp = w %*% as.matrix(alp)
    res = exp(walp)/(1 + exp(walp))
    if (r == 0)
        res = 1 - res
    return(res)
}


dpi = function(w, alp, r) {
    walp = w %*% as.matrix(alp)
    ewalp = c(exp(walp))
    res = w * ewalp/(1 + ewalp)^2
    if (r != 1)
        res = -res
    return(res)
}


terms.inner <- function(x) {
    if (inherits(x, "formula")) {
        if (length(x) ==3) c(terms.inner(x[[2]]), terms.inner(x[[3]]))
        else terms.inner(x[[2]])
    }
    else if (inherits(x, "call") && 
             (x[[1]] != as.name("$") && x[[1]] != as.name("["))) {
        if (x[[1]] == '+' || x[[1]]== '*' || x[[1]] == '-') {
            # terms in a model equation, unary minus only has one argument
            if (length(x)==3) c(terms.inner(x[[2]]), terms.inner(x[[3]]))
            else terms.inner(x[[2]])
        }
        else if (x[[1]] == as.name("Surv"))
                 unlist(lapply(x[-1], terms.inner))
        else terms.inner(x[[2]])
    }
    else(deparse(x))
}


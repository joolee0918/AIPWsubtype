
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

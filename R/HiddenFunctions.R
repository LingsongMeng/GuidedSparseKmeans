# Hidden functions in GuidedSparseKmeans function and GuidedSparseKmeans.R2out function

l1n <- function(vec) sum(abs(vec))


l2n <- function(vec) sqrt(sum(vec^2))


update.w <- function(x, c, s, lam, R2.per) {
    SS <- getSS(x, c)
    obj.per <- SS$bcss.per/SS$tss.per + lam * R2.per
    delta <- SearchDelta(obj.per, s)
    w <- soft(obj.per, delta)/l2n(soft(obj.per, delta))
    R2.w <- sum(R2.per * w)
    bc.t.prop.w <- sum(SS$bcss.per/SS$tss.per * w)
    return(list(w = w, R2.w = R2.w, bc.t.prop.w = bc.t.prop.w))
}


getSS <- function(x, c, w = NULL) {
    wcss.per <- rep(0, ncol(x))
    for (j in unique(c)) {
        ind <- which(c == j)
        if (length(ind) > 1) 
            wcss.per <- wcss.per + apply(scale(x[ind, ], center = T, scale = F)^2, 2, sum) else wcss.per <- wcss.per + scale(x[ind, ], center = T, scale = F)^2
    }
    tss.per <- apply(scale(x, center = T, scale = F)^2, 2, sum)
    bcss.per <- tss.per - wcss.per
    if (is.null(w)) 
        return(list(wcss.per = wcss.per, tss.per = tss.per, bcss.per = bcss.per, wcss = sum(wcss.per), tss = sum(tss.per), 
            bcss = sum(bcss.per))) else return(list(wcss.per = wcss.per, tss.per = tss.per, bcss.per = bcss.per, wcss = sum(wcss.per), tss = sum(tss.per), 
        bcss = sum(bcss.per), wcss.w = sum(wcss.per * w), bcss.w = sum(bcss.per * w), bc.t.prop.w = sum(bcss.per/tss.per * 
            w)))
}


SearchDelta <- function(obj.per, s) {
    delta1 <- 0
    delta2 <- max(abs(obj.per))
    while (delta2 - delta1 > 1e-04) {
        w <- soft(obj.per, (delta1 + delta2)/2)/l2n(soft(obj.per, (delta1 + delta2)/2))
        if (l1n(w) < s) 
            delta2 <- (delta1 + delta2)/2 else delta1 <- (delta1 + delta2)/2
    }
    delta <- (delta1 + delta2)/2
    return(delta)
}


soft <- function(obj.per, delta) sign(obj.per) * pmax(0, abs(obj.per) - delta)

##' Guided Sparse K-means (R-square or pseudo R-square is from the outside)
##'
##' Guided Sparse K-means integrating clinical dataset with gene expression dataset,
##' R-square or pseudo R-square is from the outside, not calculated in the function.
##' @title GuidedSparseKmeans.R2out
##' @param x Gene expression matrix, n*p (rows for subjects and columns for genes).
##' @param R2.per R-squared or pseudo R-squared between phenotypic variable and expression value of each gene, a vector.
##' @param K Number of clusters.
##' @param s The boundary of l1n weights, a vector.
##' @param lam The intensity of guidance.
##' @param nstart Specify the number of starting point for K-means.
##' @param maxiter Maximum number of iteration.
##' @param silence Output progress or not.
##' @return m lists, m is the length of parameter s. Each list is consisting of
##' \item{weights}{weight for each feature, zero weight means the feature is not selected.}
##' \item{clusters}{cluster results.}
##' \item{object}{objective value.}
##' \item{bound}{a boundary of l1n weights}
##' @export
##' @author Lingsong Meng


GuidedSparseKmeans.R2out <- function(x, R2.per, K, s, lam, nstart = nstart, maxiter = 15, silence = F) {
    
    if (is.null(K)) 
        stop("Must provide K.")
    if (is.null(s)) 
        stop("Must provide s.")
    if (min(s) <= 1) 
        stop("s should be greater than 1.")
    if (is.null(lam)) 
        stop("Must provide lam.")
    
    ind.init.zero <- order(R2.per, decreasing = T)[401:(ncol(x))]
    R2.per.new <- R2.per
    R2.per.new[ind.init.zero] <- 0
    R2.per.prop <- abs(R2.per.new/sum(R2.per.new))
    out <- list()
    for (i in 1:length(s)) {
        if (silence == F) 
            print(paste0("s=", s[i]))
        w <- R2.per.prop * s[i]
        w.old <- rep(1, ncol(x))
        obj.w <- NULL
        niter <- 0
        while (l1n(w - w.old)/l1n(w.old) > 1e-04 && niter < maxiter) {
            niter <- niter + 1
            if (silence == F) 
                print(niter)
            w.old <- w
            y <- sweep(x[, w != 0], 2, sqrt(w[w != 0]), "*")
            c <- kmeans(y, K, nstart = nstart)$cluster
            update <- update.w(x, c, s[i], lam, R2.per)
            w <- update$w
            obj.w <- c(obj.w, update$bc.t.prop.w + lam * update$R2.w)
        }
        out[[i]] <- list(weights = w, clusters = c, object = obj.w, bound = s[i])
    }
    return(out)
}

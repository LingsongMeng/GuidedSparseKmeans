##' Selection of Tuning Parameter s in Guided Sparse K-means
##'
##' Select tuning parameter s via permutation in Guided Sparse K-means integrating clinical dataset with gene expression dataset.
##' @title GuidedSparseKmeans.S
##' @param x Gene expression matrix, n*p (rows for subjects and columns for genes).
##' @param R2.per R-squared or pseudo R-squared between phenotypic variable and expression value of each gene, a vector.
##' @param K Number of clusters.
##' @param s The boundary of l1n weights, a vector.
##' @param lam The intensity of guidance.
##' @param nstart Specify the number of starting point for K-means.
##' @param maxiter Maximum number of iteration.
##' @param nperms Number of permutation times
##' @param silence Output progress or not.
##' @return A list consisting of
##' \item{nnonzerows}{number of nonzero weights.}
##' \item{gaps}{gap statistics.}
##' \item{segaps}{statard error of gap statistics.}
##' \item{s}{candidates of s.}
##' \item{s.best}{the best s with largest gap.}
##' @export
##' @author Lingsong Meng


GuidedSparseKmeans.S.R2out <- function(x, R2.per, K, s, lam, nstart = 20, maxiter = 15, nperms = 50, silence = F) {
    
    if (is.null(s)) 
        stop("Must provide either K or centers.")
    if (min(s) <= 1) 
        stop("s should be greater than 1, since otherwise only one weight will be nonzero.")
    if (length(s) < 2) 
        stop("s should be a vector of at least two elements.")
    if (is.null(K)) 
        stop("Must provide either K or centers.")
    
    permx <- list()
    nnonzerows <- NULL
    for (i in 1:nperms) {
        permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
        for (j in 1:ncol(x)) permx[[i]][, j] <- sample(x[, j])
    }
    tots <- NULL
    out <- GuidedSparseKmeans.R2out(x, R2.per, K, s, lam, nstart, maxiter, silence)
    for (i in 1:length(out)) {
        nnonzerows <- c(nnonzerows, sum(out[[i]]$weights != 0))
        tots <- c(tots, getSS(x = x, c = out[[i]]$clusters, w = out[[i]]$weights)$bcss.w)
    }
    permtots <- matrix(NA, nrow = length(s), ncol = nperms)
    for (k in 1:nperms) {
        if (silence == F) 
            print(paste0("nperms=", k))
        perm.out <- GuidedSparseKmeans.R2out(permx[[k]], R2.per, K, s, lam, nstart, maxiter, silence)
        for (i in 1:length(perm.out)) {
            permtots[i, k] <- getSS(x = permx[[k]], c = perm.out[[i]]$clusters, w = perm.out[[i]]$weights)$bcss.w
        }
    }
    gaps <- (log(tots) - apply(log(permtots), 1, mean))
    segaps <- apply(log(permtots), 1, sd)/sqrt(nperms)
    
    out <- list(nnonzerows = nnonzerows, gaps = gaps, segaps = segaps, s = s, s.best = s[which.max(gaps)])
    return(out)
}

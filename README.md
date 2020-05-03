# GuidedSparseKmeans
Github repository for Outcome-guided Sparse K-means method (GuidedSparseKmeans)


## Install This Package from Github
In R console

```{R}
library(devtools)
install_github("LingsongMeng/GuidedSparseKmeans") 
```


## Short Tutorial
This short tutorial is about how to apply GuidedSparseKmeans to obtain subtype clustering result and gene selection result.

```{R}
library(GuidedSparseKmeans)
n <- 100 ## number of subjects
g <- 5000 ## number of gene features
set.seed(123)
epxression <- matrix(rnorm(n*g), nrow=g, ncol=n) ## briefly simulate gene expression data
outcome <- rnorm(n) ## simulate a continuous outcome variable


# Estimate tuning parameter K and lambda
KLam <- GuidedSparseKmeans.KLam(x=t(epxression), z=outcome, pre.K = NULL, s.one=10, model="linear", nstart = 20, maxiter = 15, silence = F)

# Estimate tuning parameter s
s <- GuidedSparseKmeans.S(x=t(epxression), z=outcome, K=3, s=c(8:12), lam=1, model="linear", nstart = 20, maxiter = 15, nperms = 50, silence = F)

# Obtain gene selection and clustering results
results <- GuidedSparseKmeans(x=t(epxression), outcome, K=3, s=10, lam=1, model="linear", nstart=20, maxiter=15)


## Users can also obtain guidance first, then use the guidance in the next several steps.
# Obtain guidance 
R2 <- getR2(t(epxression), outcome, model="linear")

# Estimate tuning parameter K and lambda
KLam <- GuidedSparseKmeans.KLam.R2out(x=t(epxression), R2.per=R2, pre.K = NULL, s.one=10, nstart = 20, maxiter = 15, silence = F)

# Estimate tuning parameter s
s <- GuidedSparseKmeans.S.R2out(x=t(epxression), R2.per=R2, K=3, s=c(8:12), lam=1, nstart = 20, maxiter = 15, nperms = 50, silence = F)

# Obtain gene selection and clustering results
results <- GuidedSparseKmeans.R2out(x=t(epxression), R2.per=R2, K=3, s=10, lam=1, nstart=20, maxiter=15)

```

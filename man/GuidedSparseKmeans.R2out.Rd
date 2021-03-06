% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GuidedSparseKmeans.R2out.R
\name{GuidedSparseKmeans.R2out}
\alias{GuidedSparseKmeans.R2out}
\title{GuidedSparseKmeans.R2out}
\usage{
GuidedSparseKmeans.R2out(x, R2.per, K, s, lam, nstart = nstart,
  maxiter = 15, silence = F)
}
\arguments{
\item{x}{Gene expression matrix, n*p (rows for subjects and columns for genes).}

\item{R2.per}{R-squared or pseudo R-squared between phenotypic variable and expression value of each gene, a vector.}

\item{K}{Number of clusters.}

\item{s}{The boundary of l1n weights, a vector.}

\item{lam}{The intensity of guidance.}

\item{nstart}{Specify the number of starting point for K-means.}

\item{maxiter}{Maximum number of iteration.}

\item{silence}{Output progress or not.}
}
\value{
m lists, m is the length of parameter s. Each list is consisting of
\item{weights}{weight for each feature, zero weight means the feature is not selected.}
\item{clusters}{cluster results.}
\item{object}{objective value.}
\item{bound}{a boundary of l1n weights}
}
\description{
Guided Sparse K-means (R-square or pseudo R-square is from the outside)
}
\details{
Guided Sparse K-means integrating clinical dataset with gene expression dataset,
R-square or pseudo R-square is from the outside, not calculated in the function.
}
\author{
Lingsong Meng
}

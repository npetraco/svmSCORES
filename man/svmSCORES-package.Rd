\name{svmSCORES-package}
\alias{svmSCORES-package}
\alias{svmSCORES}
\docType{package}
\title{
Compute a vector of Platt scores with a bootstrap process. Intended for use with fdrID package.
}
\description{
The fdrID package requires null (non-match) and non-null (match) scores as input for estimating fdr. This package provides a set of scores derived from an SVM (Platt scores) which can be used as the requisite input. The scores are obtained from a group-wise bootstrapping process.
\details{
\tabular{ll}{
Package: \tab svmSCORES\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-26\cr
License: \tab GPL (>= 2)\cr
}
The fdrID package requires null (non-match) and non-null (match) scores as input for estimating fdr. This package provides a set of scores derived from an SVM (Platt scores: Advances in large margin classifiers 10(3):61–74 (1999)) which can be used as the requisite input. The scores are obtained from a group-wise bootstrapping process which is designed to minimize the amount of dependence between the scores. The bootstrapping process is similar to the permutation method used to generate an estimate of the null score distribution in Storey and Tibshirani: Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003).}
\author{
Nick Petraco

Maintainer: <npetraco@gmail.com>
}
\references{
Platt J. Probabilistic outputs for support vector machines and comparisons to regularized likelihood methods. Advances in large margin classifiers 10(3):61–74 (1999)

Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
}
\keyword{ package }
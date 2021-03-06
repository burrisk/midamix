% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute.R
\name{impute}
\alias{impute}
\title{Multiple imputation for data analysis via mixtures}
\usage{
impute(data, imputations = 10, max_clusters = 15, n_iter = 1000,
  burnin = 100, validator = NULL, cap = NULL, seed = 1)
}
\arguments{
\item{data}{A data frame, consisting of numeric, integer,
or ordered factor columns.}

\item{imputations}{The number of imputed datasets to return, defaults to 10.}

\item{max_clusters}{The maximum number of clusters for the mixture model.}

\item{n_iter}{Number of iterations for the MCMC sampler.}

\item{burnin}{Number of iterations for initial burn-in period.}

\item{validator}{A function that takes in an observation and determines whether
it is feasible.}

\item{cap}{Maximum number of rejected proposals allowed for the constrained sampler}

\item{seed}{Random seed.}
}
\value{
A \code{tibble} consisting of multiply imputed data sets.
}
\description{
Uses a Dirichlet Process Mixture Transformation Model to sample missing data values
from the corresponding posterior distribution.
}

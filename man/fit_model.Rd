% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pool_inferences.R
\name{fit_model}
\alias{fit_model}
\title{Fit multiple models}
\usage{
fit_model(imputations, model_function)
}
\arguments{
\item{imputations}{A nested tibble, representing multiply imputed datasets.}

\item{model_function}{A function that takes a dataset as input and returns a fitted
model object as output.}
}
\value{
A tibble with the fits of each model.
}
\description{
Fits a user-specified model to each imputed dataset.
}

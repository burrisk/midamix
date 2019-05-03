
<!-- README.md is generated from README.Rmd. Please edit that file -->

# midamix <a href='https://github.com/burrisk/midamix'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/burrisk/midamix.svg?branch=master)](https://travis-ci.org/burrisk/midamix)
<!-- badges: end -->

## Overview

midamix is an R package dedicated to flexible, tidy handling of missing
data. For numeric and/or ordinal datasets, midamix uses a transformation
model coupled with a truncated Dirichlet process mixture model to
generate multiple completed datasets. It contains the following methods
for performing principled, little-hassle, multiple imputation.

  - `impute()` generates multiply imputed datasets
  - `fit_models()` fits a model to each of the multiply imputed
    datasets.
  - `pool_inferences()` combines inferences for linear models or
    generalized linear models fit on each of the datasets via the
    combining rules of Rubin (1987).

### Development version

To use the development version, you can install `midamix` from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("burrisk/midamix")
```

## Usage

``` r
library(midamix)
```

## Getting help

If you encounter a clear bug, please file a minimal reproducible example
on [github](https://github.com/burrisk/midamix/issues).

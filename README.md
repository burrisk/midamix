<!-- badges: start -->
[![Travis build status](https://travis-ci.org/burrisk/midamix.svg?branch=master)](https://travis-ci.org/burrisk/midamix)
<!-- badges: end -->

## Overview

midamix is an R package dedicated to flexible, tidy handling of missing data.  For numeric and/or ordinal datasets, midamix uses a transformation model coupled with a truncated Dirichlet process mixture model to generate multiple completed datasets.   C

  - `impute()` generates multiply imputed datasets
  - `fit_models()` fits a model to each of the multiply imputed datasets.
  - `pool_inferences()` combines inferences on each of the datasets via the combining rules of Rubin (1987).

### Development version

To use the development version, you can
install `midamix` from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("burrisk/midamix")
```

## Usage

``` r
library(midamix)

data(diamonds)

imputations <- impute(diamonds) %>%
  print()
```

## Getting help

If you encounter a clear bug, please file a minimal reproducible example
on [github](https://github.com/burrisk/midamix/issues).

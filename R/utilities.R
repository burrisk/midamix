#' Multivariate normal distribution
#'
#' Random number generator and density for the multivariate normal distribution with mean equal
#' to \code{mean} and covariance matrix \code{sigma}
#'
#' @param n The number of observations
#' @param mean Mean vector, default is 0
#' @param sigma Covariance matrix, default is identity.
#' @param prec Precision matrix, default is identity
#' @param method string specifying the matrix decomposition used to determine the
#' matrix root of sigma
#' @param x A vector or matrix of quantiles.  If x is a matrix, each row is taken to
#' be a quantile
#' @param log logical; if \code{TRUE}, log-densities are returned
#' @param pre0.9_9994 logical; version of rmvnorm used, inherited from mvtnorm package
rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                     method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) {
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma))
    stop("mean and sigma have non-conforming size")
  method <- match.arg(method)
  R <- if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,
                                                0))))
  }
  else if (method == "svd") {
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  }
  else if (method == "chol") {
    R <- chol(sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%
    R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}

#' @rdname rmvnorm
dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) {
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
    if (length(mean) != p)
      stop("mean and sigma have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE))
      stop("sigma must be a symmetric matrix")
  }
  dec <- tryCatch(chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 *
                                                        pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log)
    logretval
  else exp(logretval)
}

#' @rdname rmvnorm
dmvnorm_prec <- function(x, mean = rep(0, p), prec = diag(p), log = FALSE){
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
  }
  if (!missing(prec)) {
    if (p != ncol(prec))
      stop("x and prec have non-conforming size")
  }
  logretval <- -0.5 * p * log(2 * pi) + 0.5 * log(abs(det(prec))) -
    0.5 * diag(t(t(x) - mean) %*% prec %*% (t(x) - mean))
  if (log)
    logretval
  else exp(logretval)
}


#' Random Wishart distributions
#'
#' Randomly generates positive definite matrices from Wishart and Inverse-Wishart distributions
#'
#' @param v degrees of freedom
#' @param S inverse scale matrix
#' @return A positive definite matrix with the same dimension as S
rwish <- function (v, S) {
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p *
                                                                  (p - 1)/2)
  }
  return(crossprod(Z %*% CC))
}

#' @rdname rwish
riwish <- function (v, S) {
  return(solve(rwish(v, solve(S))))
}


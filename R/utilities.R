#' Random Wishart distributions
#'
#' Randomly generates positive definite matrices from Wishart and Inverse-Wishart distributions
#'
#' @param v degrees of freedom
#' @param S inverse scale matrix
#' @return A positive definite matrix with the same dimension as S
rwish <- function(v, S) {
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
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * (p - 1)/2)
    }
    return(crossprod(Z %*% CC))
}

#' @rdname rwish
riwish <- function(v, S) {
    return(solve(rwish(v, solve(S))))
}



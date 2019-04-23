update_alpha <- function(model_params, hyperpars){
  a_alpha <- hyperpars$a_alpha
  b_alpha <- hyperpars$b_alpha
  cluster_probs <- model_params$cluster_probs
  K <- length(cluster_probs)
  alpha <- rgamma(1, a_alpha + K - 1, b_alpha - log(cluster_probs[K]))
  alpha
}

update_v <- function(model_params, hyperpars){
  clusters <- model_params$clusters
  alpha <- model_params$alpha
  K <- length(model_params$cluster_probs)
  V <- rep(1, K)
  for (index in 1:(K-1)){
    in_cluster <- sum(clusters == index)
    later_cluster <- sum(clusters > index)
    V[index] <- rbeta(1, 1 + in_cluster, alpha + later_cluster)
  }
  V
}

calculate_cluster_probs <- function(model_params){
  V <- model_params$V
  K <- length(V)
  cluster_probs <- rep(0, K)
  cum_prod <- 1
  for (k in 1:K){
    cluster_probs[k] <- V[k] * cum_prod
    cum_prod <- cum_prod * (1 - V[k])
  }
  cluster_probs
}

update_sigma <- function(model_params, hyperpars){
  nu <- hyperpars$nu
  omega <- hyperpars$omega
  Sigma0 <- hyperpars$Sigma0
  cluster_covs <- model_params$cluster_covs
  cluster_precs <- cluster_covs
  p <- dim(cluster_cov)[2]
  K <- dim(cluster_cov)[3]
  for (i in 1:K){
    cluster_precs[, , i] <- solve(cluster_covs[, , i])
  }
  S <- solve(solve(Sigma0) + apply(cluster_precs, c(1, 2), sum))
  df <- nu * K + omega
  Sigma <- rwish(df, S)
  Sigma
}

update_tau <- function(model_params, hyperpars){
  a_tau <- hyperpars$a_tau
  b_tau <- hyperpars$b_tau
  B0 <- model_params$B0
  cluster_coefs <- model_params$cluster_coefs
  q_star <- nrow(B0)
  p <- ncol(B0)
  K <- dim(cluster_coefs)[3]
  tau <- rep(1, p)

  BmB0 <- sweep(cluster_coefs, c(1, 2), B0, "-")
  ss <- rowSums(apply(BmB0, c(3), function(X){
    colSums(X ^ 2)
  }))
  tau <- rgamma(p, a_tau + q_star * K / 2, b_tau + ss / 2)
  tau
}

update_b0 <- function(model_params, hyperpars){
  cluster_coefs <- model_params$cluster_coefs
  q <- dim(cluster_coefs)[1]
  p <- dim(cluster_coefs)[2]
  K <- dim(cluster_coefs)[3]
  sigma2b <- hyperpars$sigma2b
  tau <- model_params$tau
  sum_cluster <- apply(cluster_coefs, c(1, 2), sum)
  Omega <- diag(1 / tau)
  covar <- solve(Omega * K + diag(p) / sigma2b)
  B0 <- sum_cluster %*% Omega %*% covar + rmvnorm(q, rep(0, p), covar)
  B0
}

update_clusters <- function(model_params, hyperpars){
  cluster_probs <- model_params$cluster_probs
  cluster_coefs <- model_params$cluster_coefs
  cluster_covs <- model_params$cluster_covs
  Z <- model_params$Z
  X <- model_params$X
  n <- nrow(Z)
  p <- dim(cluster_coefs)[2]
  K <- dim(cluster_coefs)[3]
  clusters <- rep(0, n)
  for (i in 1:n){
    Xi <- matrix(X[i, ], nrow = 1)
    Zi <- Z[i, ]
    log_density <- rep(0, K)
    for (k in 1:K){
      log_density[k] <- dmvnorm(Zi, Xi %*% cluster_coefs[, , k], cluster_covs[, , k],
                                log = T)
    }
    log_probs <- log(cluster_probs) + log_density
    log_probs <- log_probs - max(log_probs)
    probs <- exp(log_probs) / sum(exp(log_probs))
    clusters[i] <- sample(1:K, 1, replace = F, prob = probs)
  }
  clusters
}

update_cluster_coefs <- function(model_params, hyperpars, conditionals){

}

update_cluster_covs <- function(model_params, hyperpars){

}

update_zobs <- function(model_params, hyperpars, conditionals){

}

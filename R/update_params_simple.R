#' @importFrom magrittr "%>%"

update_alpha <- function(model_params, hyperpars){
  log_cluster_probs <- model_params$log_cluster_probs
  K <- length(log_cluster_probs)
  model_params$alpha <- rgamma(1, hyperpars$a_alpha + K - 1,
                               hyperpars$b_alpha - log_cluster_probs[K])
  model_params
}

update_v <- function(model_params){
  alpha <- model_params$alpha
  clusters <- model_params$clusters
  K <- length(model_params$log_cluster_probs)
  V <- rep(1, K)
  for (k in 1:(K - 1)){
    in_cluster <- sum(clusters == k)
    later_cluster <- sum(clusters > k)
    V[k] <- min(1 - 1e-6, rbeta(1, 1 + in_cluster, alpha + later_cluster))
  }
  model_params$V <- V
  model_params
}

calculate_log_cluster_probs <- function(model_params){
  V <- model_params$V
  K <- length(V)
  log_cum_probs <- c(0, cumsum(log(1 - V[-K])))
  model_params$log_cluster_probs <- log(V) + log_cum_probs
  model_params
}

update_sigma <- function(model_params, hyperpars){
  nu <- hyperpars$nu
  omega <- hyperpars$omega
  Sigma0 <- hyperpars$Sigma0
  cluster_covs <- model_params$cluster_covs
  cluster_precs <- cluster_covs
  K <- dim(cluster_covs)[3]
  for (k in 1:K){
    cluster_precs[, , k] <- solve(cluster_covs[, , k])
  }
  S <- solve(solve(Sigma0) + apply(cluster_precs, c(1, 2), sum))
  df <- nu * K + omega
  model_params$Sigma <- rwish(df, S)
  model_params
}

update_clusters <- function(model_params){
  cluster_means <- model_params$cluster_means
  cluster_covs <- model_params$cluster_covs
  log_cluster_probs <- model_params$log_cluster_probs
  Z <- model_params$Z
  K <- dim(cluster_covs)[3]
  n <- nrow(Z)
  prob_mat <- matrix(log_cluster_probs, nrow = n, ncol = K, byrow = T)
  for (k in 1:K){
    prob_mat[, k] <- prob_mat[, k] + mvnfast::dmvn(Z, cluster_means[k, ],
                                                   cluster_covs[, , k], log = T)
  }
  model_params$clusters <- apply(exp(prob_mat), 1, function(r){sample(1:K, 1, prob = r)})
  model_params
}

update_cluster_means <- function(model_params, hyperpars){
  mu0 <- hyperpars$mu0
  h <- hyperpars$h
  cluster_covs <- model_params$cluster_covs
  Z <- model_params$Z
  clusters <- model_params$clusters
  p <- length(hyperpars$mu0)
  K <- dim(cluster_covs)[3]
  cluster_means <- matrix(nrow = K, ncol = p)
  for (k in 1:K){
    cluster_indices <- which(clusters == k)
    n_k <- length(cluster_indices)
    if (n_k == 0){
      cluster_means[k, ] <- mvnfast::rmvn(1, mu0, cluster_covs[, , k] / h)
    } else{
      Z_cluster <- Z[cluster_indices, , drop = F]
      mu <- (h * mu0 + apply(Z_cluster, 2, sum)) / (n_k + h)
      sigma <- cluster_covs[, , k] / (n_k + h)
      cluster_means[k, ] <- mvnfast::rmvn(1, mu, sigma)
    }
  }
  model_params$cluster_means <- cluster_means
  model_params
}

update_cluster_covs <- function(model_params, hyperpars){
  nu <- hyperpars$nu
  mu0 <- hyperpars$mu0
  h <- hyperpars$h
  Sigma <- model_params$Sigma
  cluster_means <- model_params$cluster_means
  Z <- model_params$Z
  clusters <- model_params$clusters
  K <- nrow(cluster_means)
  p <- ncol(Z)
  cluster_covs <- model_params$cluster_covs
  for (k in 1:K){
    cluster_indices <- which(clusters == k)
    n_k <- length(cluster_indices)
    if (n_k == 0){
      riwish(nu, Sigma)
    } else{
      Z_cluster <- Z[cluster_indices, , drop = F]
      Z_mean <- apply(Z_cluster, 2, mean)
      S <- matrix(0, nrow = p, ncol = p)
      for (i in 1:n_k){
        Zmz <- matrix(Z_cluster[i, ] - Z_mean, ncol = 1)
        S <- S + Zmz %*% t(Zmz)
      }
      mkm0 <- matrix(Z_mean - mu0, ncol = 1)
      M <- mkm0 %*% t(mkm0)
      cluster_covs[, , k] <- riwish(nu + n_k, (nu * n_k) / (nu + n_k) * M + S + Sigma)
    }
  }
  model_params$cluster_covs <- cluster_covs
  model_params
}


update_zobs <- function(model_params, transformations){
  Y <- model_params$Y
  E <- is.na(Y)
  Z <- model_params$Z
  clusters <- model_params$clusters
  cluster_means <- model_params$cluster_means
  cluster_covs <- model_params$cluster_covs
  n <- nrow(Z)
  p <- ncol(Z)
  K <- nrow(cluster_means)
  for (k in 1:K){
    meank <- cluster_means[k, , drop = F]
    cov <- cluster_covs[, , k]
    for (j in 1:p){
      obs_indices <- which(clusters == k & E[, j] == 0)
      if (length(obs_indices) == 0){
        next
      }
      Y_obs <- Y[obs_indices, , drop = F]
      Z_obs <- Z[obs_indices, , drop = F]
      boundsZ <- transformations$inverse_funs[[j]](Y_obs[, j])
      obs_discrete <- which(boundsZ[,1] != boundsZ[,2])
      if (length(obs_discrete) == 0){
        Z[obs_indices, j] <- boundsZ[, 1]
        next
      }
      Z_obs_disc <- Z_obs[obs_discrete, , drop = F]
      Sigma12 <- cov[j, -j, drop = F]
      Sigma22_inv <- solve(cov[-j, -j, drop = F])
      cond_var <- (cov[j, j] - Sigma12 %*% Sigma22_inv %*% t(Sigma12))[1, 1]
      cond_mean <- meank[j] + t(Sigma12 %*% Sigma22_inv %*%
                                                             t(Z_obs_disc[, -j, drop = F] -
                                                                 meank[-j]))
      Z[obs_indices[obs_discrete], j] <- truncnorm::rtruncnorm(length(obs_discrete),
                                                               a = boundsZ[obs_discrete, 1],
                                                               b = boundsZ[obs_discrete, 2],
                                                               mean = cond_mean,
                                                               sd = sqrt(cond_var))
    }
  }
  model_params$Z <- Z
  model_params
}

update_model_params <- function(model_params, hyperpars, transformations){
  model_params <- model_params %>%
    update_zobs(transformations) %>%
    update_cluster_covs(hyperpars) %>%
    update_cluster_means(hyperpars) %>%
    update_clusters() %>%
    update_sigma(hyperpars) %>%
    update_v() %>%
    calculate_log_cluster_probs() %>%
    update_alpha(hyperpars)
  model_params
}

midamix_mcmc <- function(inits, hyperpars, transformations, n_iter = 1000,
                         burnin = 100, monitor = NULL){
  model_params <- inits
  output <- list()
  for (variable in monitor){
    if (variable %in% names(model_params)){
      output[[variable]] = array(dim = c(0, dim(as.array(model_params[[variable]]))))
    }
  }
  if (burnin > 0){
    pb_burnin <- progress::progress_bar$new(
      format = "  burn-in [:bar] :percent eta: :eta",
      total = burnin, clear = FALSE, width= 60)
    pb_burnin$tick(0)
    for (i in 1:burnin){
      pb_burnin$tick()
      Sys.sleep(0.01)
      model_params <- update_model_params(model_params, hyperpars, transformations)
    }
  }
  pb_sampling <- progress::progress_bar$new(
    format = "  sampling [:bar] :percent eta: :eta",
    total = n_iter, clear = FALSE, width= 60)
  pb_sampling$tick(0)
  for (i in 1:n_iter){
    stop <- T
    tryCatch({
      model_params <- update_model_params(model_params, hyperpars, transformations)
      stop <- F
    }, error = function(e){

    })
    if (stop){
      message("\n Sampler failed: check output to diagnose the problem.")
      return(output)
    }
    for (variable in monitor){
      if (variable %in% names(model_params)){
        as.array(model_params[[variable]])
        output[[variable]] = abind::abind(output[[variable]], as.array(model_params[[variable]]),
                                          along = 1)
      }
    }
    pb_sampling$tick()
    Sys.sleep(0.01)
  }
  output
}


midamix_mcmc <- function(inits, hyperpars, transformations, validator = NULL,
                         cap = NULL, n_iter = 1000,
                         burnin = 100, monitor = NULL) {
  model_params <- inits
  output <- list()
  for (variable in monitor) {
    if (variable %in% names(model_params)) {
      output[[variable]] = rep(list(0), n_iter)
    }
  }
  if (burnin > 0) {
    pb_burnin <- progress::progress_bar$new(format = "  burn-in [:bar] :percent eta: :eta",
                                            total = burnin, clear = FALSE, width = 60)
    pb_burnin$tick(0)
    for (i in 1:burnin) {
      pb_burnin$tick()
      Sys.sleep(0.01)
      model_params <- update_model_params(model_params, hyperpars, transformations,
                                          validator, cap)
    }
  }
  pb_sampling <- progress::progress_bar$new(format = "  sampling [:bar] :percent eta: :eta",
                                            total = n_iter, clear = FALSE, width = 60)
  pb_sampling$tick(0)
  for (i in 1:n_iter) {
    stop <- T
    tryCatch({
      model_params <- update_model_params(model_params, hyperpars, transformations,
                                          validator, cap)
      stop <- F
    }, error = function(e) {

     })
    if (stop) {
      message("\n Sampler failed: check output to diagnose the problem.")
      return(model_params)
    }
    for (variable in monitor) {
      if (variable %in% names(model_params)) {
        output[[variable]][[i]] <- as.array(model_params[[variable]])
      }
    }
    pb_sampling$tick()
    Sys.sleep(0.01)
  }
  output <- lapply(output, simplify2array)
  output
}

initialize_pars <- function(data, max_clusters){
  p <- ncol(data)
  n <- nrow(data)
  K <- max_clusters
  Y <- data.matrix(data)
  funs <- lapply(1:p, function(index) {
    y <- na.omit(Y[, index])
    f <- function(z) {
      unname(quantile(y, probs = pnorm(z), type = 1))
    }
  })

  inverse_funs <- lapply(1:p, function(index) {
    y <- na.omit(Y[, index])
    y_aug <- c(y, -Inf, Inf)
    y_sort <- sort(unique(y_aug))
    cdf <- ecdf(y)
    f <- function(x) {
      x_inds <- match(x, y_sort)
      if (anyNA(x_inds)) {
        stop("Invalid inverse transformation.")
      }
      qnorm(cbind(cdf(y_sort[x_inds - 1]), cdf(x)))
    }
  })

  transformations <- list(funs = funs, inverse_funs = inverse_funs)
  hyperpars <- list(a_alpha = 0.5, b_alpha = 0.5, Sigma = diag(p), nu = p + 5, h = 4/3, mu0 = rep(0,
                                                                                                  p))
  alpha <- 1
  V <- c(rep(0.5, K - 1), 1)
  log_cluster_probs <- calculate_log_cluster_probs(list(V = V))$log_cluster_probs
  Z <- matrix(nrow = n, ncol = p)
  for (j in 1:p) {
    yj_obs <- which(!(is.na(Y[, j])))
    yj_mis <- which(is.na(Y[, j]))
    boundsZ_obs <- transformations$inverse_funs[[j]](Y[yj_obs, j])
    Z[yj_obs, j] <- truncnorm::rtruncnorm(length(yj_obs), a = boundsZ_obs[, 1], b = boundsZ_obs[,
                                                                                                2])
    Z[yj_mis, j] <- rnorm(length(yj_mis), 0, 1)
  }

  cluster_means <- matrix(0, nrow = K, ncol = p)
  clusters <- sample(1:3, n, replace = T)
  cluster_covs <- array(rep(diag(p), K), dim = c(p, p, K))
  cluster_precs <- array(rep(diag(p), K), dim = c(p, p, K))
  Zc <- matrix(nrow = 0, ncol = p)
  obs_map <- c()
  model_params <- list(Y = Y, Z = Z, clusters = clusters, cluster_means = cluster_means, cluster_covs = cluster_covs,
                       cluster_precs = cluster_precs, log_cluster_probs = log_cluster_probs, V = V,
                       alpha = alpha, Zc = Zc, obs_map = obs_map)
  list(model_params = model_params, transformations = transformations,
       hyperpars = hyperpars)
}

# update_alpha <- function(model_params, hyperpars){ a_alpha <- hyperpars$a_alpha b_alpha <-
# hyperpars$b_alpha cluster_probs <- model_params$cluster_probs K <- length(cluster_probs)
# alpha <- rgamma(1, a_alpha + K - 1, b_alpha - log(cluster_probs[K])) model_params$alpha <-
# alpha model_params } update_v <- function(model_params, hyperpars){ clusters <-
# model_params$clusters alpha <- model_params$alpha K <- length(model_params$cluster_probs)
# V <- rep(1, K) for (index in 1:(K-1)){ in_cluster <- sum(clusters == index) later_cluster
# <- sum(clusters > index) V[index] <- rbeta(1, 1 + in_cluster, alpha + later_cluster) }
# model_params$V <- V model_params } calculate_cluster_probs <- function(model_params){ V <-
# model_params$V K <- length(V) cluster_probs <- rep(0, K) log_cum_prod <- 0 for (k in 1:K){
# cluster_probs[k] <- log(V[k]) + log_cum_prod log_cum_prod <- log_cum_prod + log(1 - V[k])
# } model_params$cluster_probs <- exp(cluster_probs) model_params } update_sigma <-
# function(model_params, hyperpars){ nu <- hyperpars$nu omega <- hyperpars$omega Sigma0 <-
# hyperpars$Sigma0 cluster_precs <- model_params$cluster_precs p <- dim(cluster_precs)[2] K
# <- dim(cluster_precs)[3] S <- solve(solve(Sigma0) + apply(cluster_precs, c(1, 2), sum)) df
# <- nu * K + omega Sigma <- rwish(df, S) Sigma model_params$Sigma <- Sigma model_params }
# update_tau <- function(model_params, hyperpars){ a_tau <- hyperpars$a_tau b_tau <-
# hyperpars$b_tau B0 <- model_params$B0 cluster_coefs <- model_params$cluster_coefs q_star
# <- nrow(B0) p <- ncol(B0) K <- dim(cluster_coefs)[3] tau <- rep(1, p) BmB0 <-
# sweep(cluster_coefs, c(1, 2), B0, '-') ss <- rowSums(apply(BmB0, c(3), function(X){
# colSums(X ^ 2) })) tau <- rgamma(p, a_tau + q_star * K / 2, b_tau + ss / 2)
# model_params$tau <- tau model_params } update_b0 <- function(model_params, hyperpars){
# cluster_coefs <- model_params$cluster_coefs q <- dim(cluster_coefs)[1] p <-
# dim(cluster_coefs)[2] K <- dim(cluster_coefs)[3] sigma2b <- hyperpars$sigma2b tau <-
# model_params$tau sum_cluster <- apply(cluster_coefs, c(1, 2), sum) Omega <- diag(1 / tau)
# covar <- solve(Omega * K + diag(p) / sigma2b) B0 <- sum_cluster %*% Omega %*% covar +
# rmvnorm(q, rep(0, p), covar) model_params$B0 <- B0 model_params } update_clusters <-
# function(model_params, hyperpars){ cluster_probs <- model_params$cluster_probs
# cluster_coefs <- model_params$cluster_coefs cluster_precs <- model_params$cluster_precs Z
# <- model_params$Z X <- model_params$X n <- nrow(Z) p <- dim(cluster_coefs)[2] K <-
# dim(cluster_coefs)[3] prec_log_dets <- apply(cluster_precs, c(3), function(S)
# log(abs(det(S)))) clusters <- rep(0, n) for (i in 1:n){ Xi <- matrix(X[i, ], nrow = 1) Zi
# <- Z[i, ] log_density <- rep(0, K) for (k in 1:K){ log_density[k] <- dmvnorm_prec(Zi, Xi
# %*% cluster_coefs[, , k], cluster_precs[, , k], prec_log_det = prec_log_dets[k], log = T)
# } log_probs <- log(cluster_probs) + log_density log_probs <- log_probs - max(log_probs)
# probs <- exp(log_probs) / sum(exp(log_probs)) clusters[i] <- sample(1:K, 1, replace = F,
# prob = probs) } model_params$clusters <- clusters model_params } update_cluster_coefs <-
# function(model_params, hyperpars){ B0 <- model_params$B0 clusters <- model_params$clusters
# Z <- model_params$Z cluster_coefs <- model_params$cluster_coefs cluster_covs <-
# model_params$cluster_covs cluster_precs <- model_params$cluster_precs tau <-
# model_params$tau Omega <- diag(1 / tau) X <- model_params$X K <- dim(cluster_precs)[3] p
# <- dim(cluster_precs)[2] q <- ncol(X) for (k in 1:K){ cluster_indices <- which(clusters ==
# k) if (length(cluster_indices) == 0){ cluster_coefs[, , k] <- B0 + rmvnorm(q, rep(0, p),
# Omega) } else{ Z_cluster <- Z[cluster_indices, , drop = F] X_cluster <- X[cluster_indices,
# , drop = F] cov <- cluster_covs[, , k] prec <- cluster_precs[, , k] for (j in 1:p){
# Sigma12 <- cov[j, -j, drop = F] Sigma22_inv <- solve(cov[-j, -j, drop = F]) sigma_j_tilde
# <- (cov[j, j] - Sigma12 %*% Sigma22_inv %*% t(Sigma12))[1, 1] z_tilde <- Z_cluster[, j,
# drop = F] - t(Sigma12 %*% Sigma22_inv %*% t(Z_cluster[, -j, drop = F] - X_cluster %*%
# cluster_coefs[, -j, k, drop = F])) var_kj <- solve(tau[j] * diag(q) + t(X_cluster) %*%
# X_cluster / sigma_j_tilde) mean_kj <- var_kj %*% (tau[j] * B0[, j, drop = F] +
# t(X_cluster) %*% z_tilde / sigma_j_tilde) cluster_coefs[, j , k] <- rmvnorm(q, mean =
# mean_kj, sigma = var_kj) } } } cluster_coefs model_params$cluster_coefs <- cluster_coefs
# model_params } update_cluster_covs <- function(model_params, hyperpars){ cluster_coefs <-
# model_params$cluster_coefs cluster_covs <- model_params$cluster_covs cluster_precs <-
# model_params$cluster_precs p <- dim(cluster_covs)[2] Sigma <- model_params$Sigma X <-
# model_params$X Z <- model_params$Z nu <- hyperpars$nu for (k in 1:K){ cluster_indices <-
# which(clusters == k) if (length(cluster_indices) == 0){ cluster_covs[, , k] <- riwish(nu,
# Sigma) cluster_precs[, , k] <- solve(cluster_covs[, , k]) } else{ n_k <-
# length(cluster_indices) Z_cluster <- Z[cluster_indices, , drop = F] X_cluster <-
# X[cluster_indices, , drop = F] M <- X_cluster %*% cluster_coefs[, , k] Zmm <- Z_cluster -
# M S <- matrix(0, nrow = p, ncol = p) for (i in 1:n_k){ S = S + t(Zmm[i, , drop = F]) %*%
# Zmm[i, , drop = F] } cluster_covs[, , k] <- riwish(nu + n_k, Sigma + S) cluster_precs[, ,
# k] <- solve(cluster_covs[, , k]) } } model_params$cluster_covs <- cluster_covs
# model_params$cluster_precs <- cluster_precs model_params } update_zobs <-
# function(model_params, hyperpars, transformations){ Y <- model_params$Y E <- is.na(Y) Z <-
# model_params$Z X <- model_params$X clusters <- model_params$clusters cluster_coefs <-
# model_params$cluster_coefs cluster_covs <- model_params$cluster_covs n <- nrow(Z) p <-
# ncol(Z) K <- dim(cluster_coefs)[3] for (k in 1:K){ coefs <- matrix(cluster_coefs[, , k],
# ncol = p) cov <- cluster_covs[, , k] prec <- cluster_precs[, , k] for (j in 1:p){
# obs_indices <- which(clusters == k & E[, j] == 0) if (length(obs_indices) == 0){ next }
# Y_obs <- Y[obs_indices, , drop = F] Z_obs <- Z[obs_indices, , drop = F] X_obs <-
# X[obs_indices, , drop = F] boundsZ <- transformations$inverse_funs[[j]](Y_obs[, j])
# obs_discrete <- which(boundsZ[,1] != boundsZ[,2]) if (length(obs_discrete) == 0){
# Z[obs_indices, j] <- boundsZ[, 1] next } X_obs_disc <- X_obs[obs_discrete, , drop = F]
# Z_obs_disc <- Z_obs[obs_discrete, , drop = F] Sigma12 <- cov[j, -j, drop = F] Sigma22_inv
# <- solve(cov[-j, -j, drop = F]) cond_var <- (cov[j, j] - Sigma12 %*% Sigma22_inv %*%
# t(Sigma12))[1, 1] cond_mean <- X_obs_disc %*% coefs[, j, drop = F] - t(Sigma12 %*%
# Sigma22_inv %*% t(Z_obs_disc[, -j, drop = F] - X_obs_disc %*% coefs[, -j, drop = F]))
# Z[obs_indices[obs_discrete], j] <- truncnorm::rtruncnorm(length(obs_discrete), a =
# boundsZ[obs_discrete, 1], b = boundsZ[obs_discrete, 2], mean = cond_mean, sd =
# sqrt(cond_var)) } } model_params$Z <- Z model_params } update_model_params <-
# function(model_params, hyperpars, transformations){ model_params <- model_params %>%
# update_zobs(hyperpars, transformations) %>% update_cluster_covs(hyperpars) %>%
# update_cluster_coefs(hyperpars) %>% update_clusters(hyperpars) %>% update_b0(hyperpars)
# %>% update_sigma(hyperpars) %>% update_v(hyperpars) %>% calculate_cluster_probs() %>%
# update_tau(hyperpars) %>% update_alpha(hyperpars) model_params } midamix_mcmc <-
# function(inits, hyperpars, transformations, n_iter = 1000, burnin = 100, monitor = NULL){
# model_params <- inits output <- list() for (variable in monitor){ if (variable %in%
# names(model_params)){ output[[variable]] = array(dim = c(0,
# dim(as.array(model_params[[variable]])))) } } if (burnin > 0){ pb_burnin <-
# progress::progress_bar$new( format = ' burn-in [:bar] :percent eta: :eta', total = burnin,
# clear = FALSE, width= 60) pb_burnin$tick(0) for (i in 1:burnin){ pb_burnin$tick()
# Sys.sleep(0.01) model_params <- update_model_params(model_params, hyperpars,
# transformations) } } pb_sampling <- progress::progress_bar$new( format = ' sampling [:bar]
# :percent eta: :eta', total = n_iter, clear = FALSE, width= 60) pb_sampling$tick(0) for (i
# in 1:n_iter){ stop <- T tryCatch({ model_params <- update_model_params(model_params,
# hyperpars, transformations) stop <- F }, error = function(e){ }) if (stop){ message('\n
# Sampler failed: check output to diagnose the problem.') return(output) } for (variable in
# monitor){ if (variable %in% names(model_params)){ as.array(model_params[[variable]])
# output[[variable]] = abind::abind(output[[variable]], as.array(model_params[[variable]]),
# along = 1) } } pb_sampling$tick() Sys.sleep(0.01) } output }

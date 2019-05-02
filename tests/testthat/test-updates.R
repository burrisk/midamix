# p <- 3
# K <- 10
# q <- 1
# n <- 10000
#
# hyperpars <- list(
#   a_alpha = 0.5,
#   b_alpha = 0.5,
#   a_tau = 10,
#   b_tau = 10,
#   Sigma0 = diag(p),
#   nu = p + 1,
#   omega = p + 2,
#   sigma2b = 0.1
# )
#
# set.seed(314)
# tau <- rgamma(p, hyperpars$a_tau, hyperpars$b_tau)
# Omega <- diag(1 / tau)
# alpha <- rgamma(1, hyperpars$a_alpha, hyperpars$b_alpha)
# V <- c(rbeta(K - 1, 1, alpha), 1)
# cluster_probs <- calculate_cluster_probs(list(V = V))$cluster_probs
# B0 <- rmvnorm(q, rep(0, p), hyperpars$sigma2b * diag(p))
# Sigma <- rwish(hyperpars$omega, hyperpars$Sigma0)
# cluster_coefs <- array(0, dim = c(q, p, K))
# cluster_covs <- array(0, dim = c(p, p, K))
# cluster_precs <- array(0, dim = c(p, p, K))
#
# for (i in 1:K){
#   cluster_coefs[, , i] <- B0 + rmvnorm(q, rep(0, p), Omega)
#   cluster_precs[, , i] <- rwish(hyperpars$nu, Sigma)
#   cluster_covs[, , i] <- solve(cluster_precs[, , i])
# }
#
# clusters <- sample(1:K, n, replace = T, prob = cluster_probs)
# X <- matrix(1, nrow = n, ncol = 1)
# Z <- matrix(nrow = n, ncol = p)
# for (i in 1:n){
#   mu_i <- X[i, , drop = F] %*% cluster_coefs[, , clusters[i]]
#   Sigma_i <- cluster_covs[, , clusters[i]]
#   Z[i, ] <- rmvnorm(1, mu_i, Sigma_i)
# }
#
# transformations <- list(
#   funs = list(
#   f1 = function(z){
#     pmin(100, qpois(pnorm(z), lambda = 1))
#   },
#   f2 = function(z){
#     pmin(100, qpois(pnorm(z), lambda = 3))
#   },
#   f3 = function(z){
#     pmin(100, qpois(pnorm(z), lambda = 5))
#   }
# ),
#   inverse_funs = list(
#     f1_inv = function(y){
#       lower <- y - 1
#       upper <- y
#       cbind(qnorm(ppois(lower, lambda = 1)), qnorm(ppois(upper, lambda = 1)))
#     },
#     f2_inv <- function(y){
#       lower <- y - 1
#       upper <- y
#       cbind(qnorm(ppois(lower, lambda = 3)), qnorm(ppois(upper, lambda = 3)))
#     },
#     f3_inv <- function(y){
#       lower <- y - 1
#       upper <- y
#       cbind(qnorm(ppois(lower, lambda = 5)), qnorm(ppois(upper, lambda = 5)))
#     }
#   )
# )
#
# Y <- applyTransformations(Z, transformations$funs, indices = NULL)
#
# model_params <- list(
#   Y = Y,
#   Z = Z,
#   X = X,
#   clusters = clusters,
#   cluster_coefs = cluster_coefs,
#   cluster_covs = cluster_covs,
#   cluster_precs = cluster_precs,
#   B0 = B0,
#   Sigma = Sigma,
#   cluster_probs = cluster_probs,
#   V = V,
#   alpha = alpha,
#   tau = tau
# )
#
# update_model_params(model_params, hyperpars, transformations)

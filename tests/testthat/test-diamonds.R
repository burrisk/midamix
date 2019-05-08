# data(alabama)
#
# alabama <- purrr::map_df(alabama, function(x) {x[sample(c(TRUE, NA),
#                                                           prob = c(0.8, 0.2),
#                                                           size = length(x),
#                                                           replace = TRUE)]}) %>%
#   dplyr::sample_n(1000)
#
# validator <- function(y){
#   impossible <- (y[1] < 14 & y[7] > 0) |
#     (y[1] < 18 & y[3] > 40) |
#     (y[1] < 65 & y[2] == 1 & y[4] == 0) |
#     (y[6] == 0 & y[7] > 0)
#   !(impossible)
# }
#
# imps <- impute(alabama, n_iter = 100, burnin = 20, validator = validator)
#
# lapply(imps$data, function(df){
#   mean(apply(df, 1, validator))
# })
# df1 <- imps$data[[8]]
#
#
# mean(apply(df1, 1, validator))
# n <- nrow(Y)
# p <- ncol(Y)
# K <- 15
#
# funs <- lapply(1:p, function(index) {
#   y <- na.omit(Y[, index])
#   f <- function(z) {
#     unname(quantile(y, probs = pnorm(z), type = 1))
#   }
# })
#
# inverse_funs <- lapply(1:p, function(index) {
#   y <- na.omit(Y[, index])
#   y_aug <- c(y, -Inf, Inf)
#   y_sort <- sort(unique(y_aug))
#   cdf <- ecdf(y)
#   f <- function(x) {
#     x_inds <- match(x, y_sort)
#     if (anyNA(x_inds)) {
#       stop("Invalid inverse transformation.")
#     }
#     qnorm(cbind(cdf(y_sort[x_inds - 1]), cdf(x)))
#   }
# })
#
# transformations <- list(funs = funs, inverse_funs = inverse_funs)
#
# hyperpars <- list(
#   a_alpha = 0.5,
#   b_alpha = 0.5,
#   Sigma = diag(p),
#   nu = p + 5,
#   h = 4/3 ,
#   mu0 = rep(0, p)
# )
#
# # Initial values
# alpha <- 1
# V <- c(rep(0.5, K - 1), 1)
# log_cluster_probs <- calculate_log_cluster_probs(list(V = V))$log_cluster_probs
# Z <- matrix(nrow = n, ncol = p)
#
# hyperpars <- list(a_alpha = 0.5, b_alpha = 0.5, Sigma = diag(p), nu = p + 5, h = 4/3, mu0 = rep(0,
#                                                                                                 p))
# alpha <- 1
# V <- c(rep(0.5, K - 1), 1)
# log_cluster_probs <- calculate_log_cluster_probs(list(V = V))$log_cluster_probs
# Z <- matrix(nrow = n, ncol = p)
# for (j in 1:p) {
#   yj_obs <- which(!(is.na(Y[, j])))
#   yj_mis <- which(is.na(Y[, j]))
#   boundsZ_obs <- transformations$inverse_funs[[j]](Y[yj_obs, j])
#   Z[yj_obs, j] <- truncnorm::rtruncnorm(length(yj_obs), a = boundsZ_obs[, 1], b = boundsZ_obs[,
#                                                                                               2])
#   Z[yj_mis, j] <- rnorm(length(yj_mis), 0, 1)
# }
#
# cluster_means <- rmvn(K, mu = hyperpars$mu0, sigma = hyperpars$Sigma)
# clusters <- sample(1:3, n, replace = T)
# cluster_covs <- array(rep(diag(p), K), dim = c(p, p, K))
# cluster_precs <- array(rep(diag(p), K), dim = c(p, p, K))
# for (i in 1:K) {
#   cluster_covs[, , i] <- riwish(hyperpars$nu, hyperpars$Sigma)
#   cluster_precs[, , i] <- solve(cluster_covs[, , i])
# }
# model_params <- list(Y = Y, Z = Z, clusters = clusters, cluster_means = cluster_means, cluster_covs = cluster_covs,
#                      cluster_precs = cluster_precs, log_cluster_probs = log_cluster_probs, V = V, alpha = alpha)
#
# n_iter <- 100
#
# samps <- midamix_mcmc(model_params, hyperpars, transformations, n_iter = n_iter,
#                       validator = validator,
#                       burnin = 20, monitor = c("Z", "cluster_means", "clusters",
#                                                "log_cluster_probs", "cluster_covs",
#                                                "alpha", "Sigma", "V", "Zc"))
# #
# # plot(samps$alpha, type = 'l')
# #
# # Z_df <- as.data.frame(model_params$Z)
# # Z_df$cluster <- factor(model_params$clusters, levels = 1:K)
# # Z_df$cluster <- factor(samps$clusters[n_iter, ], levels = 1:K)
# #
# # ggplot2::ggplot(Z_df, ggplot2::aes(x = V1, y = V2, colour = cluster)) +
# #   ggplot2::geom_point()
# #
# # sigma_list <- rep(list(1), K)
# # for (i in 1:K){
# #   sigma_list[[i]] <- samps$cluster_covs[n_iter, , , i]
# # }
# #
# # Z_mix <- mvnfast::rmixn(n, mu = samps$cluster_means[n_iter, , ],
# #                                  sigma =  sigma_list,
# #                                  w = exp(samps$log_cluster_probs[n_iter, ]))
# #
# # plot(Z_mix[,1], Z_mix[,2])
# # plot(model_params$Z[,1], model_params$Z[,2])
# # Y_mix <- applyTransformations(Z_mix, transformations$funs)
# # plot(Y_mix[,1], Y_mix[,2])
# # plot(Y[,1], Y[,2])
# #
# # diamonds_mice <- mice(diamonds, m = 5)
# # fit <- with(data = diamonds_mice, exp = lm(log(price) ~ color + clarity + carat))
# # summary(pool(fit))

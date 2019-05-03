# diamonds <- ggplot2::diamonds
#
# set.seed(314)
# diamonds_missing <- purrr::map_df(diamonds, function(x) {x[sample(c(TRUE, NA),
#                                                                   prob = c(0.8, 0.2),
#                                                                   size = length(x),
#                                                                   replace = TRUE)]}) %>%
#   dplyr::sample_n(1000)
#
# imputations <- impute(diamonds_missing, burnin = 0, n_iter = 100, seed = 314)
#
# my_model <- function(data){
#   lm(log(price) ~ carat + depth + table, data = data)
# }
#
# model_fits <- imputations %>%
#   fit_model(my_model)
#
# model_results <- model_fits %>%
#   pool_inferences()
# #
# mice_imps <- mice(diamonds_missing, m = 10)
#
# comp_1 <- mice::complete(mice_imps, action = 1)
#
# plot(comp_1[,1], comp_1[,7], col  = "red")
# points(imps[[1, 2]]$carat, imps[[1, 2]]$price, col = "blue")
# points(diamonds_missing$carat, diamonds_missing$price, col = "black")
#
#
# fit <- with(data = mice_imps, exp = lm(price ~ carat))
# summary(pool(fit))
#
# str(imps[[1, 2]])
#
# df1 <- imps[[1, 2]]
# df2 <- imps[[3, 2]]
#
# plot(df1$carat, df1$x, col = "red")
# points(df2$carat, df2$x, col = "blue")
# plot(diamonds_missing$carat, diamonds_missing$x)
#
# last_imp <- imps[[10]]
# plot(diamonds$carat, diamonds$price)
# plot(diamonds_missing$carat, diamonds_missing$price)
# points(last_imp[,1], last_imp[, 7], col = "red")
#
# # set.seed(313)
# # cluster_means <- rbind(c(0, 0), c(0, 0), c(0, 0))
# # cluster_covs <- rep(list(0.8 + 0.2 * diag(2)), 3)
# # cluster_probs <- c(1/2, 0, 1/2)
# # n <- 1000
# # Z <- mvnfast::rmixn(n, mu = cluster_means,
# #                     sigma = cluster_covs,
# #                     w = cluster_probs)
# #
# # Y <- qpois(pnorm(Z), lambda = 10)
#
# n <- nrow(Y)
# p <- ncol(Y)
# K <- 15
#
# funs <- lapply(1:p, function(index){
#   y <- Y[, index]
#   f <- function(z){
#     unname(quantile(y, probs = pnorm(z), type = 1))
#   }
# })
#
# inverse_funs <- lapply(1:p, function(index){
#   y <- Y[, index]
#   y_aug <- c(y, -Inf, Inf)
#   y_sort <- sort(unique(y_aug))
#   cdf <- ecdf(y)
#   f <- function(x){
#     x_inds <- match(x, y_sort)
#     if (anyNA(x_inds)){
#       stop("Invalid inverse transformation.")
#     }
#     qnorm(cbind(cdf(y_sort[x_inds - 1]), cdf(x)))
#   }
# })
#
# transformations <- list(
#   funs = funs,
#   inverse_funs = inverse_funs
# )
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
# for (j in 1:p){
#   boundsZ <- transformations$inverse_funs[[j]](Y[, j])
#   Z[, j] <- truncnorm::rtruncnorm(n, a = boundsZ[,1], b = boundsZ[,2])
# }
#
# cluster_means <- rmvnorm(K, mean = hyperpars$mu0, sigma = hyperpars$Sigma)
# clusters <- sample(1:3, n, replace = T)
# cluster_covs <- array(rep(diag(p), K), dim = c(p, p, K))
# cluster_precs <- array(rep(diag(p), K), dim = c(p, p, K))
# for (i in 1:K){
#   cluster_covs[, , i] <- riwish(hyperpars$nu, hyperpars$Sigma)
#   cluster_precs[, , i] <- solve(cluster_covs[, , i])
# }
#
#
# model_params <- list(
#   Y = Y,
#   Z = Z,
#   clusters = clusters,
#   cluster_means = cluster_means,
#   cluster_covs = cluster_covs,
#   cluster_precs = cluster_precs,
#   log_cluster_probs = log_cluster_probs,
#   V = V,
#   alpha = alpha
# )
#
# n_iter <- 100
#
# samps <- midamix_mcmc(model_params, hyperpars, transformations, n_iter = n_iter,
#                       burnin = 0, monitor = c("Z", "cluster_means", "clusters",
#                                                "log_cluster_probs", "cluster_covs",
#                                                "alpha", "Sigma", "V"))
#
# plot(samps$alpha, type = 'l')
#
# Z_df <- as.data.frame(model_params$Z)
# Z_df$cluster <- factor(model_params$clusters, levels = 1:K)
# Z_df$cluster <- factor(samps$clusters[n_iter, ], levels = 1:K)
#
# ggplot2::ggplot(Z_df, ggplot2::aes(x = V1, y = V2, colour = cluster)) +
#   ggplot2::geom_point()
#
# sigma_list <- rep(list(1), K)
# for (i in 1:K){
#   sigma_list[[i]] <- samps$cluster_covs[n_iter, , , i]
# }
#
# Z_mix <- mvnfast::rmixn(n, mu = samps$cluster_means[n_iter, , ],
#                                  sigma =  sigma_list,
#                                  w = exp(samps$log_cluster_probs[n_iter, ]))
#
# plot(Z_mix[,1], Z_mix[,2])
# plot(model_params$Z[,1], model_params$Z[,2])
# Y_mix <- applyTransformations(Z_mix, transformations$funs)
# plot(Y_mix[,1], Y_mix[,2])
# plot(Y[,1], Y[,2])
#
# diamonds_mice <- mice(diamonds, m = 5)
# fit <- with(data = diamonds_mice, exp = lm(log(price) ~ color + clarity + carat))
# summary(pool(fit))

#' Multiple imputation for data analysis via mixtures
#'
#' Uses a Dirichlet Process Mixture Transformation Model to sample missing data values
#' from the corresponding posterior distribution.
#'
#' @param data A tibble, data frame, or matrix, consisting of numeric, integer,
#' or ordered factor columns.
#' @param imputations The number of imputed datasets to return, defaults to 10.
#' @param max_clusters The maximum number of clusters for the mixture model.
#' @param n_iter Number of iterations for the MCMC sampler.
#' @param burnin Number of iterations for initial burn-in period.
#' @param validator A function that indicates the validity of imputations.  If \code{NULL},
#' the data are assumed unconstrained.
#' @param transformations Transformations of response variables to latent scale, defaults
#' to empirical CDFs
#' @param seed Random seed.
#' @return A \code{tibble} consisting of multiply imputed data sets.
#' @export
midamix <- function(data, imputations = 10, max_clusters = 15, n_iter = 1000, burnin = 100,
    validator = NULL, transformations = NULL, seed = NA) {
    data <- dplyr::as_tibble(data)
    doubles <- sapply(data, is.double)
    integers <- sapply(data, is.integer)
    ordered_factors <- sapply(data, is.ordered)
    p <- sum(doubles) + sum(integers) + sum(ordered_factors)
    n <- nrow(data)
    K <- max_clusters
    if (p != ncol(data)) {
        stop("All variables must either be numeric, integer, or ordered factors.")
    }
    if (is.na(seed)) {
        set.seed(314)
    } else {
        set.seed(seed)
    }
    if (is.null(validator)) {
        validator <- function(y) {
            TRUE
        }
    }
    Y <- data.matrix(data)
    if (is.null(transformations)) {
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
    }
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

    cluster_means <- rmvn(K, mu = hyperpars$mu0, sigma = hyperpars$Sigma)
    clusters <- sample(1:3, n, replace = T)
    cluster_covs <- array(rep(diag(p), K), dim = c(p, p, K))
    cluster_precs <- array(rep(diag(p), K), dim = c(p, p, K))
    for (i in 1:K) {
        cluster_covs[, , i] <- riwish(hyperpars$nu, hyperpars$Sigma)
        cluster_precs[, , i] <- solve(cluster_covs[, , i])
    }
    model_params <- list(Y = Y, Z = Z, clusters = clusters, cluster_means = cluster_means, cluster_covs = cluster_covs,
        cluster_precs = cluster_precs, log_cluster_probs = log_cluster_probs, V = V, alpha = alpha)
    samps <- midamix_mcmc(model_params, hyperpars, transformations, n_iter = n_iter, burnin = burnin,
        monitor = c("Z"))
    imputation_indices <- seq(0, n_iter, length = imputations + 1)[-1]
    imp_list <- lapply(imputation_indices, function(index) {
        applyTransformations(samps$Z[, , index], transformations$funs)
    })
    imp_tibbles <- lapply(imp_list, function(df) {
        df <- df %>% dplyr::as_tibble() %>% dplyr::mutate_at(which(integers), as.integer) %>%
            dplyr::mutate_at(which(ordered_factors), as.ordered)
        colnames(df) <- colnames(data)
        for (j in which(ordered_factors)) {
            levels(df[[j]]) <- levels(data[[j]])
        }
        df
    })
    result <- tibble::tibble(imputation = 1:imputations, data = imp_tibbles)
    return(result)
}

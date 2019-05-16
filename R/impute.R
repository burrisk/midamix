#' Multiple imputation for data analysis via mixtures
#'
#' Uses a Dirichlet Process Mixture Transformation Model to sample missing data values
#' from the corresponding posterior distribution.
#'
#' @param data A data frame, consisting of numeric, integer,
#' or ordered factor columns.
#' @param imputations The number of imputed datasets to return, defaults to 10.
#' @param max_clusters The maximum number of clusters for the mixture model.
#' @param n_iter Number of iterations for the MCMC sampler.
#' @param burnin Number of iterations for initial burn-in period.
#' @param validator A function that takes in an observation and determines whether
#' it is feasible.
#' @param seed Random seed.
#' @return A \code{tibble} consisting of multiply imputed data sets.
#' @export
impute <- function(data, imputations = 10, max_clusters = 15, n_iter = 1000, burnin = 100,
                   validator = NULL, seed = 1) {
    if (!(is.data.frame(data))) {
        stop("Argument data must be a data frame.")
    }
    if (!(is.numeric(imputations)  | imputations <= 0)) {
        stop("The number of imputed data sets must be a positive integer.")
    }
    if (!(is.numeric(max_clusters)) | max_clusters <= 0) {
        stop("The maximum number of clusters must be a positive integer.")
    }
    if (!(is.numeric(n_iter)) | n_iter < imputations) {
        stop("The number of iterations must be an integer at least equal to number
             of imputed data sets.")
    }
    if (!(is.numeric(burnin)) | burnin < 0) {
        stop("The number of burn-in iterations must be a non-nnegative integer.")
    }
    data <- dplyr::as_tibble(data)
    doubles <- sapply(data, is.double)
    integers <- sapply(data, is.integer)
    ordered_factors <- sapply(data, is.ordered)
    non_categorical <- which(doubles | integers | ordered_factors)
    data_numeric <- data[, non_categorical]
    any_observed <- sapply(data_numeric, function(x) !(all(is.na(x))))
    if (sum(any_observed) != ncol(data_numeric)) {
        warning("There is a columnn with all NA values.  This is excluded when
                performing multiple imputation.")
        data_numeric <- data_numeric[, any_observed]
    }
    data_categorical <- data[, -(non_categorical)]
    # p <- ncol(data_numeric)
    # n <- nrow(data_numeric)
    # K <- max_clusters
    # set.seed(seed)
    # Y <- data.matrix(data_numeric)
    # funs <- lapply(1:p, function(index) {
    #     y <- na.omit(Y[, index])
    #     f <- function(z) {
    #         unname(quantile(y, probs = pnorm(z), type = 1))
    #     }
    # })
    #
    # inverse_funs <- lapply(1:p, function(index) {
    #     y <- na.omit(Y[, index])
    #     y_aug <- c(y, -Inf, Inf)
    #     y_sort <- sort(unique(y_aug))
    #     cdf <- ecdf(y)
    #     f <- function(x) {
    #         x_inds <- match(x, y_sort)
    #         if (anyNA(x_inds)) {
    #           stop("Invalid inverse transformation.")
    #         }
    #         qnorm(cbind(cdf(y_sort[x_inds - 1]), cdf(x)))
    #     }
    # })
    #
    # transformations <- list(funs = funs, inverse_funs = inverse_funs)
    # hyperpars <- list(a_alpha = 0.5, b_alpha = 0.5, Sigma = diag(p), nu = p + 5, h = 4/3, mu0 = rep(0,
    #     p))
    # alpha <- 1
    # V <- c(rep(0.5, K - 1), 1)
    # log_cluster_probs <- calculate_log_cluster_probs(list(V = V))$log_cluster_probs
    # Z <- matrix(nrow = n, ncol = p)
    # for (j in 1:p) {
    #     yj_obs <- which(!(is.na(Y[, j])))
    #     yj_mis <- which(is.na(Y[, j]))
    #     boundsZ_obs <- transformations$inverse_funs[[j]](Y[yj_obs, j])
    #     Z[yj_obs, j] <- truncnorm::rtruncnorm(length(yj_obs), a = boundsZ_obs[, 1], b = boundsZ_obs[,
    #         2])
    #     Z[yj_mis, j] <- rnorm(length(yj_mis), 0, 1)
    # }
    #
    # cluster_means <- matrix(0, nrow = K, ncol = p)
    # clusters <- sample(1:3, n, replace = T)
    # cluster_covs <- array(rep(diag(p), K), dim = c(p, p, K))
    # cluster_precs <- array(rep(diag(p), K), dim = c(p, p, K))
    # Zc <- matrix(nrow = 0, ncol = p)
    # obs_map <- c()
    # model_params <- list(Y = Y, Z = Z, clusters = clusters, cluster_means = cluster_means, cluster_covs = cluster_covs,
    #     cluster_precs = cluster_precs, log_cluster_probs = log_cluster_probs, V = V,
    #     alpha = alpha, Zc = Zc, obs_map = obs_map)
    set.seed(seed)
    inits <- initialize_pars(data_numeric, max_clusters)
    samps <- midamix_mcmc(inits = inits$model_params,
                          hyperpars = inits$hyperpars,
                          transformations = inits$transformations,
                          validator = validator,
                          n_iter = n_iter, burnin = burnin,
        monitor = c("Z"))
    imputation_indices <- seq(0, n_iter, length = imputations + 1)[-1]
    imp_list <- lapply(imputation_indices, function(index) {
        applyTransformations(samps$Z[, , index], inits$transformations$funs)
    })
    imp_tibbles <- lapply(imp_list, function(df) {
        colnames(df) <- colnames(data_numeric)
        df <- df %>% dplyr::as_tibble() %>% dplyr::mutate_at(which(integers), as.integer) %>%
            dplyr::mutate_at(which(ordered_factors), as.ordered)
        for (j in which(ordered_factors)) {
            levels(df[[j]]) <- levels(data_numeric[[j]])
        }
        dplyr::bind_cols(data_categorical, df)
    })
    result <- tibble::tibble(imputation = 1:imputations, data = imp_tibbles)
    return(result)
}

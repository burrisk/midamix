update_alpha <- function(model_params, hyperpars) {
    log_cluster_probs <- model_params$log_cluster_probs
    K <- length(log_cluster_probs)
    model_params$alpha <- rgamma(1, hyperpars$a_alpha + K - 1, hyperpars$b_alpha - log_cluster_probs[K])
    model_params
}

update_v <- function(model_params) {
    alpha <- model_params$alpha
    obs_map <- model_params$obs_map
    clusters <- c(model_params$clusters, model_params$clusters[obs_map])
    K <- length(model_params$log_cluster_probs)
    V <- rep(1, K)
    for (k in 1:(K - 1)) {
        in_cluster <- sum(clusters == k)
        later_cluster <- sum(clusters > k)
        V[k] <- min(1 - 1e-06, rbeta(1, 1 + in_cluster, alpha + later_cluster))
    }
    model_params$V <- V
    model_params
}

calculate_log_cluster_probs <- function(model_params) {
    V <- model_params$V
    K <- length(V)
    log_cum_probs <- c(0, cumsum(log(1 - V[-K])))
    model_params$log_cluster_probs <- log(V) + log_cum_probs
    model_params
}

update_clusters <- function(model_params) {
    cluster_means <- model_params$cluster_means
    cluster_covs <- model_params$cluster_covs
    log_cluster_probs <- model_params$log_cluster_probs
    Z <- model_params$Z
    Zc <- model_params$Zc
    K <- dim(cluster_covs)[3]
    n <- nrow(Z)
    prob_mat <- matrix(log_cluster_probs, nrow = n, ncol = K, byrow = T)
    if (nrow(model_params$Zc) == 0){
        for (k in 1:K) {
            prob_mat[, k] <- prob_mat[, k] + dmvn(Z, cluster_means[k, ], cluster_covs[,
                                                                                      , k], log = T)
        }
        model_params$clusters <- apply(exp(prob_mat), 1, function(r) {
            sample(1:K, 1, prob = r)
        })
        return(model_params)
    } else{
        obs_map <- model_params$obs_map
        for (i in 1:n){
            Zc_obs <- Zc[which(obs_map == i), ]
            Z_obs <- rbind(Z[i, ], Zc_obs)
            for (k in 1:K) {
                prob_mat[i, k] <- prob_mat[i, k] + sum(dmvn(Z_obs, cluster_means[k, ], cluster_covs[,
                                                                       , k], log = T))
            }
        }
        model_params$clusters <- apply(exp(prob_mat), 1, function(r) {
            sample(1:K, 1, prob = r)
        })
        model_params
    }
}

update_cluster_means <- function(model_params, hyperpars) {
    mu0 <- matrix(hyperpars$mu0, ncol = 1)
    h <- hyperpars$h
    cluster_precs <- model_params$cluster_precs
    Z <- rbind(model_params$Z, model_params$Zc)
    obs_map <- model_params$obs_map
    clusters <- c(model_params$clusters, model_params$clusters[obs_map])
    p <- length(hyperpars$mu0)
    K <- dim(cluster_precs)[3]
    cluster_means <- matrix(nrow = K, ncol = p)
    for (k in 1:K) {
        cluster_indices <- which(clusters == k)
        n_k <- length(cluster_indices)
        if (n_k == 0) {
            cluster_means[k, ] <- rmvn(1, mu0, diag(p)/h)
        } else {
            prec <- cluster_precs[, , k]
            Z_cluster <- Z[cluster_indices, , drop = F]
            z_sums <- matrix(apply(Z_cluster, 2, sum), ncol = 1)
            sigma <- solve(n_k * prec + h * diag(p))
            mu <- sigma %*% (prec %*% z_sums + h * mu0)
            cluster_means[k, ] <- rmvn(1, mu, sigma)
        }
    }
    model_params$cluster_means <- cluster_means
    model_params
}

update_cluster_covs <- function(model_params, hyperpars) {
    nu <- hyperpars$nu
    mu0 <- hyperpars$mu0
    h <- hyperpars$h
    Sigma <- hyperpars$Sigma
    cluster_means <- model_params$cluster_means
    Z <- rbind(model_params$Z, model_params$Zc)
    obs_map <- model_params$obs_map
    clusters <-  c(model_params$clusters, model_params$clusters[obs_map])
    K <- nrow(cluster_means)
    p <- ncol(Z)
    cluster_covs <- model_params$cluster_covs
    cluster_precs <- model_params$cluster_precs
    for (k in 1:K) {
        cluster_indices <- which(clusters == k)
        n_k <- length(cluster_indices)
        if (n_k == 0) {
            cluster_covs[, , k] <- riwish(nu, Sigma)
            cluster_precs[, , k] <- solve(cluster_covs[, , k])
        } else {
            Z_cluster <- Z[cluster_indices, , drop = F]
            Z_mean <- apply(Z_cluster, 2, mean)
            S <- matrix(0, nrow = p, ncol = p)
            for (i in 1:n_k) {
                Zmz <- matrix(Z_cluster[i, ] - Z_mean, ncol = 1)
                S <- S + Zmz %*% t(Zmz)
            }
            Cov_k <- riwish(nu + n_k, S + Sigma)
            cluster_covs[, , k] <- Cov_k
            cluster_precs[, , k] <- solve(Cov_k)
        }
    }
    model_params$cluster_covs <- cluster_covs
    model_params$cluster_precs <- cluster_precs
    model_params
}


update_zobs <- function(model_params, transformations) {
    Y <- model_params$Y
    E <- is.na(Y)
    Z <- model_params$Z
    clusters <- model_params$clusters
    cluster_means <- model_params$cluster_means
    cluster_covs <- model_params$cluster_covs
    n <- nrow(Z)
    p <- ncol(Z)
    K <- nrow(cluster_means)
    for (k in 1:K) {
        meank <- cluster_means[k, , drop = F]
        cov <- cluster_covs[, , k]
        for (j in 1:p) {
            obs_indices <- which(clusters == k & E[, j] == 0)
            if (length(obs_indices) == 0) {
                next
            }
            Y_obs <- Y[obs_indices, , drop = F]
            Z_obs <- Z[obs_indices, , drop = F]
            boundsZ <- transformations$inverse_funs[[j]](Y_obs[, j])
            obs_discrete <- which(boundsZ[, 1] != boundsZ[, 2])
            if (length(obs_discrete) == 0) {
                Z[obs_indices, j] <- boundsZ[, 1]
                next
            }
            Z_obs_disc <- Z_obs[obs_discrete, , drop = F]
            Sigma12 <- cov[j, -j, drop = F]
            Sigma22_inv <- solve(cov[-j, -j, drop = F])
            cond_var <- (cov[j, j] - Sigma12 %*% Sigma22_inv %*% t(Sigma12))[1, 1]
            cond_mean <- meank[j] + t(Sigma12 %*% Sigma22_inv %*% t(sweep(Z_obs_disc[, -j, drop = F],
                2, meank[-j])))
            Z[obs_indices[obs_discrete], j] <- truncnorm::rtruncnorm(length(obs_discrete), a = boundsZ[obs_discrete,
                1], b = boundsZ[obs_discrete, 2], mean = cond_mean, sd = sqrt(cond_var))
        }
    }
    model_params$Z <- Z
    model_params
}

update_zmis <- function(model_params, transformations, validator) {
    Y <- model_params$Y
    Z <- model_params$Z
    n <- nrow(Y)
    p <- ncol(Y)
    E <- is.na(Y)
    missing_observations <- which(rowSums(E) != 0)
    for (i in missing_observations) {
        missing_values <- which(E[i, ] == 1)
        obs_values <- which(E[i, ] == 0)
        cluster <- model_params$clusters[i]
        cluster_mean <- t(model_params$cluster_means[cluster, , drop = F])
        cluster_cov <- model_params$cluster_covs[, , cluster]
        if (length(obs_values) == 0) {
            cond_mean <- cluster_mean
            cond_var <- cluster_cov
        } else {
            Sigma22_inv <- solve(cluster_cov[obs_values, obs_values])
            Sigma12 <- cluster_cov[missing_values, obs_values, drop = F]
            cond_mean <- cluster_mean[missing_values, drop = F] + Sigma12 %*% Sigma22_inv %*%
                (Z[i, obs_values] - cluster_mean[obs_values])
            cond_var <- cluster_cov[missing_values, missing_values] - Sigma12 %*% Sigma22_inv %*%
                t(Sigma12)
        }
        if (is.null(validator)){
            z_mis <- rmvn(1, cond_mean, cond_var)
        } else{
            in_region <- F
            z_star <- Z[i, , drop = F]
            while(!(in_region)){
                z_mis <- rmvn(1, cond_mean, cond_var)
                z_star[1, missing_values] <- z_mis
                y_star <- applyTransformations(z_star, transformations$funs)
                in_region <- validator(y_star)
            }
        }
        Z[i, missing_values] <- z_mis
    }
    model_params$Z <- Z
    model_params
}

update_zc <- function(model_params, transformations, validator){
    Z <- model_params$Z
    n <- nrow(Z)
    p <- ncol(Z)
    Zc <- matrix(nrow = 0, ncol = p)
    obs_map <- c()
    if (is.null(validator)){
        model_params$obs_map <- obs_map
        model_params$Zc <- Zc
        return(model_params)
    }
    clusters <- model_params$clusters
    cluster_means <- model_params$cluster_means
    cluster_covs <- model_params$cluster_covs
    for (i in 1:n){
        cl <- clusters[i]
        mu <- cluster_means[cl, ]
        sigma <- cluster_covs[, , cl]
        accepted <- F
        while(!(accepted)){
            Z_proposed <- rmvn(1, mu, sigma)
            Y_proposed <- applyTransformations(Z_proposed, transformations$funs)
            accepted <- apply(Y_proposed, 1, validator)
            if (!(accepted)){
                Zc <- rbind(Zc, Z_proposed)
                obs_map <- c(obs_map, i)
            }
        }
    }

    model_params$Zc <- Zc
    model_params$obs_map <- obs_map
    model_params
}

update_model_params <- function(model_params, hyperpars, transformations, validator) {
    model_params <- model_params %>%
        update_zmis(transformations, validator) %>%
        update_zobs(transformations) %>%
        update_cluster_covs(hyperpars) %>%
        update_cluster_means(hyperpars) %>%
        update_zc(transformations, validator) %>%
        update_clusters() %>%
        update_v() %>%
        calculate_log_cluster_probs() %>%
        update_alpha(hyperpars)
    model_params
}

midamix_mcmc <- function(inits, hyperpars, transformations, validator = NULL, n_iter = 1000,
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
                                                validator)
        }
    }
    pb_sampling <- progress::progress_bar$new(format = "  sampling [:bar] :percent eta: :eta",
        total = n_iter, clear = FALSE, width = 60)
    pb_sampling$tick(0)
    for (i in 1:n_iter) {
        stop <- T
        tryCatch({
            model_params <- update_model_params(model_params, hyperpars, transformations,
                                                validator)
            stop <- F
        }, error = function(e) {

        })
        if (stop) {
            message("\n Sampler failed: check output to diagnose the problem.")
            return(output)
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


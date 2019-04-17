# Generates data from a mixture of normal distributions
rMixtureNormal <- function(n,
                           mean = matrix(0, nrow = dim(sigma)[1], ncol = num_clusters),
                           sigma = array(rep(diag(nrow(mean)), num_clusters),
                                         dim = c(nrow(mean), nrow(mean), num_clusters)),
                           num_clusters = 1,
                           prob = rep(1, num_clusters) / num_clusters,
                           method = c("eigen", "svd", "chol")){

  # Conforming size of mean, sigma, num_clusters, and probs
  if (nrow(mean) != dim(sigma)[1] | nrow(mean) != dim(sigma)[2]){
    stop("mean and sigma have non-conforming size.")
  }

  if (ncol(mean) != num_clusters){
    stop("mean and num_clusters have non-conforming size.")
  }

  if (length(prob) != num_clusters){
    stop("prob and num_clusters have non-conforming size")
  }

  cluster_memberships <- sample(1:num_clusters, size = n, replace = TRUE, prob = prob)
  sim_list <- lapply(1:num_clusters, function(cluster){
    num_members <- sum(cluster_memberships == cluster)
    mvtnorm::rmvnorm(num_members, mean = mean[, cluster], sigma = sigma[, , cluster],
                     method = method)
  })
  sim <- do.call(rbind, sim_list)
  sim
}

rMixtureNormalTransform <- function(n,
                                    mean = matrix(0, nrow = dim(sigma)[1], ncol = num_clusters),
                                    sigma = array(rep(diag(nrow(mean)), num_clusters),
                                                  dim = c(nrow(mean), nrow(mean), num_clusters)),
                                    num_clusters = 1,
                                    prob = rep(1, num_clusters) / num_clusters,
                                    transformations = rep(list(function(y) y), nrow(mean)),
                                    method = c("eigen", "svd", "chol")){

  # Conforming size of transformations and mean
  if (length(transformations) != nrow(mean)){
    stop("The number of transformations must be equal to the number of variables.")
  }
  latent_variables <- rMixtureNormal(n, mean = mean, sigma = sigma, num_clusters = num_clusters,
                                     prob = prob, method = method)
  output <- applyTransformations(latent_variables, transformations)
  output

}

# Examples
n <- 100
num_clusters <- 3
p <- 2
mean <- matrix(c(rep(0, p), c(-3, 3), rep(3, p)), ncol = num_clusters)
sigma <- array(rep(diag(nrow(mean)), num_clusters),
               dim = c(nrow(mean), nrow(mean), num_clusters))
prob <- c(0.25, 0.5, 0.25)
transformations <- list(
  function(z){
    qpois(pnorm(z), lambda = 7)
  },
  function(z){
    qpois(pnorm(z), lambda = 5)
  }
)
x <- rMixtureNormalTransform(n, mean, sigma, num_clusters, prob, transformations)
plot(x[,1], x[,2])


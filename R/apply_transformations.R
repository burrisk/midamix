applyTransformations <- function(x, transformations, indices = NULL){
  x_mat <- as.matrix(x)
  p <- ncol(x_mat)
  if (p != length(transformations)){
    stop("The number of transformations must be equal to the number of variables.")
  }
  if (is.null(indices)){
    indices <- 1:p
  }
  y_mat <- sapply(indices, function(index){
    transformations[[index]](x_mat[, index])
  })
  y_mat
}

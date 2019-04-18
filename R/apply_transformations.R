#' Apply a list of transformations across a matrix.
#'
#' @param X A matrix.
#' @param transformations A list of functions to apply to corresponding columns.
#' @param indices The columns of the matrix for which to apply the transformation.
#' @return A matrix (or vector) of transformed values.
applyTransformations <- function(X, transformations, indices = NULL){
  x_mat <- as.matrix(X)
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

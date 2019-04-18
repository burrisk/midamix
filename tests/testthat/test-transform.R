context("transform")

vec <- 1:3
row_mat <- matrix(1:3, nrow = 1)
col_mat <- matrix(1:3, ncol = 1)
full_mat <- matrix(1:6, nrow = 2, ncol = 3)
trans1 <- list(function(y){2 * y})
trans3 <- list(function(y){2 * y}, function(y){3 * y}, function(y){4 * y})

test_that("Transforming a vector returns a matrix with one column", {
  expect_equal(ncol(applyTransformations(vec, trans1)), 1)
})

test_that("Simple multiplication transformation", {
  expect_equal(applyTransformations(row_mat, trans3), matrix(c(2, 6, 12), nrow = 1))
  expect_equal(applyTransformations(col_mat, trans1), matrix(c(2, 4, 6), nrow = 3))
  expect_equal(applyTransformations(full_mat, trans3), matrix(c(2, 4, 9, 12, 20, 24), nrow = 2))
})

test_that("Indices subset works properly", {
  expect_equal(applyTransformations(full_mat, trans3, indices = c(1:2)),
               matrix(c(2, 4, 9, 12), nrow = 2))
  expect_equal(applyTransformations(full_mat, trans3, indices = c(2)),
               matrix(c(9, 12), nrow = 2))
})

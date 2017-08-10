context("Test compute_** functions")

test_that("compute_sum_of_list works with and without weights", {
  hold <- matrix(4, ncol = 4, nrow = 4)
  hold <- replicate(4, hold, simplify = FALSE)
  weights1 <- rep(10, 4)
  weights2 <- rep(10, 5)

  m1 <- compute_sum_of_list(.l = hold)
  m2 <- compute_sum_of_list(.l = hold, .w = weights1)
  expect_is(m1, "matrix")
  expect_is(m2, "matrix")

  expect_equal(m1, matrix(16, nrow = 4, ncol = 4))
  expect_equal(m2, matrix(160, nrow = 4, ncol = 4))
  expect_error(compute_sum_of_list(.l = hold, .w = weights2))
})

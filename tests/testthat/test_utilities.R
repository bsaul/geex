context("Geex utilies")

test_that("check_array function works for arrays and vectors", {
  vec <- rep(0, 5)
  arr <- array(rep(0, 21), dim = c(3, 3, 3))

  expect_output(str(check_array(vec)), '[1, 1, 1:5]')

  expect_output(str(check_array(arr)), '[1:3, 1:3, 1:3]')

})

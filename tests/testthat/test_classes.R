context("S4 classes")

test_that("estimation_function S4 class",{
  expect_error(new('estimating_function', .estFUN = function(x, z) {function(y) y + x}))
  expect_error(new('estimating_function', .estFUN = function(data, z) {function(y) y + data}))
  expect_error(new('estimating_function', .estFUN = function(x, data) {function(theta) theta + data}))
  expect_silent(new('estimating_function', .estFUN = function(data, z) {function(theta) theta + data}))
})


test_that("m_estimation_basis S4 class",{
  expect_error(new('m_estimation_basis',
    .estFUN = function(x, z) {function(theta) theta + x},
    .units = "group",
    .data  = geexex))

})

test_that("root_control S4 class",{
  expect_error(new('root_control', .FUN = stats::dnorm))
  expect_error(new('root_control', .options = list(test = 2)))
  expect_silent(new('root_control', .options = list(start = 2)))
})




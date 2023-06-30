context("S4 classes")

test_that("estimation_function S4 class validates correctly", {
  expect_error(new("estimating_function", .estFUN = function(x, z) {
    function(y) y + x
  }))
  # expect_error(new('estimating_function', .estFUN = function(data, z) {function(y) y + data}))
  expect_error(new("estimating_function", .estFUN = function(x, data) {
    function(theta) theta + data
  }))
  expect_silent(new("estimating_function", .estFUN = function(data, z) {
    function(theta) theta + data
  }))
})


test_that("m_estimation_basis S4 class validates correctly", {
  expect_error(new("m_estimation_basis",
    .estFUN = function(x, z) {
      function(theta) theta + x
    },
    .units = "group",
    .data = geexex
  ))

  expect_error(create_basis(
    estFUN = function(x, z) {
      function(theta) theta + x
    },
    units = "group",
    data = geexex
  ))

  expect_silent(create_basis(
    estFUN = function(data, z) {
      function(theta) theta + data
    },
    data = geexex
  ))
})

test_that("root_control S4 class validates correctly", {
  expect_error(new("root_control", .FUN = stats::dnorm))
  expect_error(new("root_control", .options = list(test = 2)))
  expect_silent(new("root_control", .options = list(start = 2)))
})

test_that("correct_control S4 class validates correctly", {
  testFUN1 <- function(components, z) {}
  testFUN2 <- function(components, b) {}
  testFUN3 <- function(component) {}
  expect_error(new("correct_control",
    .FUN = testFUN1,
    .options = list(b = .75)
  ))
  expect_silent(new("correct_control",
    .FUN = testFUN2,
    .options = list(b = .75)
  ))
  expect_error(new("correct_control",
    .FUN = testFUN3
  ))
})

test_that("sandwich_components S4 class validates correctly", {
  expect_error(new("sandwich_components",
    .A = matrix(NA, nrow = 4, ncol = 2)
  ))
  expect_error(new("sandwich_components",
    .A = matrix(NA, nrow = 4, ncol = 4),
    .B = matrix(NA, nrow = 4, ncol = 2)
  ))
  expect_silent(new("sandwich_components",
    .A = matrix(NA, nrow = 4, ncol = 4),
    .B = matrix(NA, nrow = 4, ncol = 4)
  ))
})

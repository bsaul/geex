library(geepack, quietly = TRUE)
data('ohio')
context("Text generic grab() function")
test_that("grab() function return expected objects and fail safely", {
  test_binomial <- glm(resp ~ age, data = ohio,
                       weights = rep(2, nrow(ohio)),
                       family = binomial(link = 'logit'))
  ff <- grab(from = test_binomial, what = "fixed_formula")
  rf <- grab(from = test_binomial, what = "response_formula")
  rp <- grab(from = ohio, what = "response", formula = resp ~ age)
  dm <- grab(from = ohio, what = "design_matrix", rhs_formula = ff)

  expect_is(ff, 'formula')
  expect_is(rf, 'formula')
  expect_is(dm, 'matrix')
  expect_error(  grab(from = ohio, what = "xxx"))

  ee <- grab(from = test_binomial, what = 'psiFUN', data = ohio)
  expect_is(ee, 'function')
})

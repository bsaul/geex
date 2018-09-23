context("Test extraction model eefun utilies for glm objects")
library(geepack, quietly = TRUE)
library(sandwich, quietly = TRUE)
data('ohio')
ohio$row_num <- 1:nrow(ohio)

test_binomial_logit <- glm(resp ~ age, data = ohio,
                     family = binomial(link = 'logit'))


test_gaussian_identity <- glm(Y4 ~ X1 + X2, data = geexex,
                     family = gaussian(link = 'identity'))

glm_eefun <- function(data, model){
  f <- grab_psiFUN(model, data)
  function(theta){
    f(theta)
  }
}

test_that("make_eefun returns functions", {
  expect_is(grab_psiFUN(test_binomial_logit, data = subset(ohio, id == 1)),
            'function')
})

# Removing weight functionality for time being, until it can be more robust
# test_that("make_eefun returns error when length of weight vector is not equal to 1 or # in cluster", {
#   expect_error(grab_psiFUN(test_binomial, data = subset(ohio, id == 1), weights = c(2, 2)))
# })

test_that("estimate equations obtains correct values for parameters and standard errors for logit link", {
  x <- m_estimate(glm_eefun,
                  data = ohio,
                  units = 'id',
                  root_control  = new('root_control', .options = list(start = c(-1.7, -.11))),
                  outer_args = list(model = test_binomial_logit))
  expect_equivalent(x@estimates, coef(test_binomial_logit))

  # Check nobs()
  expect_equal(nobs(x), 537)

  # Form sandwich estimator "by hand" with the help of sandwich
  psi <- apply(estfun(test_binomial_logit), 2, function(x) tapply(x, test_binomial_logit$data[['id']], sum))
  n <- length(unique(test_binomial_logit$data[['id']]))
  mmeat <- crossprod(as.matrix(psi))/n
  bbread <- vcov(test_binomial_logit) * n
  Sig <- (bbread %*% mmeat %*% t(bbread))/n

  expect_equal(Sig, x@vcov, tolerance = 1e-5, check.attributes = FALSE)
})

test_that("estimate equations obtains correct values for parameters and standard errors for identity link", {
  x <- m_estimate(glm_eefun,
                  data = geexex,
                  root_control  = new('root_control', .options = list(start = c(0, 0, 0))),
                  outer_args = list(model = test_gaussian_identity))
  expect_equivalent(x@estimates, coef(test_gaussian_identity))

  # Compare to sandwich
  Sig <- sandwich(test_gaussian_identity)

  expect_equal(Sig, x@vcov, tolerance = 1e-5, check.attributes = FALSE)
})


test_that("estimate equations obtains correct values for parameters and standard errors under rowwise conditions", {
  x <- m_estimate(glm_eefun,
                  data = ohio,
                  units = 'row_num',
                  root_control  = new('root_control', .options = list(start = c(-1.7, -.11))),
                  outer_args = list(model = test_binomial_logit))
  expect_equivalent(x@estimates, coef(test_binomial_logit))

  x2 <- m_estimate(glm_eefun,
                  data = ohio,
                  # units = 'row_num',
                  root_control  = new('root_control', .options = list(start = c(-1.7, -.11))),
                  outer_args = list(model = test_binomial_logit))
  expect_equivalent(x2@estimates, coef(test_binomial_logit))

  # Form sandwich estimator "by hand" with the help of sandwich
  psi <- apply(estfun(test_binomial_logit), 2, function(x) tapply(x, test_binomial_logit$data[['row_num']], sum))
  n <- length(unique(test_binomial_logit$data[['row_num']]))
  mmeat <- crossprod(as.matrix(psi))/n
  bbread <- vcov(test_binomial_logit) * n
  Sig <- (bbread %*% mmeat %*% t(bbread))/n

  expect_equal(Sig, x@vcov, tolerance = 1e-5, check.attributes = FALSE)
  expect_equal(Sig, x2@vcov, tolerance = 1e-5, check.attributes = FALSE)
})

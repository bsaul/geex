context("Test extraction model eefun utilies for glm objects")
library(geepack, quietly = TRUE)
library(sandwich, quietly = TRUE)
data('ohio')

test_binomial <- glm(resp ~ age, data = ohio,
                     weights = rep(2, nrow(ohio)),
                     family = binomial(link = 'logit'))

glm_eefun <- function(data, model){
  make_eefun(model, data, weights = 2)
}

test_that("make_eefun returns functions", {
  expect_is(make_eefun(test_binomial, data = subset(ohio, id == 1), weights = 2),
            'function')
})

test_that("make_eefun returns error when length of weight vector is not equal to 1 or # in cluster", {
  expect_error(make_eefun(test_binomial, data = subset(ohio, id == 1), weights = c(2, 2)))
})

test_that("estimate equations obtains correct values for parameters and standard errors ", {
  x <- estimate_equations(glm_eefun,
                          data = ohio,
                          units = 'id',
                          rootFUN_control  = list(start = c(-1.7, -.11)),
                          outer_eeargs = list(model = test_binomial))
  expect_equivalent(x$estimates, coef(test_binomial))

  # Form sandwich estimator "by hand" with the help of sandwich
  psi <- apply(estfun(test_binomial), 2, function(x) tapply(x, test_binomial$data[['id']], sum))
  n <- length(unique(test_binomial$data[['id']]))
  mmeat <- crossprod(as.matrix(psi))/n
  bbread <- vcov(test_binomial) * n
  Sig <- (bbread %*% mmeat %*% t(bbread))/n

  expect_equal(Sig, x$vcov, tolerance = 1e-5, check.attributes = FALSE)
})

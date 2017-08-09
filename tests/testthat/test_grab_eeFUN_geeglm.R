context("Test extraction model eefun utilies for geepack objects")
library(geepack, quietly = TRUE)
data('ohio')

test_binomial <- geeglm(resp ~ age, id = id, data = ohio,
                        family = binomial(link = 'logit'))

gee_eefun <- function(data, model){
  f <- grab_psiFUN(model, data)
  function(theta){
    f(theta)
  }
}

test_that("grab_psiFUN returns functions", {
  expect_is(grab_psiFUN(test_binomial, data = subset(ohio, id == 1)),
            'function')
})

test_that("estimate equations obtains correct values for parameters and standard errors ", {
  x <- m_estimate(gee_eefun,
                          data = ohio,
                          units = 'id',
                          root_control = new('root_control', .options = list(start = coef(test_binomial))),
                          outer_args  = list(model = test_binomial))
  expect_equal(x@estimates, coef(test_binomial))
  expect_equal(sqrt(diag(x@vcov)), summary(test_binomial)$coefficients[, 2])
})

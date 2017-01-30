context("Test extraction model eefun utilies for geepack objects")
library(geepack, quietly = TRUE)
data('ohio')

test_binomial <- geeglm(resp ~ age, id = id, data = ohio,
                        family = binomial(link = 'logit'))

gee_eefun <- function(data, model){
  make_eefun(model, data)
}

test_that("make_eefun returns functions", {
  expect_is(make_eefun(test_binomial, data = subset(ohio, id == 1)),
            'function')
})

test_that("estimate equations obtains correct values for parameters and standard errors ", {
  x <- estimate_equations(gee_eefun,
                          data = ohio,
                          units = 'id',
                          roots = coef(test_binomial),
                          model = test_binomial)
  expect_equal(x$parameters, coef(test_binomial))
  expect_equal(sqrt(diag(x$vcov)), summary(test_binomial)$coefficients[, 2])
})

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
                          rootFUN_control = list(start = coef(test_binomial)),
                          outer_eeargs  = list(model = test_binomial))
  expect_equal(x$estimates, coef(test_binomial))
  expect_equal(sqrt(diag(x$vcov)), summary(test_binomial)$coefficients[, 2])
})

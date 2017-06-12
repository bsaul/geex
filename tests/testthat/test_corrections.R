context("Finite sample corrections")

set.seed(1)
n  <- 30
mu <- 5
sigma <- 2
dt <- data.frame(Y = rnorm(n, mean = mu, sd = sigma), id = 1:n)

test_eefun <- function(data){
  function(theta){
    with(data,
         c(Y - theta[1],
           (Y - theta[1])^2 - theta[2] )
    )
  }
}

bias_estimates <- estimate_equations(
  eeFUN = test_eefun,
  data  = dt, units = 'id',
  roots = c(1,1),
  corrections_list = list(test = list(fun = fay_bias_correction, options = list(b = 0.75))))

df1_estimates <- estimate_equations(
  eeFUN = test_eefun,
  data  = dt, units = 'id',
  roots = c(1,1),
  corrections_list = list(test = list(fun = fay_df_correction, options = list(b = 0.75, L = c(1, 1), version = 1))))

df2_estimates <- estimate_equations(
  eeFUN = test_eefun,
  data  = dt, units = 'id',
  roots = c(1,1),
  corrections_list = list(test = list(fun = fay_df_correction, options = list(b = 0.75, L = c(1, 1), version = 2))))

test_that("Bias correction returns matrix", {
  expect_is(bias_estimates$vcov, 'matrix')
})

test_that("DF1 correction returns matrix", {
  expect_is(df1_estimates$vcov, 'matrix')
})

test_that("DF2 correction returns matrix", {
  expect_is(df1_estimates$vcov, 'matrix')
})

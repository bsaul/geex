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

rooter <- new('root_control', .options = list(start = c(1, 1)))

bias_estimates <- m_estimate(
  estFUN = test_eefun,
  data  = dt,
  units = 'id',
  root_control = rooter,
  corrections = list(test = correction(fay_bias_correction,
                                            list(b = 0.75))))

df1_estimates <- m_estimate(
  estFUN = test_eefun,
  data  = dt,
  units = 'id',
  root_control = rooter,
  corrections = list(test = correction(correctFUN = fay_df_correction,
                            correctFUN_options = list(b = 0.75, L = c(1, 1), version = 1))))

df2_estimates <- m_estimate(
  estFUN = test_eefun,
  data  = dt,
  units = 'id',
  root_control = rooter,
  corrections = list(test = correction(correctFUN = fay_df_correction,
                                      correctFUN_options = list(b = 0.75, L = c(1, 1), version = 2))))



test_that("Bias correction returns matrix", {
  expect_is(bias_estimates@corrections$test, 'matrix')
})

test_that("DF1 correction returns a scalar", {
  expect_is(df1_estimates@corrections$test, 'numeric')
  expect_equal(length(df1_estimates@corrections$test), 1)
})

test_that("DF2 correction returns a scalar", {
  expect_is(df2_estimates@corrections$test, 'numeric')
  expect_equal(length(df2_estimates@corrections$test), 1)
})

test_that("get_corrections() accessor returns list", {
  expect_is(get_corrections(bias_estimates), 'list')
})



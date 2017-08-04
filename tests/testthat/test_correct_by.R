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
  corrections_list = list(test = list(correctFUN = correct_by_fay_bias,
                                      correctFUN_control = list(b = 0.75))))

df1_estimates <- m_estimate(
  estFUN = test_eefun,
  data  = dt,
  units = 'id',
  root_control = rooter,
  corrections_list = list(test = list(correctFUN = correct_by_fay_df,
                                      correctFUN_control = list(b = 0.75, L = c(1, 1), version = 1))))

df2_estimates <- m_estimate(
  estFUN = test_eefun,
  data  = dt,
  units = 'id',
  root_control = rooter,
  corrections_list = list(test = list(correctFUN = correct_by_fay_df,
                                      correctFUN_control = list(b = 0.75, L = c(1, 1), version = 2))))



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

test_that("check_corrections picks up missing correctFUN", {
  correction_tester <- list(test1 = list(correctFUN = correct_by_fay_df),
                            test2 = list(fun = correct_by_fay_df))
  expect_error(check_corrections(correction_tester))
})

test_that("check_corrections picks up additional arguments", {
  correction_tester <- list(test1 = list(correctFUN = correct_by_fay_df,
                                         correctFUN_control = list(x = 2),
                                         errormaker = 2))
  expect_warning(check_corrections(correction_tester))
})

test_that("check_corrections does not throw error when correction list is correct", {
  correction_tester <- list(test1 = list(correctFUN = correct_by_fay_df,
                                         correctFUN_control = list(x = 1)),
                            test2 = list(correctFUN = correct_by_fay_bias,
                                         correctFUN_control = list(x = 1)))
  expect_silent(check_corrections(correction_tester))
})

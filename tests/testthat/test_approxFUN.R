context("approxFUN computations")

quantile_eefun <- function(data){
  function(theta){
    0.5  - (data$Y1 <= theta[1])
  }
}

spline_approx <- function(psi, eval_theta){
  ### Use splinefun ####
  y <- Vectorize(psi)(eval_theta)
  function(theta) splinefun(x = eval_theta, y = y)(theta)
}

test_that("approxFUN works", {
  x <- estimate_equations(
    eeFUN = quantile_eefun,
    data  = geexex,
    roots = 4.7,
    approxFUN = spline_approx,
    approxFUN_control = list(eval_theta = seq(3, 6, by = .4))
  )
  expect(!is.null(x$parameters), 'Estimate_equations did not return parameters using approxFUN')
  expect(!is.null(x$vcov), 'Estimate_equations did not return vcov using approxFUN')
})





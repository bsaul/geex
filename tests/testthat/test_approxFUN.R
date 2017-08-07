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

rooter <- new('root_control', .options = list(start = c(1)))
approx <- new('approx_control',
              .FUN = spline_approx,
              .options = list(eval_theta = seq(3, 6, by = .4)))
gbasis <- new('m_estimation_basis',
              .estFUN = quantile_eefun,
              .data   = geexex)


test_that("approx_control is working", {
  expect_silent({psii <- create_psiFUN_list(.basis = gbasis,
                                    .approx_control = approx)})

  x <- m_estimate(
    estFUN = quantile_eefun,
    data  = geexex,
    root_control = rooter,
    approx_control = approx )

  expect(!is.null(x@estimates), 'Estimate_equations did not return parameters using approxFUN')
  expect(!is.null(x@vcov), 'Estimate_equations did not return vcov using approxFUN')
})





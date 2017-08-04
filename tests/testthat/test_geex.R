context("Geex computations")


test_eefun1 <- function(data){
  function(theta){
    with(data,
         c(Y1 - theta[1],
           (Y1 - theta[1])^2 - theta[2] )
    )
  }
}

#This is nonsensical used for testing the units= option
test_eefun2 <- function(data){
  function(theta){
    with(data,
         c(sum(Y1)- theta[1],
           (sum(Y1) - theta[1])^2 - theta[2] )
    )
  }
}
n <- nrow(geexex)
theta_hat <- c(mean(geexex$Y1), var(geexex$Y1) * (n - 1) / n)

gbasis_good <- new('m_estimation_basis',
               .estFUN = test_eefun1,
               .data   = geexex)

gbasis_bad <- new('m_estimation_basis',
                   .estFUN = test_eefun2,
                   .data   = geexex)

test_that("estimate_GFUN_roots is working", {
  expect_silent({psii <- create_psi(.basis = gbasis_good)})
  expect_is(psii, 'list')
  expect_silent({gtest <- create_GFUN(.psi_list  = psii, .basis = gbasis_good)})
  expect_is(gtest, 'function')

  # Check running estimate_GFUN_roots without root_control gives error
  expect_error(estimate_GFUN_roots(.GFUN = gtest))

  rooter_test <- new("root_control", .options = list(start = c(3, 3)))
  # Check running estimate_GFUN_roots with proper root_control does not give error
  expect_silent({root_test <- estimate_GFUN_roots(.GFUN = gtest,
                                                 .root_control = rooter_test )})

  # Check that the first elements of the root is at least close to the
  # closed form
  expect_equal(root_test$root[1], theta_hat[1], tolerance = 1e-4)
})

test_that("estimate_sandwich_matrices working", {
  expect_silent({psii <- create_psi(.basis = gbasis_good)})
  expect_is(psii, 'list')

  # Check running estimate_sandwich_matrices without necessary arguments gives error
  expect_error(estimate_sandwich_matrices(.psi_list = psii))

  # Check running estimate_GFUN_roots with proper root_control does not give error
  expect_silent({mat_test <- estimate_sandwich_matrices(.psi_list = psii,
                                                        .basis    = gbasis_good,
                                                        .theta    = theta_hat)})

  expect_is(mat_test, 'list')
})


test_that("m_estimation computations are working", {


  rooter_test <- new("root_control", .options = list(start = c(3, 3)))

  # should give error if units is not in the data.frame
  expect_error(m_estimate(estFUN  = test_eefun1,
                            data  = geexex, units = 'id',
                            root_control = rooter_test))

  estimates_2 <- m_estimate(estFUN = test_eefun1,
                            data  = geexex,
                            root_control = rooter_test)

  expect_is(estimates_2, 'geex')

})

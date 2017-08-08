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

# For testing inner and outer args
test_eefun3 <- function(data, psi){
  function(theta, alpha){
    with(data,
         c(Y1 - theta[1]/psi,
           (Y1 - theta[1])^2 - theta[2]*alpha)
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
  expect_is(gbasis_good, 'm_estimation_basis')
  expect_silent({gtest <- create_GFUN(.basis = gbasis_good)})
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
  # Check running estimate_sandwich_matrices without necessary arguments gives error
  expect_error(estimate_sandwich_matrices(.basis = gbasis_good))

  # Check running estimate_GFUN_roots with proper root_control does not give error
  expect_silent({mat_test <- estimate_sandwich_matrices(.basis = gbasis_good,
                                                        .theta    = theta_hat)})

  expect_is(mat_test, 'sandwich_components')
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


test_that("m_estimation works when given outer and/or inner args", {

  # should work
  expect_is(m_estimate(estFUN  = test_eefun3,
                          data  = geexex,
                          outer_args = list(psi = 2),
                          inner_args = list(alpha = 3),
                          root_control = setup_root_solver(start = c(3, 3))),
            'geex')

  # should fail
  expect_error(m_estimate(estFUN  = test_eefun3,
                       data  = geexex,
                       outer_args = list(psi = 2),
                       root_control = setup_root_solver(start = c(3, 3))))
  # should fail
  expect_error(m_estimate(estFUN  = test_eefun3,
                       data  = geexex,
                       inner_args = list(alpha = 3),
                       root_control = setup_root_solver(start = c(3, 3))))

})

test_that("geex class accessors work", {

  # should work
  expect_is({hold <- m_estimate(estFUN  = test_eefun1,
                       data  = geexex,
                       root_control = setup_root_solver(start = c(3, 3)))},
            'geex')

  expect_is(vcov(hold), 'matrix')
  expect_is(coef(hold), 'numeric')
  expect_is(roots(hold), 'numeric')
})

test_that("m_estimate() runs when compute_roots = FALSE", {

  expect_error(m_estimate(estFUN  = test_eefun1,
                                data  = geexex,
                                compute_roots = FALSE))

  expect_is(m_estimate(estFUN  = test_eefun1,
                       data  = geexex,
                       roots = c(1, 1),
                       compute_roots = FALSE), 'geex')
})

test_that("m_estimate() runs when compute_vcov = FALSE", {
  hold <- m_estimate(estFUN  = test_eefun1,
                     data  = geexex,
                     compute_vcov = FALSE,
                     root_control = setup_root_solver(start = c(3, 3)))

  expect_is(hold, 'geex')
  expect_is(vcov(hold), 'matrix')
  expect_equal(dim(vcov(hold)), c(0, 0))
})

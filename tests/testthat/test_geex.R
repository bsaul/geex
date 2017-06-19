context("Geex computations")

test_that("Basic computations are working", {
  set.seed(1)
  n  <- 30
  mu <- 5
  sigma <- 2
  dt <- data.frame(Y = rnorm(n, mean = mu, sd = sigma), id = 1:n, id2 = rep(1:(n/2), each = 2))

  test_eefun1 <- function(data){
    function(theta){
      with(data,
           c(Y - theta[1],
             (Y - theta[1])^2 - theta[2] )
      )
    }
  }

  #This is nonsensical used for testing the units= option
  test_eefun2 <- function(data){
    function(theta){
      with(data,
           c(sum(Y)- theta[1],
             (sum(Y) - theta[1])^2 - theta[2] )
      )
    }
  }

  glist <- list(eeFUN = test_eefun1, splitdt = split(dt, dt$id))
  theta_hat <- c(mean(dt$Y), var(dt$Y) * (n - 1) / n)
  mats      <- compute_matrices(glist,
                                theta = theta_hat)
  roots     <- eeroot(glist, start = c(1, 1))

  estimates_1 <- estimate_equations(eeFUN = test_eefun1,
                                  data  = dt, units = 'id',
                                  roots = c(1,1))

  estimates_2 <- estimate_equations(eeFUN = test_eefun1,
                                  data  = dt,
                                  roots = c(1,1))

  estimates_3 <- estimate_equations(eeFUN = test_eefun2,
                                    data  = dt,
                                    units = 'id2',
                                    roots = c(1,1))
  # Check running eeroot without start gives error
  expect_error(eeroot(glist))

  # Check mats is list
  expect_is(mats, 'list')

  # Check length of mats
  expect_length(mats, 4)

  # Check accuracy of roots
  expect_equal(theta_hat, roots$root)

  # Check the estimates
  expect_equal(estimates_1$parameters, c(theta_hat))

  # Check that estimate_equations gives the same answer when units left NULL
  expect_equal(estimates_1, estimates_2)

  # This should be obvious but just checking
  expect_false(isTRUE(all.equal(estimates_1, estimates_3)))

})

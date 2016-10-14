context("Geex computations")

test_that("Basic computations are working", {
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
  glist <- list(eeFUN = test_eefun, splitdt = split(dt, dt$id))
  theta_hat <- c(mean(dt$Y), var(dt$Y) * (n - 1) / n)
  mats      <- compute_matrices(glist,
                                theta = theta_hat)

  roots     <- eeroot(glist, start = c(1, 1))
  estimates <- estimate_equations(eeFUN = test_eefun,
                                  data  = dt, units = 'id',
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
  expect_equal(estimates$parameters, c(theta_hat))
})

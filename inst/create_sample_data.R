#------------------------------------------------------------------------------#
# Create the example data set for geex: geexex
#------------------------------------------------------------------------------#
set.seed(100)
n <- 100
# Oracle parms
alpha <- 2
beta  <- 3
gamma <- 2
delta <- 1.5

e1 <- e2 <- e3 <- rnorm(n)
sigma_e <- 1
sigma_U <- .25
sigma_tau <- 1

tau <- c(0.1, 0.1, 0.5)

geexex <- as.data.frame(dplyr::data_frame(
  # Used in Examples 1-3, 5, 6, 7
  Y1 = rnorm(n, mean = 5, sd = 4),
  Y2 = rnorm(n, mean = 2, sd = 1),

  # Used in Example 4
  X1 = rgamma(n, shape = 5),
  Y3 = alpha + (beta * X1) + (sigma_e * e1),
  W1 = X1 + (sigma_U * e2),
  Z1 = gamma + (delta * X1) + (sigma_tau * e3),

  # Used in Example 8
  X2 = rep(0:1, each = n/2),
  Y4 = as.numeric(cbind(1, X1, X2) %*% tau) + e1,

  # Used in Example 9
  Y5 = rbinom(n, 1, prob = as.numeric(plogis(cbind(1, X1, X2) %*% tau)))
))


save(geexex, file = 'data/geexex.rda')
rm(geexex, n, alpha, beta, gamma, delta, e1, e2, e3, sigma_e, sigma_U, sigma_tau, tau)

## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
library(geex)
library(knitr)
opts_knit$set(progress = TRUE, verbose = TRUE)

## ----SB1_setup, echo=FALSE-----------------------------------------------
n  <- 100
mu <- 5
sigma <- 2
dt <- data.frame(Y = rnorm(n, mean = mu, sd = sigma), id = 1:n)

## ----SB1_eefun, echo=TRUE, results='hide'--------------------------------
SB1_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
       (Y - theta[1])^2 - theta[2] )
    )
  }
}

## ----SB1_run, echo=TRUE--------------------------------------------------
estimates <- estimate_equations(
  eeFUN = SB1_eefun, 
  data  = dt,
  rootFUN_control = list(start = c(1,1)))

## ----SB1_clsform, echo=FALSE---------------------------------------------
## Compare to closed form ##

A <- diag(1, nrow = 2)

B <- with(dt, {
  Ybar <- mean(Y)
  B11 <- mean( (Y - Ybar)^2 )
  B12 <- mean( (Y - Ybar) * ((Y - Ybar)^2 - B11) )
  B22 <- mean( ((Y - Ybar)^2 - B11)^2 )
  matrix(
    c(B11, B12,
      B12, B22), nrow = 2
  )
})

## closed form roots
# note that var() divides by n - 1, not n
theta_cls <- dplyr::summarize(dt, p1 = mean(Y), p2 = var(Y) * (n() - 1)/ n() )

# closed form
Sigma_cls <- (solve(A) %*% B %*% t(solve(A))) / n

## ----SB1_results, echo = FALSE-------------------------------------------
results <- list(geex = estimates[c('estimates', 'vcov')], 
                cls = list(parameters = theta_cls, vcov = Sigma_cls))
results

## ----SB2_setup, echo=FALSE-----------------------------------------------
n  <- 100
muY <- 5
sigmaY <- 2
muX <- 2
sigmaX <- 0.2
dt <- data.frame(Y  = rnorm(n, mean = muY, sd = sigmaY), 
                 X  = rnorm(n, mean = muX, sd = sigmaX),
                 id = 1:n)

## ----SB2_eefun, echo = TRUE----------------------------------------------
SB2_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
        X - theta[2],
        theta[1] - (theta[3] * theta[2]) )
    )
  }
}

## ----SB2_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(
  eeFUN = SB2_eefun, 
  data  = dt, 
  rootFUN_control = list(start = c(1, 1, 1)))

## ----SB2_clsform, echo = FALSE-------------------------------------------
## Compare to closed form ##

A <- with(dt, {
 matrix(
  c(1 , 0, 0,
    0 , 1, 0,
    -1, mean(Y)/mean(X), mean(X)),
  byrow = TRUE, nrow = 3)
})

B <- with(dt, {
  matrix(
    c(var(Y)   , cov(Y, X), 0,
      cov(Y, X), var(X)   , 0,
      0, 0, 0),
    byrow = TRUE, nrow = 3)
})

## closed form roots
theta_cls <- dplyr::summarize(dt, p1 = mean(Y), p2 = mean(X), p3 = p1/p2)

## closed form covariance
Sigma_cls <- (solve(A) %*% B %*% t(solve(A))) / n

## ----SB2_results, echo = FALSE-------------------------------------------
results <- list(geex = estimates[c('estimates', 'vcov')], 
                cls = list(parameters = theta_cls, vcov = Sigma_cls))
results

## ----SB3_setup, echo=FALSE-----------------------------------------------
n  <- 100
mu <- 5
sigma <- 4
set.seed(100) # running into issue where sqrt(theta2) and log(theta2) return NaN for some seeds
dt <- data.frame(Y  = rnorm(n, mean = mu, sd = sigma), 
                 id = 1:n)

## ----SB3_eefun, echo = TRUE----------------------------------------------
SB3_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
       (Y - theta[1])^2 - theta[2],
       sqrt(theta[2]) - theta[3],
       log(theta[2]) - theta[4])
    )
  }
}

## ----SB3_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(
  eeFUN= SB3_eefun, 
  data  = dt,
  rootFUN_control = list(start = c(1, 1, 1, 1)))

## ----SB3_clsform, echo = FALSE-------------------------------------------
## closed form roots
theta_cls <- dplyr::summarize(dt, p1 = mean(Y), p2 = sum((Y - p1)^2)/n(), p3 = sqrt(p2), p4 = log(p2))

## Compare to closed form ##
theta2 <- theta_cls$p2
mu3 <- moments::moment(dt$Y, order = 3, central = TRUE)
mu4 <- moments::moment(dt$Y, order = 4, central = TRUE)
# A <- matrix(c(1, 0, 0, 0,
#               0, 1, 0, 0,
#               0, -1/(2 * sqrt(theta2)), 1, 0,
#               0, -1/theta2, 0, 1), 
#             byrow = TRUE, nrow = 4)
# B <- matrix(c(1/theta2, mu3/(2 * theta2^3), 0, 0,
#               mu3/(2 * theta2^3), (mu4 - theta2^2)/(4 * theta2^4), 0, 0,
#               0, 0, 0, 0,
#               0, 0, 0, 0),
#             byrow = TRUE, nrow = 4)

## closed form covariance
Sigma_cls <- matrix(
  c(theta2, mu3, mu3/(2*sqrt(theta2)), mu3/theta2,
    mu3, mu4 - theta2^2, (mu4 - theta2^2)/(2*sqrt(theta2)), (mu4 - theta2^2)/theta2,
    mu3/(2 * sqrt(theta2)), (mu4 - theta2^2)/(2*sqrt(theta2)), (mu4 - theta2^2)/(4*theta2), (mu4 - theta2^2)/(2*theta2^(3/2)),
    mu3/theta2, (mu4 - theta2^2)/theta2, (mu4 - theta2^2)/(2*theta2^(3/2)), (mu4/theta2^2) - 1) ,
  nrow = 4, byrow = TRUE) / n
## closed form covariance
# Sigma_cls <- (solve(A) %*% B %*% t(solve(A))) / n

## ----SB3_results, echo = FALSE-------------------------------------------
results <- list(geex = estimates[c('estimates', 'vcov')], 
                cls = list(parameters = theta_cls, vcov = Sigma_cls))
results


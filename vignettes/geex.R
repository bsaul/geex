## ---- echo = TRUE, message = FALSE, warning=FALSE------------------------
library(eex)
library(dplyr)
library(inferference)
library(sandwich)
# library(microbenchmark)

## ----eefun, echo=TRUE----------------------------------------------------
eefun <- function(data, model){
  X <- model.matrix(model, data = data)
  Y <- model.response(model.frame(model, data = data))
  function(theta){
    lp  <- X %*% theta
    rho <- plogis(lp)

    score_eqns <- apply(X, 2, function(x) sum((Y - rho) * x))
    score_eqns
  }
}

## ----example1------------------------------------------------------------
vaccinesim$ID <- 1:nrow(vaccinesim)
mglm    <- glm(A ~ X1, data = vaccinesim, family = binomial)
split_data  <- split(vaccinesim, vaccinesim$ID)
# The list needed for the compute_matrices
# For now, theta needs to be passed since geex does not do point estimates yet
example <- list(eeFUN = eefun, splitdt = split_data)
root <- eeroot(example, start = c(-.35, 0), model = mglm )$root

mats <- compute_matrices(obj = example, model = mglm,
                         theta = coef(mglm),
                         numDeriv_options = list(method = 'Richardson'))
# Compare point estimates
root # from GEEX
coef(mglm) # from the GLM function

# Compare variance estimates
compute_sigma(mats)
sandwich::sandwich(mglm)

## ----eefun2, echo=TRUE---------------------------------------------------
eefun2 <- function(data, model, alpha){
  X <- model.matrix(model, data = data)
  A <- model.response(model.frame(model, data = data))
  
  function(theta){
    p  <- length(theta)
    p1 <- length(coef(model))
    lp  <- X %*% theta[1:p1]
    rho <- plogis(lp)

    hh  <- (rho/alpha)^A * ((1-rho)/(1-alpha))^(1 - A)
    IPW <- 1/(exp(sum(log(hh))))

    score_eqns <- apply(X, 2, function(x) sum((A - rho) * x))
    with(data, {
      ce0 <- mean(Y * (A == 0)) * IPW
      ce1 <- mean(Y * (A == 1)) * IPW
      
      c(score_eqns,
        ce0 - theta[p - 1],
        ce1 - theta[p])
    })
  }
}

## ----example2, echo =TRUE------------------------------------------------
test <- interference(Y | A ~ X1 | group, 
                     data = vaccinesim,
                     model_method = 'glm',
                     allocations = c(.35, .4))

mglm        <- glm(A ~ X1, data = vaccinesim, family = binomial)
split_data  <- split(vaccinesim, vaccinesim$group)

# The list needed for the compute_matrices
example <- list(eeFUN   = eefun2, 
                splitdt = split_data)

root <- eeroot(example, start =c(coef(mglm), .4,  .13), model = mglm, alpha = .35 )$root
                # Plugging in the values of CE_hat(0) and CE_hat(1)
                # since geex does not currently do point estimates)

mats <- compute_matrices(obj = example, 
                         numDeriv_options = list(method = 'Richardson'),
                         theta   = c(coef(mglm), 0.42186669,  0.15507946),
                # Plugging in the values of CE_hat(0) and CE_hat(1)
                # since geex returns incorrect point estimates for this problem at the moment
                         model = mglm, 
                         alpha = .35)
# conpare SE estimates
L <- c(0, 0, -1, 1)
Sigma <- compute_sigma(mats)
sqrt(t(L) %*% Sigma %*% L)  # from GEEX
direct_effect(test, allocation = .35)$std.error # from inferference

## ----  echo = TRUE-------------------------------------------------------
### Point Estimate via inferference
c(coef(mglm), 0.42186669,  0.15507946)

### Points Estimates with GEEX
root

## ---- echo = TRUE--------------------------------------------------------
eeroot(example, start =c(coef(mglm), 0.42186669,  0.15507946), model = mglm, alpha = .35 )$root

## ----SB_example1, echo=TRUE----------------------------------------------
n  <- 100
dt <- data.frame(Y = rnorm(n, mean = 5, sd = 4), id = 1:n)
split_data <- split(dt, dt$id)

SB_ex1_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
       (Y - theta[1])^2 - theta[2] )
    )
  }
}

example <- list(eeFUN   = SB_ex1_eefun, 
                splitdt = split_data)

root <- eeroot(example, start = c(1,1))$root

mats  <- compute_matrices(obj = example,
                          theta = root,
                          numDeriv_options = list(method = 'Richardson'))
Sigma <- compute_sigma(mats)

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

## geex roots 
root

## closed form roots
# note that var() divides by n - 1, not n
summarize(dt, p1 = mean(Y), p2 = var(Y))

# geex variance estimate
Sigma

# closed form
(solve(A) %*% B %*% t(solve(A))) / n


## ----SB_example2, echo=TRUE----------------------------------------------
n  <- 100
dt <- data.frame(Y  = rnorm(n, mean = 5, sd = 4), 
                 X  = rnorm(n, mean = 2, sd = .09),
                 id = 1:n)
split_data <- split(dt, dt$id)

SB_ex2_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
        X - theta[2],
        theta[1] - (theta[3] * theta[2]) )
    )
  }
}

example <- list(eeFUN   = SB_ex2_eefun, 
                splitdt = split_data)
root <- eeroot(obj = example, start = c(1, 1, 1))$root
mats  <- compute_matrices(obj = example,
                          theta = root,
                          numDeriv_options = list(method = 'Richardson'))
Sigma <- compute_sigma(mats)

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

## geex roots 
root

## closed form roots
summarize(dt, p1 = mean(Y), p2 = mean(X), p3 = p1/p2)

# geex variance estimate
Sigma

# closed form
(solve(A) %*% B %*% t(solve(A))) / n


## ----SB_example3, echo=TRUE----------------------------------------------
set.seed(100) # running into issue where sqrt(theta2) and log(theta2) return NaN for some seeds
n  <- 100
dt <- data.frame(Y  = rnorm(n, mean = 5, sd = 4), 
                 id = 1:n)
split_data <- split(dt, dt$id)

SB_ex3_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
       (Y - theta[1])^2 - theta[2],
       sqrt(theta[2]) - theta[3],
       log(theta[2]) - theta[4])
    )
  }
}

example <- list(eeFUN   = SB_ex3_eefun, 
                splitdt = split_data)
root <- eeroot(obj = example, start = c(1, 1, 1, 1))$root


mats  <- compute_matrices(obj = example,
                          theta = root,
                          numDeriv_options = list(method = 'Richardson'))
Sigma <- compute_sigma(mats)

## geex roots 
root

## closed form roots
summarize(dt, p1 = mean(Y), p2 = sum((Y - p1)^2)/n(), p3 = sqrt(p2), p4 = log(p2))

# geex variance estimate
Sigma



## ----SB_example4, echo=TRUE----------------------------------------------
n  <- 10000

# Oracle parms
alpha <- 2
beta  <- 3
gamma <- 2
delta <- 1.5
e1 <- e2 <- e3 <- rnorm(n)
sigma_e <- 1
sigma_U <- .25
sigma_tau <- 1
### Random variables

X <- rgamma(n, shape = 5)
X <- rnorm(n, sd = 1)

dt <- data.frame(Y  = alpha + (beta * X) + (sigma_e * e1), 
                 W  = X + (sigma_U * e2),
                 T_  = gamma + (delta * X) + (sigma_tau * e3),
                 id = 1:n)
split_data <- split(dt, dt$id)

SB_ex4_eefun <- function(data){
  function(theta){
    with(data,
      c(theta[1] - T_,
        theta[2] - W,
        (Y - (theta[3] * W)) * (theta[2] - W),
        (Y - (theta[4] * W)) * (theta[1] - T_))
    )
  }
}



example <- list(eeFUN   = SB_ex4_eefun, 
                splitdt = split_data)
root <- eeroot(obj = example, start = c(1, 1, 1, 1))$root
root

## compare to closed form
c(theta1 = mean(dt$T_),
  theta2 = mean(dt$W),
  theta3 = coef(lm(Y ~ W, data = dt))[2],
  theta4 = coef(lm(Y ~ T_, data = dt))[2]/coef(lm(W ~ T_, data = dt))[2])


mats  <- compute_matrices(obj = example,
                          theta = root,
                          numDeriv_options = list(method = 'Richardson'))
Sigma <- compute_sigma(mats)

Sigma


## compare to closed form
# TODO


## ----SB_example5, echo=TRUE----------------------------------------------

n <- 100
theta0 <- 0
dt <- data.frame(X = rnorm(n, mean = 2),
                 id = 1:n)
split_data <- split(dt, dt$id)

distr <- function(x, theta0){
  FF <- ecdf(x - theta0)
  approxfun(x - theta0, y = FF(x - theta0), method = "linear",
            0, 1, rule = 1, f = 0, ties = mean)
}

FF <- distr(dt$X, 0)
dens <- density(dt$X)
ff2 <- approxfun(dens$x, dens$y, yleft = 0, yright = 0)
integrand <- function(y){
  ff2(y)^2
}

IC_denom <- integrate(integrand, lower = min(dt$X), upper = max(dt$X))$value

SB_ex5_eefun <- function(data, theta0 = 0){
  Xi <- data$X
  function(theta){
     IC_HL <- (FF(Xi - theta0) - 0.5)/IC_denom
     c(IC_HL - (theta[1] - theta0),
       Xi - theta[2])
  }
}

split_data <- split(dt, dt$id)

example <- list(eeFUN   = SB_ex5_eefun, 
                splitdt = split_data)

root <- eeroot(obj = example, start = c(2.2, 1))$root
root

X <- dt$X
pair_means <- numeric(length(dt$X) - 1)
for(i in 1:(length(X) - 1)){
 pair_means[i] <-  (X[i] + X[i + 1])/2
}

c(median(pair_means), mean(dt$X))

# Asymp correlation btn HL estimator and mean
mats  <- compute_matrices(obj = example,
                          theta = root,
                          numDeriv_options = list(method = 'Richardson'))
Sigma <- compute_sigma(mats)
Sigma



## ----SB_example6, echo=TRUE----------------------------------------------



## ----SB_example7, echo=TRUE----------------------------------------------




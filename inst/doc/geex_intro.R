## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
library(geex)
library(dplyr)
library(inferference)
library(sandwich)
library(xtable)
library(moments)
library(MASS)
library(knitr)
library(rprojroot)
# child.path <- normalizePath(paste0(find_package_root_file(), '/vignettes/examples/'))
opts_knit$set(progress = TRUE, verbose = TRUE, child.path = 'examples/')
# library(microbenchmark)

## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
library(geex)
library(dplyr)
library(inferference)
library(sandwich)
library(xtable)
# library(microbenchmark)

## ----functions_results, echo = FALSE-------------------------------------
print_pmatrix <- function(object, digits = 4){
  if(!is.matrix(object)){
    object <- matrix(object, nrow = 1)
  }
  
  paste0('$', print(xtable(object, align=rep("",ncol(object)+1), digits =digits), comment = FALSE,
        floating=FALSE, tabular.environment="pmatrix", hline.after=NULL, 
        include.rownames=FALSE, include.colnames=FALSE, print.results = FALSE), '$')
}

first_diff_dec <- function(x){
  -floor(log10(abs(x)))
}

print_results <- function(results, label, caption){
  r <- results
  cat('\\begin{table}[H] \n',
      '\\centering \n',
      '\\label{', label, '} \n',
      '\\caption{"', caption, '"} \n',
      '\\begin{tabular}{lcc} \n',
      ' & $\\hat{\\theta}$ & $\\hat{\\Sigma}$  \\\\ \n',
      'Closed form &', print_pmatrix(r$cls$parameters),  '&', print_pmatrix(r$cls$vcov), '\\\\ \n',
      'geex &',  print_pmatrix(r$geex$parameters),  '&', print_pmatrix(r$geex$vcov), '\\\\ \n',
      'Decimal of difference &',  print_pmatrix(first_diff_dec(r$cls$parameters - r$geex$parameters), d = 0),  '&',
                                  print_pmatrix(first_diff_dec(r$cls$vcov - r$geex$vcov), d = 0), '\\\\ \n',
      '\\end{tabular} \n', 
      '\\end{table}')
}

## ----SB1_setup, echo=FALSE-----------------------------------------------
n  <- 100
mu <- 5
sigma <- 2
dt <- data.frame(Y = rnorm(n, mean = mu, sd = sigma), id = 1:n)

## ----SB1_eefun, echo=FALSE, results='hide'-------------------------------
SB1_eefun <- function(data){
  function(theta){
    with(data,
      c(Y - theta[1],
       (Y - theta[1])^2 - theta[2] )
    )
  }
}

## ----SB1_run, echo=TRUE--------------------------------------------------
estimates <- estimate_equations(eeFUN = SB1_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1,1))

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
theta_cls <- summarize(dt, p1 = mean(Y), p2 = var(Y) * (n() - 1)/ n() )

# closed form
Sigma_cls <- (solve(A) %*% B %*% t(solve(A))) / n

## ----SB1_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))

print_results(results, 'test', 'test')

## ----SB_setup, echo=FALSE------------------------------------------------
n  <- 100
muY <- 5
sigmaY <- 2
muX <- 2
sigmaX <- 0.2
dt <- data.frame(Y  = rnorm(n, mean = muY, sd = sigmaY), 
                 X  = rnorm(n, mean = muX, sd = sigmaX),
                 id = 1:n)

## ----SB2_eefun, echo = FALSE---------------------------------------------
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
estimates <- estimate_equations(eeFUN = SB2_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1, 1, 1))

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
theta_cls <- summarize(dt, p1 = mean(Y), p2 = mean(X), p3 = p1/p2)

## closed form covariance
Sigma_cls <- (solve(A) %*% B %*% t(solve(A))) / n

## ----SB2_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'test', 'test')

## ----SB3_setup, echo=FALSE-----------------------------------------------
n  <- 100
mu <- 5
sigma <- 4
set.seed(100) # running into issue where sqrt(theta2) and log(theta2) return NaN for some seeds
dt <- data.frame(Y  = rnorm(n, mean = mu, sd = sigma), 
                 id = 1:n)

## ----SB3_eefun, echo = FALSE---------------------------------------------
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
estimates <- estimate_equations(eeFUN= SB3_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1, 1, 1, 1))

## ----SB3_clsform, echo = FALSE-------------------------------------------
## closed form roots
theta_cls <- summarize(dt, p1 = mean(Y), p2 = sum((Y - p1)^2)/n(), p3 = sqrt(p2), p4 = log(p2))

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

## ----SB3_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'ex3_results', 'Example 3')

## ----SB4_setup, echo=FALSE-----------------------------------------------
n  <- 100

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

## ----SB4_eefun, echo = FALSE---------------------------------------------
SB4_eefun <- function(data){
  function(theta){
    with(data,
      c(theta[1] - T_,
        theta[2] - W,
        (Y - (theta[3] * W)) * (theta[2] - W),
        (Y - (theta[4] * W)) * (theta[1] - T_))
    )
  }
}

## ----SB4_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(eeFUN = SB4_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1, 1, 1, 1))

## ----SB4_clsform, echo = TRUE--------------------------------------------
YW_model <- lm(Y ~ W, data = dt)
YT_model <- lm(Y ~ T_, data = dt)
WT_model <- lm(W ~ T_, data = dt)
## closed form roots
theta_cls <- c(theta1 = mean(dt$T_),
  theta2 = mean(dt$W),
  theta3 = coef(YW_model)[2],
  theta4 = coef(YT_model)[2]/coef(WT_model)[2])

## closed form covariance
# Not sure how compute SB's closed form since it depends on X, which is
# supposed to be unobserved.
Sigma_cls <- matrix(NA, nrow = 2, ncol = 2)

## ----SB4_results, echo = FALSE, results = 'asis'-------------------------
# primary interest lies in the lower 2 x 2 submatrix of the asymptotic variance matrix
estimates$vcov <- estimates$vcov[3:4, 3:4]
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'Example 4', 'Example 4')

## ----SB5_setup, echo=FALSE-----------------------------------------------
n <- 100
theta0 <- 0
theta_tru <- 2
sigma <- 1
dt <- data.frame(X = rnorm(n, mean = 2, sd = sigma),
                 id = 1:n)

## ----SB5_eefun, echo = TRUE----------------------------------------------
F0 <- function(y, theta0, distrFUN = pnorm){
  distrFUN(y - theta0, mean = 0)
}

f0 <- function(y, densFUN){
  densFUN(y, mean = 0)
}

integrand <- function(y, densFUN = dnorm){
  f0(y, densFUN = densFUN)^2
}

IC_denom <- integrate(integrand, lower = -Inf, upper = Inf)$value

SB5_eefun <- function(data, theta0 = 0){
  Xi <- data$X
  IC_HL <- (1/IC_denom) * (F0(Xi, theta0) - 0.5)
  function(theta){
     c(IC_HL - (theta[1] - theta0),
       Xi - theta[2]) 
  }
}

## ----SB5_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(eeFUN = SB5_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1, 1))

## ----SB5_clsform, echo = TRUE--------------------------------------------
X <- dt$X
pair_means <- numeric(length(X) - 1)
for(i in 1:(length(X) - 1)){
 pair_means[i] <-  (X[i] + X[i + 1])/2
}

theta_cls <- c(median(pair_means), mean(X))

## closed form covariance
# Not sure how compute SB's closed form since it depends on X, which is
# supposed to be unobserved.
Sigma_cls <- matrix(c(1/(12 * IC_denom^2) / n, NA, NA, NA), 
                    nrow = 2, ncol = 2, byrow = TRUE)

## ----SB5_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'Example 5', 'Example 5')

## ----SB6_setup, echo=FALSE-----------------------------------------------
n <- 100
theta_tru <- 2
sigma <- 1
dt <- data.frame(Y = rnorm(n, mean = 2, sd = sigma),
                 id = 1:n)

## ----SB6_eefun, echo = FALSE---------------------------------------------
SB6_eefun <- function(data, k = 1.5){
  function(theta){
    x <- data$Y - theta[1]
    if(abs(x) <= k) x else sign(x) * k
  }
}

## ----SB6_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(eeFUN = SB6_eefun, 
                                data  = dt, units = 'id', 
                                roots = 1)

## ----SB6_clsform, echo = TRUE--------------------------------------------
theta_cls <- MASS::huber(dt$Y)$mu

psi_k <- function(x, k = 1.5){
  if(abs(x) <= k) x else sign(x) * k
}

A <- lapply(dt$Y, function(y){
  x <- y - theta_cls
  -numDeriv::grad(psi_k, x = x)
}) %>% unlist() %>% mean()

B <- lapply(dt$Y, function(y){
  x <- y - theta_cls
  psi_k(x = x)^2
}) %>% unlist() %>% mean()

## closed form covariance
Sigma_cls <- matrix(A * B * A / n)

## ----SB6_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'Example 6', 'Example 6')

## ----SB7_setup, echo=FALSE-----------------------------------------------
n <- 100
theta_tru <- 2
sigma <- 1
dt <- data.frame(Y = rnorm(n, mean = 2, sd = sigma),
                 id = 1:n)

## ----SB7_eefun, echo = FALSE---------------------------------------------
SB7_eefun <- function(data){
  function(theta){
    with(data,
      c(0.5  - (Y <= theta[1]),
        0.65 - (Y <= theta[2]))
    )
  }
}

## ----SB7_run, echo = TRUE, eval=FALSE------------------------------------
#  estimates <- estimate_equations(eeFUN = SB7_eefun,
#                                  data  = dt, units = 'id',
#                                  roots = c(.5, .65))

## ----SB7_clsform, echo = TRUE--------------------------------------------
theta_cls <- c(quantile(dt$Y, 0.5), quantile(dt$Y, 0.65))


## ----SB7_results, echo = FALSE, results = 'asis', eval=FALSE-------------
#  results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
#  print_results(results, 'Example 7', 'Example 7')

## ----SB8_setup, echo=FALSE-----------------------------------------------
n <- 50
beta <- c(0.5, 2)
dt <- data_frame(X  = rep(0:1, each = n/2),
                 e  = rnorm(n),
                 Y  = as.numeric(cbind(1, X) %*% beta) + e,
                 id = 1:n)

## ----SB8_eefun, echo = FALSE---------------------------------------------
psi_k <- function(x, k = 1.345){
  if(abs(x) <= k) x else sign(x) * k
}

SB8_eefun <- function(data){
  Yi <- data$Y
  xi <- model.matrix(Y ~ X, data = data)
  function(theta){
    r <- Yi - xi %*% theta
    c(psi_k(r) %*% xi)
  }
}

## ----SB8_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(eeFUN = SB8_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1, 1))

## ----SB8_clsform, echo = TRUE--------------------------------------------
m <- MASS::rlm(Y ~ X, data = dt, method = 'M')
theta_cls <- coef(m)
Sigma_cls <- vcov(m)


## ----SB8_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'Example 8', 'Example 8')

## ----SB9_setup, echo=FALSE-----------------------------------------------
n <- 100
beta <- c(0.5, 2, .1)
dt <- data_frame(X1 = rep(0:1, each = n/2), 
                 X2 = rep(0:1, times = n/2),
                 Y  = rbinom(n, 1, prob = as.numeric(plogis(cbind(1, X1, X2) %*% beta))),
                 id = 1:n)

## ----SB9_eefun, echo = FALSE---------------------------------------------
SB9_eefun <- function(data){
  Yi <- data$Y
  xi <- model.matrix(Y ~ X1 + X2, data = data, drop = FALSE)
  function(theta){
    lp <- xi %*% theta
    mu <- plogis(lp)
    D  <- t(xi) %*% dlogis(lp)
    V  <- mu * (1 - mu)
    D %*% solve(V) %*% (Yi - mu)
  }
}

## ----SB9_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(eeFUN = SB9_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(1, 1, 1))

## ----SB9_clsform, echo = TRUE--------------------------------------------
m <- glm(Y ~ X1 + X2, data = dt, family = binomial(link = 'logit'))
theta_cls <- coef(m)
Sigma_cls <- sandwich(m)


## ----SB9_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'Example 9', 'Example 9')

## ----SB10_setup, echo=FALSE----------------------------------------------
shaq <- data_frame(game = 1:23,
                 ft_made = c(4, 5, 5, 5, 2, 7, 6, 9, 4, 1, 13, 5, 6, 9, 7, 3, 8, 1, 18, 3, 10, 1, 3),
                 ft_attp = c(5, 11, 14, 12, 7, 10, 14, 15, 12, 4, 27, 17, 12, 9, 12, 10, 12, 6, 39, 13, 17, 6, 12))


## ----SB10_eefun, echo = FALSE--------------------------------------------
SB10_eefun <- function(data){
  Y <- data$ft_made
  n <- data$ft_attp
  function(theta){
    p <- theta[2]
    c(((Y - (n * p))^2)/(n * p * (1 - p))  - theta[1], 
      Y - n * p)
  }
}

## ----SB10_run, echo = TRUE-----------------------------------------------
estimates <- estimate_equations(eeFUN = SB10_eefun, 
                                data  = shaq, units = 'game', 
                                roots = c(.5, .5))

## ----SB10_clsform, echo = TRUE-------------------------------------------
V11 <- function(p) {
  k <- length(nrow(shaq))
  sumn <- sum(shaq$ft_attp)
  sumn_inv <- sum(1/shaq$ft_attp)
  term2_n <- 1 - (6 * p) + (6 * p^2)
  term2_d <- p * (1 - p) 
  term2 <- term2_n/term2_d
  print(term2)
  term3 <- ((1 - 2 * p)^2)/( (sumn/k) * p * (1 - p))
  print(term3)
  2 + (term2 * (1/k) * sumn_inv)  - term3
}

### ???? I keep getting a negative value for V11

p_tilde <- sum(shaq$ft_made)/sum(shaq$ft_attp)
V <- V11(.45)
V
pnorm(estimates$parameters[1], mean = 1, sd = sqrt(V))


## ----SB10_results, echo = FALSE, results = 'asis', eval = FALSE----------
#  results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
#  print_results(results, 'Example 10', 'Example 10')


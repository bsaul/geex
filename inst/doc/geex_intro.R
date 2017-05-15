## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
library(geex)
library(dplyr)
library(inferference)
library(sandwich)
library(xtable)
library(moments)
library(MASS)
library(knitr)
# library(rprojroot)
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
  
  paste0('$', print(xtable::xtable(object, align=rep("",ncol(object)+1), digits =digits), comment = FALSE,
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
                                data  = dt,
                                units = 'id', 
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

print_results(results, 'ex1', 'Comparing estimates from closed form versus geex')

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
# X <- rnorm(n, sd = 1)

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
# sigma2_e <- var(residuals(YW_model))
# sigma2_W <- var(dt$W)
# sigma2_U <- var(residuals(YT_model))

# ((sigma_e^2 * (1 + sigma_U^2 * sigma_e^2) + beta^2 * (sigma_U^2 * 1))/(1 + sigma_U^2 * sigma_e^2)^2) / (100^2)
# ((sigma_e^2 * (1 + sigma_U^2 * sigma_e^2) + beta^2 * (sigma_U^2 * 1))/(1 + sigma_U^2 * sigma_e^2)^2) / (100^2)
Sigma_cls <- matrix(NA, nrow = 2, ncol = 2)

## ----SB4_results, echo = FALSE, results = 'asis'-------------------------
# primary interest lies in the lower 2 x 2 submatrix of the asymptotic variance matrix
estimates$vcov <- estimates$vcov[3:4, 3:4]
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))

print_results(results, 'Example 4', 'Example 4')

## ----SB5_setup, echo=FALSE-----------------------------------------------
library(ICSNP)
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

SB5_eefun <- function(data){
  Xi <- data$X
  function(theta){
     IC_HL <- (1/IC_denom) * (F0(Xi, theta[1]) - 0.5)
     c(IC_HL,
       Xi - theta[2]) 
  }
}

## ----SB5_run, echo = TRUE------------------------------------------------
estimates <- estimate_equations(eeFUN = SB5_eefun, 
                                data  = dt, units = 'id', 
                                roots = c(2, 1))

## ----SB5_clsform, echo = TRUE--------------------------------------------
theta_cls <- c(hl.loc(dt$X), mean(dt$X))

## closed form covariance
# Not sure how compute SB's closed form since it depends on X, which is
# supposed to be unobserved.
Sigma_cls <- matrix(c(1/(12 * IC_denom^2) / n, NA, NA, var(dt$X)/100), 
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
theta_cls <- MASS::huber(dt$Y, tol = 1e-10)$mu

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
Sigma_cls <- matrix(1/A * B * 1/A / n)

## ----SB6_results, echo = FALSE, results = 'asis'-------------------------
results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
print_results(results, 'Example 6', 'Example 6')

## ----SB7_setup, echo=FALSE-----------------------------------------------
library(geex)
set.seed(9)
n <- 100
theta_tru <- 2
sigma <- 1
dt <- data.frame(Y = rnorm(n, mean = theta_tru, sd = sigma), id = 1:n)

## ----SB7_eefun, echo = TRUE----------------------------------------------
SB7_eefun <- function(data){
  function(theta){
    0.5  - (data$Y <= theta[1])
  }
}

## ----check_array, echo = FALSE, eval = TRUE------------------------------
# this is not an exported function from geex, so add it here
check_array <- function(object){
  if(is.array(object)){
    object
  } else if(is.numeric(object)){
    array(object, dim = c(1, 1, length(object)))
  } else if(is.matrix(object)){
    array(object, dim = c(1, 1, length(object)))
  } else {
    stop('Object is not an array, matrix, or numeric')
  }
}

## ----eeroot_modified, echo = TRUE, eval=TRUE-----------------------------
eeroot_mod <- function(geex_list,
                   start         = NULL,
                   rootsolver    = rootSolve::multiroot,
                   root_options  = NULL,
                   apprx_fun     = NULL,
                   apprx_options = NULL,
                   ...){

  # Create estimating equation functions per group
  psi_i <- lapply(geex_list$splitdt, function(data_i){
    geex_list$eeFUN(data = data_i, ...)
  })

  # Create psi function that sums over all ee funs
  psi <- function(theta){
    psii <- lapply(psi_i, function(f) {
      do.call(f, args = append(list(theta = theta), geex_list$ee_args))
    })
    apply(check_array(simplify2array(psii)), 1, sum)
  }

  # apprx_fun is a function that manipulates the psi function and returns a new psi function
  if(!is.null(apprx_fun)){
    psi <- do.call(apprx_fun, args = append(list(psi = psi), apprx_options))
  }
  # Find roots of psi
  rargs <- append(root_options, list(f = psi, start = start))
  do.call(rootsolver, args = rargs)
}

## ----splinefun, echo = TRUE----------------------------------------------
spline_approx <- function(psi, eval_theta){
  ### Use splinefun ####
  psi2 <- Vectorize(psi)
  psi <- function(theta) splinefun(x = eval_theta, psi2(eval_theta))(theta)
  psi
}


myList <- list(eeFUN = SB7_eefun, splitdt = split(dt, f = dt$id))

root_spline1 <- eeroot_mod(
  geex_list = myList,
  start     = 2,
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(1, 5, by = .5))
)

# Compare to the truth
median(dt$Y) - root_spline1$root

# But notice that the basis of the spline matters too
root_spline2 <- eeroot_mod(
  geex_list = myList,
  start     = 2,
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(1, 5, by = .1))
)

# Compare to the truth
median(dt$Y) - root_spline2$root


# But notice that the basis of the spline matters too
root_spline3 <- eeroot_mod(
  geex_list = myList,
  start     = 2,
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(-5, 5, by = 1))
)

# Compare to the truth
median(dt$Y) - root_spline3$root


## ----gam_approx, echo = TRUE---------------------------------------------
library(mgcv)
gam_approx <- function(psi, eval_theta){
  ### Use splinefun ####
  psi2 <- Vectorize(psi)
  Y    <- psi2(eval_theta)
  gam_basis <- gam(Y ~ s(eval_theta))
  psi <- function(theta) predict(gam_basis, newdata = data.frame(eval_theta = theta))
  psi
}

root_gam1 <- eeroot_mod(
  geex_list = myList,
  start     = 2,
  apprx_fun = gam_approx,
  apprx_options = list(eval_theta = seq(1, 5, by = .3))
)

# Compare to the truth
median(dt$Y) - root_gam1$root


# Still, the basis of the spline influences the result
root_gam2 <- eeroot_mod(
  geex_list = myList,
  start     = 2,
  apprx_fun = gam_approx,
  apprx_options = list(eval_theta = seq(1.5, 3.5, length.out = 20))
)

# Compare to the truth
median(dt$Y) - root_gam2$root


## ----gam_approx_10000, echo = TRUE---------------------------------------
dt2 <- data.frame(Y = rnorm(n, mean = theta_tru, sd = sigma), id = 1:10000)
myList2 <- list(eeFUN = SB7_eefun, splitdt = split(dt2, f = dt2$id))

root_spline4 <- eeroot_mod(
  geex_list = myList2,
  start     = 2,
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(0, 5, by = .1))
)

median(dt2$Y) - root_spline4$root

root_gam3 <- eeroot_mod(
  geex_list = myList2,
  start     = 2,
  apprx_fun = gam_approx,
  apprx_options = list(eval_theta = seq(0, 5, by = .1))
)

# Compare to the truth
median(dt$Y) - root_gam3$root

## ----compute_matrices_mod, echo = TRUE-----------------------------------
compute_matrices_mod <- function(geex_list,
                             theta,
                             numDeriv_options = list(method = 'Richardson'),
                             silent = TRUE,
                             apprx_fun     = NULL,
                             apprx_options = NULL,
                             ...){
  if(is.null(geex_list$ee_args)){
    ee_args <- NULL
  }

  with(geex_list, {

    # Create list of estimating eqn functions per unit
    psi_i <- lapply(splitdt, function(data_i){
      f <- eeFUN(data = data_i, ...)
      if(!is.null(apprx_fun)){
        f <- do.call(apprx_fun, args = append(list(psi = f), apprx_options))
      }
      f
    })

    # Compute the negative of the derivative matrix of estimating eqn functions
    # (the information matrix)
    A_i <- lapply(psi_i, function(ee){
      args <- append(list(fun = ee, x = theta), numDeriv_options)
      val  <- do.call(numDeriv::jacobian, args = append(args, ee_args))
      -val
    })
    A_i_array <- check_array(simplify2array(A_i))
    A   <- apply(A_i_array, 1:2, sum)

    # Compute outer product of observed estimating eqns
    B_i <- lapply(psi_i, function(ee) {
      ee_val <- do.call(ee, args = append(list(theta = theta), ee_args))
      ee_val %*% t(ee_val)
    })
    B   <- apply(check_array(simplify2array(B_i)), 1:2, sum)

    list(A = A, A_i = A_i, B = B, B_i = B_i)
  })
}

## ----matrices, echo = TRUE-----------------------------------------------
spline_matrices1 <- compute_matrices_mod(
  geex_list = myList,
  theta     = median(dt$Y),
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(0, 4, by = .5))
)

compute_sigma(A = spline_matrices1$A, B = spline_matrices1$B) 


gam_matrices1 <- compute_matrices_mod(
  geex_list = myList,
  theta     = median(dt$Y),
  apprx_fun = gam_approx,
  apprx_options = list(eval_theta = seq(0, 4, by = .3))
)

compute_sigma(A = gam_matrices1$A, B = gam_matrices1$B) 

# V(theta) from Stefanski and Boos
((0.5 * (1 -0.5))/((dnorm(median(dt2$Y), mean = mean(dt2$Y), sd = sd(dt2$Y)))^2) ) 

## ----dens_approx, echo = TRUE--------------------------------------------
dens <- density(dt$Y)
ff <- splinefun(x = dens$x, y = dens$y)
density_apprx <- function(psi){
  function(theta) ff(theta)
}

density_matrices1 <- compute_matrices_mod(
  geex_list = myList,
  theta     = median(dt$Y),
  apprx_fun = density_apprx
)

compute_sigma(A = density_matrices1$A, B = density_matrices1$B) 


## ----eeroot_modified2, echo = TRUE, eval=TRUE----------------------------

# Put apprx_fun code in same place as compute_matrices
eeroot_mod2<- function(geex_list,
                   start         = NULL,
                   rootsolver    = rootSolve::multiroot,
                   root_options  = NULL,
                   apprx_fun     = NULL,
                   apprx_options = NULL,
                   ...){

  # Create estimating equation functions per group
  psi_i <- lapply(geex_list$splitdt, function(data_i){
      f <- geex_list$eeFUN(data = data_i, ...)
      if(!is.null(apprx_fun)){
        f <- do.call(apprx_fun, args = append(list(psi = f), apprx_options))
      }
      f
  })

  # Create psi function that sums over all ee funs
  psi <- function(theta){
    psii <- lapply(psi_i, function(f) {
      do.call(f, args = append(list(theta = theta), geex_list$ee_args))
    })
    apply(check_array(simplify2array(psii)), 1, sum)
  }
  # Find roots of psi
  rargs <- append(root_options, list(f = psi, start = start))
  do.call(rootsolver, args = rargs)
}

## ---- echo = TRUE--------------------------------------------------------
# But notice that the basis of the spline matters too
root_spline5 <- eeroot_mod2(
  geex_list = myList,
  start     = 2,
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(-5, 5, by = 1))
)

# Compare to the truth
median(dt$Y) - root_spline5$root

root_gam5 <- eeroot_mod2(
  geex_list = myList,
  start     = 2,
  apprx_fun = gam_approx,
  apprx_options = list(eval_theta = seq(0, 5, by = .1))
)

# Compare to the truth
median(dt$Y) - root_gam5$root

## ---- echo = TRUE--------------------------------------------------------
# But notice that the basis of the spline matters too
root_spline6 <- eeroot_mod2(
  geex_list = myList2,
  start     = 2,
  apprx_fun = spline_approx,
  apprx_options = list(eval_theta = seq(-5, 5, by = 1))
)

# Compare to the truth
median(dt$Y) - root_spline6$root

root_gam6 <- eeroot_mod(
  geex_list = myList2,
  start     = 2,
  apprx_fun = gam_approx,
  apprx_options = list(eval_theta = seq(0, 5, by = .1))
)

# Compare to the truth
median(dt$Y) - root_gam6$root

## ----SB7_eefun2, echo = TRUE---------------------------------------------
SB7_eefun2 <- function(data){
  function(theta){
    with(data,
    c(0.25  - (Y <= theta[1]),
      0.5   - (Y <= theta[2]),
      0.75  - (Y <= theta[3]))
    )
  }
}

## ----multi_psi, echo = TRUE, eval=FALSE----------------------------------
#  gam_approx_multi <- function(psi, eval_theta1, eval_theta2, eval_theta3){
#    ### Use splinefun ####
#    # psi2 <- Vectorize(psi)
#    Y <- matrix(NA, nrow = length(eval_theta1), ncol = 3)
#    for(i in 1:length(eval_theta1)){
#        Y[i, ]  <- psi(c(eval_theta1[i], eval_theta2[i], eval_theta3[i]))
#    }
#  
#    gam_data <- data.frame(
#      Y1 = Y[, 1],
#      Y2 = Y[, 2],
#      Y3 = Y[, 3],
#      x1 = eval_theta1,
#      x2 = eval_theta2,
#      x3 = eval_theta3
#    )
#  
#    gam_basis <- gam(list(Y1 ~ s(x1),
#                          Y2 ~ s(x2),
#                          Y3 ~ s(x3)),
#                     family = mvn(d = 3),
#                     data = gam_data)
#    print(gam_basis)
#    psi <- function(theta) {
#      predict(gam_basis, newdata = data.frame(eval_theta1 = theta[1],
#                                              eval_theta2 = theta[2],
#                                              eval_theta3 = theta[3]))
#    }
#    psi
#  }
#  
#  SB7_eefun2(dt[1, ])(c(1, 2, 3))
#  
#  myList3 <- list(eeFUN = SB7_eefun2, splitdt = split(dt, f = dt$id))
#  
#  root_gam1 <- eeroot_mod(
#    geex_list = myList3,
#    start     = c(1, 2, 3),
#    apprx_fun = gam_approx_multi,
#    apprx_options = list(eval_theta1 = seq(-5, 5, by = .1),
#                         eval_theta2 = seq(-5, 5, by = .1),
#                         eval_theta3 = seq(-5, 5, by = .1))
#  )
#  
#  quantile(dt$Y, probs = c(.25, .5, .75))
#  
#  

## ----SB7_run, echo = TRUE, eval=FALSE------------------------------------
#  estimates <- estimate_equations(eeFUN = SB7_eefun,
#                                  data  = dt, units = 'id',
#                                  findroots = FALSE,
#                                  roots = c(2))

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
                                data  = shaq, 
                                units = 'game', 
                                numDeriv_options = list(method.args = list(eps = 1e-7, r = 10, zero.tol = .Machine$double.eps)), 
                                roots = c(.5, .5))

## ----SB10_clsform, echo = TRUE-------------------------------------------
V11 <- function(p) {
  k    <- nrow(shaq)
  sumn <- sum(shaq$ft_attp)
  sumn_inv <- sum(1/shaq$ft_attp)
  term2_n  <- 1 - (6 * p) + (6 * p^2)
  term2_d <- p * (1 - p) 
  term2  <- term2_n/term2_d
  term3  <- ((1 - (2 * p))^2) / ((sumn/k) * p * (1 - p))
  2 + (term2 * (1/k) * sumn_inv) - term3
}

p_tilde <- sum(shaq$ft_made)/sum(shaq$ft_attp)
V11_hat <- V11(p_tilde)/23

# Compare variance estimates
V11_hat
estimates$vcov[1, 1]

# Compare p-values
pnorm(35.51/23, mean  = 1, sd = sqrt(V11_hat), lower.tail = FALSE)

pnorm(estimates$parameters[1], 
      mean = 1, 
      sd = sqrt(estimates$vcov[1, 1]),
      lower.tail = FALSE)


## ----SB10_results, echo = FALSE, results = 'asis', eval = FALSE----------
#  results <- list(geex = estimates, cls = list(parameters = theta_cls, vcov = Sigma_cls))
#  print_results(results, 'Example 10', 'Example 10')


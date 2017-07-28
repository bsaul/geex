## ---- echo = TRUE, eval=FALSE--------------------------------------------
#  correct_by_nothing <- function(A, A_i, B, B_i){
#    compute_sigma(A = A, B = B)
#  }

## ---- echo = TRUE--------------------------------------------------------
correct_by_bias <- function(A, A_i, B, B_i, b){
  Ainv <- solve(A)

  H_i <- lapply(A_i, function(m){
    diag( (1 - pmin(b, diag(m %*% Ainv) ) )^(-0.5) )
  })

  Bbc_i <- lapply(seq_along(B_i), function(i){
    H_i[[i]] %*% B_i[[i]] %*% H_i[[i]]
  })
  Bbc   <- apply(simplify2array(Bbc_i), 1:2, sum)

  compute_sigma(A = A, B = Bbc)
}

## ----FAY1_eefun, echo = TRUE---------------------------------------------
gee_eefun <- function(data, formula, family){
  X <- model.matrix(object = formula, data = data)
  Y <- model.response(model.frame(formula = formula, data = data))
  n <- nrow(X)
  function(theta, alpha, psi){
    mu  <- family$linkinv(X %*% theta)
    Dt  <- t(X) %*% diag(as.numeric(mu), nrow = n)
    A   <- diag(as.numeric(family$variance(mu)), nrow = n)
    R   <- matrix(alpha, nrow = n, ncol = n)
    diag(R) <- 1
    V   <- psi * (sqrt(A) %*% R %*% sqrt(A))
    Dt %*% solve(V) %*% (Y - mu)
  }
}

## ----setup_gee, echo = TRUE, message = FALSE, results = 'hide'-----------
g <- gee::gee(breaks~tension, id=wool, data=warpbreaks, corstr="exchangeable")
guo <- saws::geeUOmega(g)

## ----correction_run, echo = TRUE-----------------------------------------
library(geex)
results <- estimate_equations(
  eeFUN = gee_eefun, data  = warpbreaks, 
  units = 'wool', roots = coef(g), compute_roots = FALSE,
  outer_eeargs = list(formula = breaks ~ tension, 
                      family  = gaussian()),
  inner_eeargs = list(alpha   = g$working.correlation[1,2], 
                      psi     = g$scale), 
  corrections_list = list(
   bias_correction_.1 = list(correctFUN = correct_by_bias, 
                             correctFUN_control = list(b = .1)),
   bias_correction_.3 = list(correctFUN = correct_by_bias, 
                             correctFUN_control = list(b = .3)))) 

## ----correction_comparison, echo = FALSE, results = 'hide'---------------
saws::saws(guo, method = 'd1')$V 
results$vcov

saws::saws(guo, method = 'd4', bound  = 0.1)$V
results$corrections$bias_correction_.1

saws::saws(guo, method = 'd4', bound  = 0.3)$V
results$corrections$bias_correction_.3


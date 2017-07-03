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
library(geex)
library(inferference)
mglm    <- glm(A ~ X1, data = vaccinesim, family = binomial)
estimates <- estimate_equations(eeFUN = eefun,
                   data = vaccinesim,
                   roots = c(-.35, 0),
                   outer_eeargs = list(model = mglm))

# Compare point estimates
estimates$parameters # from GEEX
coef(mglm) # from the GLM function

# Compare variance estimates
estimates$vcov
sandwich::sandwich(mglm)


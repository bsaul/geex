## ----eefun, echo=TRUE----------------------------------------------------
eefun <- function(data, model, alpha){
  X <- model.matrix(model, data = data)
  A <- model.response(model.frame(model, data = data))
  
  function(theta){
    p  <- length(theta)
    p1 <- length(coef(model))
    lp  <- X %*% theta[1:p1]
    rho <- plogis(lp)

    hh  <- ((rho/alpha)^A * ((1-rho)/(1-alpha))^(1 - A)) 
    IPW <- 1/(exp(sum(log(hh))))

    score_eqns <- apply(X, 2, function(x) sum((A - rho) * x))
    with(data, {
      ce0 <- mean(Y * (A == 0)) * IPW / (1 - alpha)
      ce1 <- mean(Y * (A == 1)) * IPW / (alpha)
      
      c(score_eqns,
        ce0 - theta[p - 1],
        ce1 - theta[p])
    })
  }
}

## ----example2, echo =TRUE------------------------------------------------
library(geex)
library(inferference)

test <- interference(
  formula = Y | A ~ X1 | group, 
  data   = vaccinesim,
  model_method = 'glm',
  allocations = c(.35, .4))

mglm <- glm(A ~ X1, data = vaccinesim, family = binomial)

ce_estimates <- m_estimate(
  estFUN = eefun,
  data  = vaccinesim,
  units = 'group',
  root_control = setup_root_solver(start = c(coef(mglm), .4,  .13)),
  outer_args = list(alpha = .35, model = mglm)
)

roots(ce_estimates)

# Compare parameter estimates
direct_effect(test, allocation = .35)$estimate
roots(ce_estimates)[3] - roots(ce_estimates)[4]

# conpare SE estimates
L <- c(0, 0, 1, -1)
Sigma <- vcov(ce_estimates)
sqrt(t(L) %*% Sigma %*% L)  # from GEEX
direct_effect(test, allocation = .35)$std.error # from inferference


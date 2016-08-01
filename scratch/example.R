library(inferference)


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

mglm    <- glm(A ~ X1, data = vaccinesim, family = binomial)
s_data  <- split(vaccinesim, vaccinesim$group)
example <- list(eeFUN = eefun, splitdt = s_data, theta = coef(mglm))

mats <- compute_matrices(obj = example,
                         contrast = c(0, 1),
                         corrections = c('df'),
                         correction_options = list(b = 0.75),
                         model = mglm)

compute_sigma(mats)


## TESTING

psi <- lapply(s_data, function(data_i){
  psi_i <- eefun(data = data_i, model = mglm)
})

psi[[1]]
psi1 <- function(theta){
  psi_ <- lapply(psi, function(f) f(theta))
  abs(sum(unlist(psi_)))
}

psi_gr <- function(theta){
  numDeriv::jacobian(psi1, x = theta)
}

optimx::optimx(par = c(0.01, 0.01), fn = psi1, gr = psi_gr)

optimx::optimx(par = coef(mglm), fn = psi1, gr = psi_gr)

coef(mglm)

psi1(coef(mglm))




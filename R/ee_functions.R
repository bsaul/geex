#------------------------------------------------------------------------------#
#' Make Estimating Equation functions
#'
#' Converts a model object into
#'
#' @param model a model object object
#' @param data data with which to create the estimating equation function
#' @param ... passed to methods
#' @export
#------------------------------------------------------------------------------#
make_eefun <- function(model, data, ...)
{
  UseMethod("make_eefun")
}

#------------------------------------------------------------------------------#
#' geeglm Estimating Equations
#'
#' Create estimating equation function from a \code{geeglm} object
#'
#' @inheritParams make_eefun
#' @export
#------------------------------------------------------------------------------#

make_eefun.geeglm <- function(model, data)
{
  if(model$corstr != 'independence'){
    stop("only independence working correlation is supported at this time")
  }

  X <- model.matrix(model$formula, data = data)
  # DO NOT use model.matrix(geepack_obj, data = subdata)) -
  # returns entire model matrix, not just the subset
  Y  <- model.response(model.frame(model, data = data))
  n  <- length(Y)
  p  <- length(coef(model))
  phi    <- summary(model)$dispersion$Estimate
  family <- model$family$family
  link   <- model$family$link
  invlnk <- model$family$linkinv
  family_link <- paste(family, link, sep = '_')

  function(theta){
    lp <- X %*% theta # linear predictor
    f  <- as.numeric(invlnk(lp))  # fitted values
    r  <- Y - f       # residuals

    ### TODO: this is cludgy and needs to be reworked to be more general
    ### TODO: how to handle weights
    if(family_link == 'gaussian_identity'){
      D <- X
      V <- phi * diag(1, nrow = n, ncol = n)
    } else if(family_link == 'binomial_logit'){
      D <- apply(X, 2, function(x) x * exp(lp)/((1+exp(lp))^2) )
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }
    t(D) %*% solve(V) %*% (r)
  }
}

#------------------------------------------------------------------------------#
#' glmer Estimating Equations
#'
#' Create estimating equation function from a \code{merMod} object
#'
#' @inheritParams make_eefun
#' @export
#------------------------------------------------------------------------------#

make_eefun.merMod <- function(object, data, numderiv_opts = NULL)
{
  ## Warnings ##
  if(length(lme4::getME(object, 'theta')) > 1){
    stop('make_eefun.merMod currently does not handle >1 random effect')
  }

  fm     <- get_fixed_formula(object)
  X      <- get_design_matrix(fm, data)
  Y      <- get_response(formula(object), data = data)
  family <- object@resp$family
  lnkinv <- family$linkinv
  objfun <- objFun_merMod(family$family)

  function(theta){
    args <- list(func = objfun, x = theta, response = Y, xmatrix = X, linkinv = lnkinv)
    do.call(numDeriv::grad, args = append(args, numderiv_opts))
  }
}

#------------------------------------------------------------------------------#
#' glmer Objective Fundtion
#'
#'@param family distribution family of objective function
#'@param ... additional arguments pass to objective function
#'@export
#------------------------------------------------------------------------------#

objFun_merMod <- function(family, ...){
  switch(family,
         binomial = objFun_glmerMod_binomial,
         stop('Objective function not defined'))
}

#------------------------------------------------------------------------------#
#' glmer Objective Function for Logistic-Normal Likelihood
#'
#' @param parms vector of parameters
#' @param response vector of response values
#' @param xmatrix the matrix of covariates
#' @param linkinv inverse link function
#' @export
#------------------------------------------------------------------------------#

objFun_glmerMod_binomial <- function(parms, response, xmatrix, linkinv)
{
  log(integrate(binomial_integrand, lower = -Inf, upper = Inf,
                parms    = parms,
                response = response,
                xmatrix  = xmatrix,
                linkinv  = linkinv)$value)

}

#------------------------------------------------------------------------------#
#' glmer Objective Function for Logistic-Normal Likelihood
#'
#' @inheritParams objFun_glmerMod_binomial
#' @export
#------------------------------------------------------------------------------#

binomial_integrand <- function(b, response, xmatrix, parms, linkinv){
  if(class(xmatrix) != 'matrix'){
    xmatrix <- as.matrix(xmatrix)
  }
  pr  <- linkinv( drop(outer(xmatrix %*% parms[-length(parms)], b, '+') ) )
  hh  <- dbinom(response, 1, prob = pr)
  hha <- apply(hh, 2, prod)
  hha * dnorm(b, mean = 0, sd = parms[length(parms)])
}

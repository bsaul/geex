#------------------------------------------------------------------------------#
# grab_eeFUN description:
# a generic function that takes a model object and "grabs" the estimating
# functions from the object.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Grab estimating functions from a model object
#'
#' @export
#------------------------------------------------------------------------------#

grab_estFUN <- function(object, ...){
  # S3 generic, for S3 dispatch
  UseMethod("grab_estFUN")
}

#------------------------------------------------------------------------------#
#' grab_estFUN generic
#'
#' Grabs the \code{.estFUN} from an \code{\linkS4class{m_estimation_basis}} object
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_estFUN")
setMethod("grab_estFUN", "estimating_function", function(object) object@.estFUN)

#------------------------------------------------------------------------------#
#' Grab estimating functions from a glm object
#'
#' Create estimating equation function from a \code{glm} object
#'
#' @inheritParams grab_estFUN
#' @param weights a scalar or vector of weight values
#' @export
#------------------------------------------------------------------------------#

grab_estFUN.glm <- function(object, data, weights = 1, ...){

  X  <- stats::model.matrix(object$formula, data = data)
  Y  <- as.numeric(stats::model.frame(grab_response_formula(object), data = data)[[1]])
  n  <- length(Y)
  p  <- length(stats::coef(object))
  phi    <- as.numeric(summary(object)$dispersion[1])
  W      <- weights
  family <- object$family$family
  link   <- object$family$link
  invlnk <- object$family$linkinv
  family_link <- paste(family, link, sep = '_')

  stopifnot(length(W) == 1 | length(W) == n)
  if(length(W) == 1){
    W <- rep(W, n)
  }

  function(theta){
    lp <- X %*% theta # linear predictor
    f  <- as.numeric(invlnk(lp))  # fitted values
    r  <- Y - f       # residuals

    ### TODO: this is cludgy and needs to be reworked to be more general
    if(family_link == 'gaussian_identity'){
      D <- X
      V <- phi * diag(1, nrow = n, ncol = n)
    } else if(family_link == 'binomial_logit'){
      D <- apply(X, 2, function(x) x * exp(lp)/((1+exp(lp))^2) )
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }

    t(D) %*% solve(V) %*% diag(W, nrow = n, ncol = n) %*% (r)
  }
}

#------------------------------------------------------------------------------#
#' Grab estimating functions from a geeglm object
#'
#' Create estimating equation function from a \code{geeglm} object
#'
#' @inheritParams grab_estFUN
#------------------------------------------------------------------------------#

grab_estFUN.geeglm <- function(object, data, ...){
  if(object$corstr != 'independence'){
    stop("only independence working correlation is supported at this time")
  }

  X <- stats::model.matrix(object$formula, data = data)
  # DO NOT use stats::model.matrix(geepack_obj, data = subdata)) -
  # returns entire model matrix, not just the subset
  Y  <- stats::model.response(stats::model.frame(object, data = data))
  n  <- length(Y)
  p  <- length(stats::coef(object))
  phi    <- as.numeric(summary(object)$dispersion[1])
  family <- object$family$family
  link   <- object$family$link
  invlnk <- object$family$linkinv
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
#' Grab estimating functions from a merMod object
#'
#' Create estimating equation function from a \code{merMod} object
#'
#' @param numderiv_opts a list of argument passed to \code{numDeriv::grad}
#' @inheritParams grab_estFUN
#' @export
#------------------------------------------------------------------------------#

# setMethod("grab_estFUN", "merMod", function(object, data, numderiv_opts = NULL,...)
grab_estFUN.merMod <- function(object, data, numderiv_opts = NULL,...)
{
  ## Warnings ##
  if(length(lme4::getME(object, 'theta')) > 1){
    stop('make_eefun.merMod currently does not handle >1 random effect')
  }

  fm     <- grab_fixed_formula(model = object)
  X      <- grab_design_matrix(data = data, rhs_formula = fm)
  Y      <- grab_response(data = data, formula = stats::formula(object))
  family <- object@resp$family
  lnkinv <- family$linkinv
  objfun <- objFun_merMod(family$family)

  function(theta){
    args <- list(func = objfun, x = theta, response = Y, xmatrix = X, linkinv = lnkinv)
    do.call(numDeriv::grad, args = append(args, numderiv_opts))
  }
}

#------------------------------------------------------------------------------#
# Objective Function for merMod object
#
# @param family distribution family of objective function
# @param ... additional arguments pass to objective function
#------------------------------------------------------------------------------#

objFun_merMod <- function(family, ...){
  switch(family,
         binomial = objFun_glmerMod_binomial,
         stop('Objective function for this link/family not defined'))
}

#------------------------------------------------------------------------------#
# Objective Function for Logistic-Normal Likelihood
#
# @param parms vector of parameters
# @param response vector of response values
# @param xmatrix the matrix of covariates
# @param linkinv inverse link function
#------------------------------------------------------------------------------#

objFun_glmerMod_binomial <- function(parms, response, xmatrix, linkinv)
{
  log(stats::integrate(binomial_integrand, lower = -Inf, upper = Inf,
                parms    = parms,
                response = response,
                xmatrix  = xmatrix,
                linkinv  = linkinv)$value)

}

#------------------------------------------------------------------------------#
# Objective Function for Logistic-Normal Likelihood
#
# @inheritParams objFun_glmerMod_binomial
# @param b the random effect to integrate over
#------------------------------------------------------------------------------#

binomial_integrand <- function(b, response, xmatrix, parms, linkinv){
  if(class(xmatrix) != 'matrix'){
    xmatrix <- as.matrix(xmatrix)
  }
  pr  <- linkinv( drop(outer(xmatrix %*% parms[-length(parms)], b, '+') ) )
  hh  <- stats::dbinom(response, 1, prob = pr)
  hha <- apply(hh, 2, prod)
  hha * stats::dnorm(b, mean = 0, sd = parms[length(parms)])
}

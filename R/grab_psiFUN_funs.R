#------------------------------------------------------------------------------#
# grab_psiFUN description:
# a generic function that takes a model object and "grabs" the inner estFUN
# from the object.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Grab estimating functions from a model object
#'
#' @param object the object from which to extrace \code{psiFUN}
#' @param data the data to use for the estimating function
#' @param ... additonal arguments passed to other methods
#' @docType methods
#'
#' @export
#' @examples
#'
#' \dontrun{
#' library(geepack)
#' data('ohio')
#'
#' glmfit  <- glm(resp ~ age, data = ohio,
#'                family = binomial(link = "logit"))
#' geefit  <- geeglm(resp ~ age, data = ohio, id = id,
#'                   family = binomial(link = "logit"))
#' glmmfit <- glmer(resp ~ age + (1|id), data = ohio,
#'                  family = binomial(link = "logit"))
#' example_ee <- function(data, model){
#'  f <- grab_psiFUN(model, data)
#'  function(theta){
#'   f(theta)
#'  }
#' }
#'
#' m_estimate(
#'   estFUN = example_ee,
#'   data = ohio,
#'   compute_roots = FALSE,
#'   units = 'id',
#'   roots = coef(glmfit),
#'   outer_args = list(model = glmfit))
#' m_estimate(
#'   estFUN = example_ee,
#'   data = ohio,
#'   compute_roots = FALSE,
#'   units = 'id',
#'   roots = coef(geefit),
#'   outer_args = list(model = geefit))
#' m_estimate(
#'   estFUN = example_ee,
#'   data = ohio,
#'   compute_roots = FALSE,
#'   units = 'id',
#'   roots = unlist(getME(glmmfit, c('beta', 'theta'))),
#'   outer_args = list(model = glmmfit))
#'  }
#------------------------------------------------------------------------------#

grab_psiFUN <- function(object, ...){
  # S3 generic, for S3 dispatch
  UseMethod("grab_psiFUN")
}


#' @inheritParams grab_psiFUN
#' @param weights a scalar or vector of weight values
#' @describeIn grab_psiFUN Create estimating equation function from a \code{glm} object
#' @export

grab_psiFUN.glm <- function(object, data, weights = 1, ...){

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
      if (n==1) { D <- t(D) } ## apply will undesireably coerce to vector
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }

    t(D) %*% solve(V) %*% diag(W, nrow = n, ncol = n) %*% (r)
  }
}

#' @inheritParams grab_psiFUN
#' @describeIn grab_psiFUN Create estimating equation function from a \code{geeglm} object
#' @export

grab_psiFUN.geeglm <- function(object, data, ...){
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
      if (n==1) { D <- t(D) } ## apply will undesireably coerce to vector
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }
    t(D) %*% solve(V) %*% (r)
  }
}

#' @param numderiv_opts a list of arguments passed to \code{numDeriv::grad}
#' @inheritParams grab_psiFUN
#' @describeIn grab_psiFUN Create estimating equation function from a \code{merMod} object
#' @export

grab_psiFUN.merMod <- function(object, data, numderiv_opts = NULL,...)
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

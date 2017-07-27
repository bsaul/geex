#------------------------------------------------------------------------------#
# xxx_** description:
# The primary function of geex. Performs M-estimation in finding roots and/or
# estimates parameter covariance matrix
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Estimate parameters and their covariance from a set of estimating equations
#'
#' @param eeFUN a function that takes in group-level data and returns a function
#' that takes parameters as its first argument
#' @param data a data.frame
#' @param units a string identifying the grouping variable in \code{data}
#' @param outer_eeargs a list of arguments passed to the outer (data) function of \code{eeFUN}. (optional)
#' @param inner_eeargs a list of arguments passed to the inner (theta) function of \code{eeFUN}. (optional)
#' @param corrections_list an optional list of small sample corrections where each
#' list element is a list with two elements: `fun` and `options`. See details.
#' @param compute_roots whether or not to find the roots of the estimating equations.
#' Defaults to \code{TRUE}.
#' @param compute_vcov whether or not to compute the variance-covariance matrix.
#' Defaults to \code{TRUE}.
#' @param roots a numeric vector containing either starting values for the roots when using
#' the default \code{rootsolver} or roots that have been estimated elsewhere
#' @param rootFUN_object the name of the object within the output of \code{rootFUN} that contains
#' the parameters estimates. For example, the 'root' object within the output of multiroot::rootSolve
#' contains the parameter estimates. Defaults to 'root'. Can also be a set of numeric positions within
#' the object.
#' @inheritParams compute_eeroot
#' @inheritParams compute_matrices
#' @return a list with the following
#' \itemize{
#' \item \code{parameters} - a vector of estimated parameters
#' \item \code{vcov} - the variance-covariance matrix for the parameters
#' \item \code{corrections} - a list of corrected variance-covariance matrices
#' }
#'
#' @section eeFUN arguments:
#'
#' Additional arguments may be passed to both the inner and outer function of the `eeFUN`.
#' Elements in an \code{outer_eeargs} listare passed to the outer function; any elements of the \code{inner_eeargs} list
#' are passed to the inner function. For example, a practical example might be computing a
#' counterfactual mean using an IPW estimator:
#'
#' \preformatted{
#' myeeFUN <- function(data, model){
#'   X <- model.matrix(model, data = data) #covariates
#'   A <- data$A #treatment
#'   Y <- data$Y #outcome
#'   p <- ncol(X) #number of parameters in model
#'   function(theta, a){
#'     Y * (A == a) * 1/plogis(X \%*\% theta[p - 1]) - theta[p]
#'     # Here theta[p] is the target parameter.
#'   }
#' }
#' }
#'
#' Then to estimate the mean where `a == 1`:
#'
#' \preformatted{
#' estimate_equations(
#'   eeFUN = myeeFUN,
#'   data  = mydata,
#'   units = myunits,
#'   inner_eeargs = list(a = 1),
#'   outer_eeargs = list(model = mymodel)
#' )
#' }
#'
#' @export
#------------------------------------------------------------------------------#

estimate_equations <- function(eeFUN,
                               data,
                               units             = NULL,
                               weights           = NULL,
                               roots             = NULL,
                               outer_eeargs      = NULL,
                               inner_eeargs      = NULL,
                               compute_roots     = TRUE,
                               compute_vcov      = TRUE,
                               corrections_list  = NULL,
                               derivFUN          = numDeriv::jacobian,
                               derivFUN_control  = list(method = 'Richardson'),
                               rootFUN           = rootSolve::multiroot,
                               rootFUN_control   = NULL,
                               rootFUN_object    = 'root',
                               approxFUN         = NULL,
                               approxFUN_control = NULL){

  # Split data frame into data frames for each independent unit
  if(is.null(units)){
    # if units are not specified, split into one per observation
    split_data <- split(x = data, f = 1:nrow(data) )
    message('When units are not specified, each observation is considered independent.')
  } else {
    split_data <- split(x = data, f = data[[units]] )
  }

  geex_list  <- list(
    eeFUN        = eeFUN,
    splitdt      = split_data,
    inner_eeargs = inner_eeargs,
    outer_eeargs = outer_eeargs,
    weights      = weights)

  ## Checks/Warnings ##
  if(is.null(roots) & !compute_roots){
    stop('If findroots = FALSE, estimates for the roots must be specified in the roots argument.')
  }

  if(!is.null(corrections_list)){
    check_corrections(corrections_list)
  }

  if(!is.null(weights)){
    if(length(weights) != length(split_data)){
      stop("Length of the weights vector is not equal to the number of units. Check the weights, data, and units arguments.")
    }
  }

  check_eeFUN(geex_list)

  out <- list()
  ## Compute estimating equation roots ##
  if(compute_roots == TRUE){
    eesolved <- compute_eeroot(
      geex_list       = geex_list,
      rootFUN         = rootFUN,
      rootFUN_control = rootFUN_control,
      approxFUN       = approxFUN_control)
    out$rootFUN_results <- eesolved
    out$estimates <- theta_hat <- eesolved[[rootFUN_object]]
  } else {
    out$estimates <- theta_hat <- roots
  }
  if (compute_vcov == FALSE){
    return(list(estimates = theta_hat))
  }

  ## Compute core matrices ##
  mats <- compute_matrices(
    geex_list         = geex_list,
    theta             = theta_hat,
    derivFUN          = derivFUN,
    derivFUN_control  = derivFUN_control,
    approxFUN         = approxFUN,
    approxFUN_control = approxFUN_control)

  ## Compute corrections ##
  if(!is.null(corrections_list)){
    out$corrections <- make_corrections(mats, corrections_list)
  }

  ## Compute covariance estimate(s) ##
  out$vcov <- compute_sigma(A = mats$A, B = mats$B)
  out$component_matrices <- mats

  ## Additional components of output ##
  out$split_data <- split_data

  out
}


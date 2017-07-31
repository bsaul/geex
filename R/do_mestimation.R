#------------------------------------------------------------------------------#
# xxx_** description:
# The primary function of geex. Performs M-estimation in finding roots and/or
# estimates parameter covariance matrix
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Estimate parameters and their covariance from a set of estimating equations
#'
#' @description
#' M-estimation theory provides a framework for asympotic properties of estimators
#' that are solutions to estimating equations. Many R packages implement specific
#' applications of estimating equations. \pkg{geex} aims to be provide a more general
#' framework that any modelling method can use to compute point and variance estimates
#' for parameters that are solutions to estimating equations of the form:
#' \deqn{\sum_i \psi(O_i, \hat{\theta}) = 0}{\sum_i \psi(O_i, \theta) = 0}
#'
#' @param eeFUN a function that takes in group-level data and returns a function
#' that takes parameters as its first argument
#' @param data a data.frame
#' @param units an optional character string identifying the grouping variable in \code{data}
#' @param outer_eeargs a list of arguments passed to the outer (data) function of \code{eeFUN}. (optional)
#' @param inner_eeargs a list of arguments passed to the inner (theta) function of \code{eeFUN}. (optional)
#' @param corrections_list an optional list of small sample corrections where each
#' list element is a list with two elements: \code{correctFUN} and \code{correctFUN_options}.
#' See details.
#' @param compute_roots whether or not to find the roots of the estimating equations.
#' Defaults to \code{TRUE}.
#' @param compute_vcov whether or not to compute the variance-covariance matrix.
#' Defaults to \code{TRUE}.
#' @param roots a vector of parameter estimates must be provided if \code{compute_roots = FALSE}
#' @inheritParams compute_eeroot
#' @inheritParams compute_matrices
#'
#' @details The basic idea of \pkg{geex} is for the analyst to provide at least
#' two items:
#' \itemize{
#' \item data
#' \item \code{eeFUN}: (the \eqn{\psi} function), a function that takes unit-level
#' data and returns a function in terms of parameters (\eqn{\theta})
#' }
#'
#' With the \code{eeFUN}, \pkg{geex} computes the roots of the estimating equations
#' and/or the empirical sandwich variance estimator.
#'
#' The root finding algorithm defaults to \code{\link[rootSolve]{multiroot}} to
#' estimate roots though the solver algorithm can be specified in the \code{rootFUN}
#' argument. Starting values for \code{\link[rootSolve]{multiroot}} are passed via the
#' \code{rootFUN_control} argument. See \code{vignette("03_root_solvers", package = "geex")}
#' for information on customizing the root solver function.
#'
#' To compute only the covariance matrix, set \code{compute_roots = FALSE} and pass
#' estimates of \eqn{\theta} via the \code{roots} argument.
#'
#' M-estimation is often used for clustered data, and a variable by which to split
#' the data.frame  into independent units is specified by the \code{units} argument.
#' This argument defaults to \code{NULL}, in which case the number of units equals
#' the number of rows in the data.frame.
#'
#' For information on the finite-sample corrections, refer to the finite sample
#' correction API vignette: \code{vignette("05_finite_sample_corrections", package = "geex")}
#'
#' @section Writing an eeFUN:
#'
#' \subsection{Description}{
#' An \code{eeFUN} is a function that takes in *unit* level data plus possible
#' "outer" arguments (see section on eefUN argument) and returns a function
#' whose first argument is \code{theta}. See the examples below or the package
#' vignettes for more information.
#' }
#'
#' \subsection{Additional arguments}{
#' Additional arguments may be passed to both the inner and outer function of the
#' \code{eeFUN}. Elements in an \code{outer_eeargs} list are passed to the outer
#' function; any elements of the \code{inner_eeargs} list are passed to the inner
#' function. For an example, see the finite sample correction vignette [\code{
#' vignette("05_finite_sample_corrections", package = "geex")}].
#' }
#'
#'@return a list with the following
#' \itemize{
#' \item \code{parameters} - a vector of estimated parameters
#' \item \code{vcov} - the variance-covariance matrix for the parameters
#' \item \code{corrections} - a list of corrected variance-covariance matrices
#' }
#' @references Stefanski, L. A., & Boos, D. D. (2002). The calculus of M-estimation.
#' The American Statistician, 56(1), 29-38.
#'
#' @examples
#' # Estimate the mean and variance of Y1 in the geexex dataset
#' ex_eeFUN <- function(data){
#'  function(theta){
#'    with(data,
#'     c(Y1 - theta[1],
#'      (Y1 - theta[1])^2 - theta[2] ))
#' }}
#'
#' estimate_equations(
#'  eeFUN = ex_eeFUN,
#'  data  = geexex,
#'  rootFUN_control = list(start = c(1,1)))
#'
#' # A simple linear model for regressing X1 and X2 on Y4
#' lm_eefun <- function(data){
#'  X <- cbind(1, data$X1, data$X2)
#'  Y <- data$Y4
#'   function(theta){
#'     t(X) %*% (Y - X %*% theta)
#'    }
#'  }
#'
#' estimate_equations(
#'  eeFUN = lm_eefun,
#'  data  = geexex,
#'  rootFUN_control = list(start = c(0, 0, 0)))
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


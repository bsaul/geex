#------------------------------------------------------------------------------#
# estimate_** description:
# Functions that estimate quantities from an m_estimation_basis *or* are higher
# level wrappers for estimate_** functions (e.g. m_estimate())
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Estimates roots for a set of estimating equations
#'
#' Using the \code{rootFUN} specified by the user (defaults to \code{\link[rootSolve]{multiroot}}),
#' this function estimates the roots of the equations:
#' \deqn{G_m = sum_i psi(O_i, \hat{\theta}) = 0}{G_m = sum_i psi(O_i, theta) = 0}
#'
#' This is primilary an internal function used within \code{\link{estimate_equations}},
#' but it is exported for use in debugging and development.
#'
#' For an example of how to use a different \code{rootFUN},
#' see the root solver vignette, \code{vignette('geex_root_solvers', package = 'geex')}.
#'
#' @param basis an object of class \code{\linkS4class{m_estimation_basis}}
#' @param rootFUN the function used to find roots of the estimating equations.
#' Defaults to \code{\link[rootSolve]{multiroot}}.
#' @param rootFUN_control a list of options to be passed to the \code{rootsolver}
#' function
#' @inheritParams create_psi
#' @return the output of the \code{rootFUN} function
#' @export
#'
#------------------------------------------------------------------------------#

estimate_GFUN_roots <- function(basis,
                              rootFUN           = rootSolve::multiroot,
                              rootFUN_control   = NULL,
                              approxFUN         = NULL,
                              approxFUN_control = NULL){

  rootFUN <- match.fun(rootFUN)

  # Create estimating equation functions per group
  psi_i <- create_psi(.split_data   = basis@.split_data,
                      .estFUN       = grab_estFUN(basis),
                      .outer_estFUN_args = basis@.outer_args)

  # Create psi function that sums over all ee funs
  # G_m = sum_i psi(O_i, theta) in SB notation]
  GmFUN <- create_GFUN(.psi_list     = psi_i,
                       .inner_estFUN_args = basis@.inner_args,
                       .weights      = basis@.weights)

  # Find roots of psi
  rargs <- append(rootFUN_control, list(f = GmFUN))
  do.call(rootFUN, args = rargs)
}

#------------------------------------------------------------------------------#
#' Estimate component matrices of the empirical sandwich covariance estimator
#'
#' For a given set of estimating equations computes the 'meat' (\eqn{B_m}{B_m}
#' in Stefanski and Boos notation) and 'bread' (\eqn{A_m}{A_m} in Stefanski and
#'  Boos notation) matrices necessary to compute the covariance matrix.
#'
#' @param basis basis an object of class \code{\linkS4class{m_estimation_basis}}
#' @param theta vector of parameter estimates (i.e. estimated roots) passed to \code{geex_list$eeFUN}.
#' @param derivFUN the function used to take derivatives of the estimating equation functions.
#' Defaults to \code{\link[numDeriv]{jacobian}}.
#' @param derivFUN_control a list of options passed to \code{\link[numDeriv]{jacobian}}
#' (or the \code{derivFUN} function).
#' @inheritParams create_psi
#'
#' @return a list with
#' \itemize{
#' \item A - the 'bread' matrix
#' \item B - the 'meat' matrix
#' \item A_i - a list of 'bread' matrices for each group
#' \item B_i - a list of 'meat' matrices for each group
#' }
#' @details For a set of estimating equations (\eqn{\sum_i \psi(O_i, \theta) = 0}{sum_i \psi(O_i, \theta) = 0}),
#' this function computes:
#'
#' \deqn{A_i =  \partial \psi(O_i, \theta)/\partial \theta}{A_i =  partial \psi(O_i, theta)/\partial \theta}
#'
#' \deqn{A =  \sum_i A_i}{A =  \sum_i A_i}
#'
#' \deqn{B_i =  \psi(O_i, \theta)\psi(O_i, \theta)^T}{B_i = outer(\psi(O_i, \theta), \psi(O_i, \theta))}
#'
#' \deqn{B =  \sum_i B_i}{B =  \sum_i B_i}
#'
#' where all of the above are evaluated at \eqn{\hat{\theta}}{hat(\theta)}. The partial derivatives in \eqn{A_i}{A_i}
#' numerically approximated by the \code{derivFUN}.
#'
#' Note that \eqn{A =  \sum_i A_i}{A =  \sum_i A_i} and not \eqn{\sum_i A_i/m}{A =  \sum_i A_i/m}, and the same for \eqn{B}{B}.
#'
#' @export
#' @references Stefanski, L. A., & Boos, D. D. (2002). The calculus of m-estimation. The American Statistician, 56(1), 29-38.
#------------------------------------------------------------------------------#

estimate_sandwich_matrices <- function(basis,
                                       theta,
                                       derivFUN         = numDeriv::jacobian,
                                       derivFUN_control = list(method = 'Richardson'),
                                       approxFUN        = NULL,
                                       approxFUN_control = NULL){

  derivFUN <- match.fun(derivFUN)

  # Create list of estimating eqn functions per unit
  psi_i <- create_psi(.split_data = basis@.split_data,
                      .estFUN     = grab_estFUN(basis),
                      .outer_estFUN_args = basis@.outer_args,
                      approxFUN = approxFUN,
                      approxFUN_control = approxFUN_control)

  # Compute the negative of the derivative matrix of estimating eqn functions
  # (the information matrix)
  A_i <- lapply(psi_i, function(ee){
    args <- append(list(fun = ee, x = theta), derivFUN_control)
    val  <- do.call(derivFUN, args = append(args, basis@.inner_args))
    -val
  })

  w <- basis@.weights
  A_i_pre   <- if(length(w) > 0){ Map(`*`, A_i, w) } else A_i
  A_i_array <- check_array(simplify2array(A_i_pre))
  A   <- apply(A_i_array, 1:2, sum)

  # Compute outer product of observed estimating eqns
  B_i <- lapply(psi_i, function(ee) {
    ee_val <- do.call(ee, args = append(list(theta = theta), basis@.inner_args))
    ee_val %*% t(ee_val)
  })
  B_i_pre   <- if(length(w) > 0){ Map(`*`, B_i, w) } else B_i
  B_i_array <- check_array(simplify2array(B_i_pre))
  B   <- apply(B_i_array, 1:2, sum)

  list(A = A, A_i = A_i, B = B, B_i = B_i)
}


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
#' @inheritParams estimate_GFUN_roots
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

m_estimate <- function(estFUN,
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



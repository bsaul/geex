#------------------------------------------------------------------------------#
#' Creates list of psi functions
#'
#' @param splitdt list of dataframes with data per unit
#' @param eeFUN the estimating equation function
#' @param approxFUN a function that approximates the inner function of \code{eeFUN}
#' @param approxFUN_control arguments passed to \code{approxFUN}
#' @param outer_eeargs a list of arguments passed to \code{eeFUN}
#' @export
#'
#------------------------------------------------------------------------------#

create_psi <- function(splitdt,
                       eeFUN,
                       approxFUN = NULL,
                       approxFUN_control = NULL,
                       outer_eeargs = NULL){
  out <- lapply(splitdt, function(data_i){
    do.call(eeFUN, args = append(list(data = data_i), outer_eeargs))
  })

  # if user specifies an approximation function, apply the function to each
  # evaluation of psi
  if(!is.null(approxFUN)){
    lapply(out, function(f){
      do.call(approxFUN, args = append(list(psi = f), approxFUN_control))
    }) -> out
  }

  out
}

#------------------------------------------------------------------------------#
#' Creates a function that sums over psi functions
#'
#' @param psi_list list of psi functions
#' @param inner_eeargs list of arguments passed to psi
#' @export
#'
#------------------------------------------------------------------------------#

create_GFUN <- function(psi_list, inner_eeargs = NULL){
  function(theta){
    psii <- lapply(psi_list, function(f) {
      do.call(f, args = append(list(theta = theta), inner_eeargs))
    })
    # sum over unit-wise contributions to the estimating equations
    apply(check_array(simplify2array(psii)), 1, sum)
  }
}

#------------------------------------------------------------------------------#
#' Compute roots for a set of estimating equations
#'
#' @param geex_list a list containing \code{splitdt} (a \code{data.frame} that
#' has been \code{\link[base]{split}} by the grouping variable) and \code{eeFUN}
#' (see \code{\link{estimate_equations}}). Note that the function that is
#' returned by \code{eeFUN} must take \code{theta} as its first argument, where
#' \code{theta} represents the parameters.
#' @param rootFUN the function used to find roots of the estimating equations.
#' Defaults to \code{\link[rootSolve]{multiroot}}.
#' @param rootFUN_control a list of options to be passed to the \code{rootsolver}
#' function
#' @inheritParams create_psi
#' @return the output of the \code{rootsolver} function
#' @export
#------------------------------------------------------------------------------#

compute_eeroot <- function(geex_list,
                   rootFUN           = rootSolve::multiroot,
                   rootFUN_control   = NULL,
                   approxFUN         = NULL,
                   approxFUN_control = NULL){

  rootFUN <- match.fun(rootFUN)

  # Create estimating equation functions per group
  psi_i <- create_psi(splitdt      = geex_list$splitdt,
                      eeFUN        = geex_list$eeFUN,
                      outer_eeargs = geex_list$outer_eeargs)

  # Create psi function that sums over all ee funs
  # G_m = sum_i psi(O_i, theta) in SB notation]
  GmFUN <- create_GFUN(psi_list     = psi_i,
                       inner_eeargs = geex_list$inner_eeargs)

  # Find roots of psi
  rargs <- append(rootFUN_control, list(f = GmFUN))
  do.call(rootFUN, args = rargs)
}

#------------------------------------------------------------------------------#
#' Compute component matrices for covariance matrix
#'
#' For a given set of estimating equations computes the 'meat' and 'bread' matrices
#' necessary to compute the covariance matrix.
#'
#' @param geex_list a list containing \code{splitdt} (a \code{data.frame} that
#' has been \code{\link[base]{split}} by the grouping variable) and \code{eeFUN}
#' (see \code{\link{estimate_equations}})
#' @param theta vector of parameters passed to \code{geex_list$eeFUN}.
#' @param derivFUN the function used to take derivatives of the estimating equation functions.
#' Defaults to \code{numDeriv::jacobian}
#' @param derivFUN_control a list of options passed to \code{\link[numDeriv]{jacobian}}
#' (or the \code{derivFUN} function).
#' @inheritParams create_psi
#' @return a list with
#' \itemize{
#' \item A - the 'bread' matrix
#' \item B - the 'meat' matrix
#' \item A_i - the 'bread' matrix for each group
#' \item B_i - the 'meat' matrix for each group
#' }
#'
#' @export
#------------------------------------------------------------------------------#

compute_matrices <- function(geex_list,
                             theta,
                             derivFUN         = numDeriv::jacobian,
                             derivFUN_control = list(method = 'Richardson'),
                             approxFUN        = NULL,
                             approxFUN_control = NULL){

  derivFUN <- match.fun(derivFUN)

  if(is.null(geex_list$ee_args)){
    ee_args <- NULL
  }

  # Create list of estimating eqn functions per unit
  psi_i <- create_psi(splitdt = geex_list$splitdt,
                      eeFUN   = geex_list$eeFUN,
                      outer_eeargs = geex_list$outer_eeargs,
                      approxFUN = approxFUN,
                      approxFUN_control = approxFUN_control)

  # Compute the negative of the derivative matrix of estimating eqn functions
  # (the information matrix)
  A_i <- lapply(psi_i, function(ee){
    args <- append(list(fun = ee, x = theta), derivFUN_control)
    val  <- do.call(derivFUN, args = append(args, geex_list$inner_eeargs))
    -val
  })
  A_i_array <- check_array(simplify2array(A_i))
  A   <- apply(A_i_array, 1:2, sum)

  # Compute outer product of observed estimating eqns
  B_i <- lapply(psi_i, function(ee) {
    ee_val <- do.call(ee, args = append(list(theta = theta), geex_list$inner_eeargs))
    ee_val %*% t(ee_val)
  })
  B   <- apply(check_array(simplify2array(B_i)), 1:2, sum)

  list(A = A, A_i = A_i, B = B, B_i = B_i)
}

#------------------------------------------------------------------------------#
#' Compute covariance matrix for set of estimating equations
#'
#' Computes \deqn{\Sigma = A_m^{-1} B_m  (A_m^{-1})^T / m } where
#'
#' \eqn{A_m = \sum_{i = 1}^m A_i}{A_m = sum_i A_i} and \eqn{B_m = \sum_{i = 1}^m B_i}{A_m = sum_i B_i}.
#'
#' @param A the `A` matrix returned in the list of matrices from
#'   \code{\link{compute_matrices}}
#' @param B the `B` matrix returned in the list of matrices from
#'   \code{\link{compute_matrices}}
#'
#' @export
#------------------------------------------------------------------------------#

compute_sigma <- function(A, B){
  Ainv <- solve(A)
  Ainv %*% B %*% t(Ainv)
}

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
    outer_eeargs = outer_eeargs)

  ## Checks/Warnings ##
  if(is.null(roots) & !compute_roots){
    stop('If findroots = FALSE, estimates for the roots must be specified in the roots argument.')
  }

  if(!is.null(corrections_list)){
    check_corrections(corrections_list)
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


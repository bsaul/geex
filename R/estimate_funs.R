#------------------------------------------------------------------------------#
# estimate_** description:
# Functions that estimate quantities from an m_estimation_basis *or* are higher
# level wrappers for estimate_** functions (e.g. m_estimate()). The workhorses
# of performing M-estimation
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Estimate roots for a set of estimating equations
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
#' @param root_control an object of class \code{\linkS4class{root_control}}
#' @return the output of the \code{rootFUN} function
#' @export
#'
#------------------------------------------------------------------------------#

estimate_GFUN_roots <- function(.GFUN,
                                .root_control,
                                .approx_control){

  # Use root_control defaults if no options passed
  if(missing(.root_control)){
    .root_control <- new('root_control')
  }

  # Use approx_control defaults if no options passed
  # TODO: this code is in create_psi() as well, can it be in 1 place instead?
  if(missing(.approx_control)){
    .approx_control <- new('approx_control')
  }

  rootFUN <- match.fun(FUN(.root_control))

  # Find roots of psi
  rargs <- append(options(.root_control), list(f = .GFUN))
  do.call(rootFUN, args = rargs)
}

#------------------------------------------------------------------------------#
# Process a list of matrices to sum across them
#
# @param l a list of matrices
# @param w a numeric vector of weights
#------------------------------------------------------------------------------#

process_matrix_list <- function(.l, .w = numeric(0)){
  M_i_pre   <- if(length(.w) > 0){ Map(`*`, .l, .w) } else .l
  M_i_array <- check_array(simplify2array(M_i_pre))
  apply(M_i_array, 1:2, sum)
}

#------------------------------------------------------------------------------#
#' Estimate component matrices of the empirical sandwich covariance estimator
#'
#' For a given set of estimating equations computes the 'meat' (\eqn{B_m}{B_m}
#' in Stefanski and Boos notation) and 'bread' (\eqn{A_m}{A_m} in Stefanski and
#'  Boos notation) matrices necessary to compute the covariance matrix.
#'
#' @param .basis basis an object of class \code{\linkS4class{m_estimation_basis}}
#' @param .theta vector of parameter estimates (i.e. estimated roots)
#' @param .deriv_control an object of class \code{\linkS4class{deriv_control}}
#' @inheritParams create_GFUN
#'
#' @return a \code{\linkS4class{sandwich_components}} object
#'
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

estimate_sandwich_matrices <- function(.basis,
                                       .theta,
                                       .deriv_control,
                                       .approx_control){
  # Use deriv_control defaults if no options passed
  if(missing(.deriv_control)){
    .deriv_control <- new('deriv_control')
  }
  # Use approx_control defaults if no options passed
  # TODO: this code is in create_psi() as well, can it be in 1 place instead?
  if(missing(.approx_control)){
    .approx_control <- new('approx_control')
  }

  derivFUN <- match.fun(FUN(.deriv_control))
  w <- .basis@.weights
  psi_list <- grab_psiFUN_list(.basis)
  # Compute the negative of the derivative matrix of estimating eqn functions
  # (the information matrix)
  A_i <- lapply(psi_list, function(ee){
    args <- append(list(func = ee, x = .theta), options(.deriv_control))
    val  <- do.call(derivFUN, args = append(args, .basis@.inner_args))
    -val
  })

  A <- process_matrix_list(A_i, w)

  # Compute outer product of observed estimating eqns
  B_i <- lapply(psi_list, function(ee) {
    ee_val <- do.call(ee, args = append(list(theta = .theta), .basis@.inner_args))
    ee_val %*% t(ee_val)
  })

  B <- process_matrix_list(B_i, w)

  new('sandwich_components',
      .A = A, .A_i = A_i, .B = B, .B_i = B_i)
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
#' @param corrections an optional list of small sample corrections where each
#' list element is a \code{\linkS4class{correct_control}} object which contains
#' two elements: \code{correctFUN} and \code{correctFUN_options}. The function
#' \code{\link{correction}} constructs \code{\linkS4class{correct_control}} objects.
#' See details for more information.
#' @param compute_roots whether or not to find the roots of the estimating equations.
#' Defaults to \code{TRUE}.
#' @param compute_vcov whether or not to compute the variance-covariance matrix.
#' Defaults to \code{TRUE}.
#' @param roots a vector of parameter estimates must be provided if \code{compute_roots = FALSE}
#' @inheritParams estimate_GFUN_roots
#' @inheritParams estimate_sandwich_matrices
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
                       units             = character(0),
                       weights           = numeric(0),
                       outer_args        = list(),
                       inner_args        = list(),
                       roots             = NULL,
                       compute_roots     = TRUE,
                       compute_vcov      = TRUE,
                       corrections,
                       deriv_control,
                       root_control,
                       approx_control){

  # # Split data frame into data frames for each independent unit
  # if(is.null(units)){
  #   # if units are not specified, split into one per observation
  #   split_data <- split(x = data, f = 1:nrow(data) )
  #   message('When units are not specified, each observation is considered independent.')
  # } else {
  #   split_data <- split(x = data, f = data[[units]] )
  # }

  # Use deriv_control defaults if no options passed
  if(missing(deriv_control)){
    deriv_control <- new('deriv_control')
  }
  # Use root_control defaults if no options passed
  if(missing(root_control)){
    root_control <- new('root_control')
  }
  # Use approx_control defaults if no options passed
  # TODO: this code is in create_psi() as well, can it be in 1 place instead?
  if(missing(approx_control)){
    approx_control <- new('approx_control')
  }

  basis <- new("m_estimation_basis",
               .estFUN     = estFUN,
               .data       = data,
               .units      = units,
               .weights    = weights,
               .outer_args = outer_args,
               .inner_args = inner_args)

  ## Checks/Warnings ##
  if(is.null(roots) & !compute_roots){
    stop('If findroots = FALSE, estimates for the roots must be specified in the roots argument.')
  }

  if(length(basis@.weights) > 0){
    if(length(basis@.weights) != length(basis@.split_data)){
      stop("Length of the weights vector is not equal to the number of units. Check the weights, data, and units arguments.")
    }
  }

  # Create estimating equation functions per group
  basis  <- create_psiFUN_list(.basis          = basis,
                               .approx_control = approx_control)

  # Create psi function that sums over all ee funs
  # G_m = sum_i psi(O_i, theta) in SB notation]
  GmFUN <- create_GFUN(.basis = basis)

  ## Compute estimating equation roots ##
  if(compute_roots == TRUE){
    eesolved <- estimate_GFUN_roots(
      .GFUN           = GmFUN,
      .root_control   = root_control,
      .approx_control = approx_control)
    rootFUN_results <- eesolved
    theta_hat <- eesolved[[root_control@.object_name]]
  } else {
    rootFUN_results <- list()
    theta_hat <- roots
  }

  ## Compute component matrices ##
  if(compute_vcov == TRUE){
    mats <- estimate_sandwich_matrices(
      .basis             = basis,
      .theta             = theta_hat,
      .deriv_control     = deriv_control,
      .approx_control    = approx_control)

    ## Compute corrections ##
    if(!missing(corrections)){
      correction_results <- make_corrections(mats, corrections)
    } else {
      correction_results <- list()
    }

    ## Compute covariance estimate(s) ##
    vcov <- compute_sigma(A = grab_bread(mats), B = grab_meat(mats))
  } else {
    mats <- new('sandwich_components')
    correction_results <- list()
  }

  out <- new('geex',
             basis           = basis,
             root_control    = root_control,
             approx_control  = approx_control,
             deriv_control   = deriv_control,
             rootFUN_results = rootFUN_results,
             GFUN            = GmFUN,
             sandwich_components = mats,
             corrections     = correction_results,
             estimates       = theta_hat,
             vcov            = vcov)

  out
}



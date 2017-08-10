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
#' This is primilary an internal function used within \code{\link{m_estimate}},
#' but it is exported for use in debugging and development.
#'
#' For an example of how to use a different \code{rootFUN},
#' see the root solver vignette, \code{vignette('geex_root_solvers', package = 'geex')}.
#'
#' @param .basis an object of class \code{\linkS4class{m_estimation_basis}}
#' @return the output of the \code{rootFUN} function
#' @export
#'
#------------------------------------------------------------------------------#

estimate_GFUN_roots <- function(.basis){
  GFUN    <- grab_GFUN(.basis)
  rootFUN <- match.fun(grab_FUN(.basis@.control, "root"))
  rargs   <- append(grab_options(.basis@.control, "root"), list(f = GFUN))
  do.call(rootFUN, args = rargs)
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
#'
#' @return a \code{\linkS4class{sandwich_components}} object
#'
#' @details For a set of estimating equations (\eqn{\sum_i \psi(O_i, \theta) = 0}{sum_i \psi(O_i, \theta) = 0}),
#' this function computes:
#'
#' \deqn{A_i =  \partial \psi(O_i, \theta)/\partial \theta}{A_i =  \partial \psi(O_i, \theta)/\partial \theta}
#'
#' \deqn{A =  \sum_i A_i}{A =  \sum_i A_i}
#'
#' \deqn{B_i =  \psi(O_i, \theta)\psi(O_i, \theta)^T}{B_i = outer(\psi(O_i, \theta), \psi(O_i, \theta))}
#'
#' \deqn{B =  \sum_i B_i}{B =  \sum_i B_i}
#'
#' where all of the above are evaluated at \eqn{\hat{\theta}}{hat(\theta)}. The partial derivatives in \eqn{A_i}{A_i}
#' numerically approximated by the function defined in \code{\linkS4class{deriv_control}}.
#'
#' Note that \eqn{A =  \sum_i A_i}{A =  \sum_i A_i} and not \eqn{\sum_i A_i/m}{A =  \sum_i A_i/m}, and the same for \eqn{B}{B}.
#'
#' @export
#' @references Stefanski, L. A., & Boos, D. D. (2002). The calculus of m-estimation. The American Statistician, 56(1), 29-38.
#------------------------------------------------------------------------------#

estimate_sandwich_matrices <- function(.basis, .theta){

  derivFUN <- match.fun(grab_FUN(.basis@.control, "deriv"))
  derivOPT <- grab_options(.basis@.control, "deriv")
  w        <- .basis@.weights
  psi_list <- grab_psiFUN_list(.basis)

  # Compute the negative of the derivative matrix of estimating eqn functions
  # (the information matrix)
  A_i <- lapply(psi_list, function(ee){
    args <- append(list(func = ee, x = .theta), derivOPT)
    val  <- do.call(derivFUN, args = append(args, .basis@.inner_args))
    -val
  })

  A <- compute_sum_of_list(A_i, w)

  # Compute outer product of observed estimating eqns
  B_i <- lapply(psi_list, function(ee) {
    ee_val <- do.call(ee, args = append(list(theta = .theta), .basis@.inner_args))
    ee_val %*% t(ee_val)
  })

  B <- compute_sum_of_list(B_i, w)

  methods::new('sandwich_components',
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
#' @param estFUN a function that takes in group-level data and returns a function
#' that takes parameters as its first argument
#' @param data a data.frame
#' @param units an optional character string identifying the grouping variable in \code{data}
#' @param weights an optional vector of weights. See details.
#' @param outer_args a list of arguments passed to the outer (data) function of \code{estFUN}. (optional)
#' @param inner_args a list of arguments passed to the inner (theta) function of \code{estFUN}. (optional)
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
#' @param deriv_control a \code{\linkS4class{deriv_control}} object
#' @param root_control a \code{\linkS4class{root_control}} object
#' @param approx_control a \code{\linkS4class{approx_control}} object
#'
#' @details The basic idea of \pkg{geex} is for the analyst to provide at least
#' two items:
#' \itemize{
#' \item data
#' \item \code{estFUN}: (the \eqn{\psi} function), a function that takes unit-level
#' data and returns a function in terms of parameters (\eqn{\theta})
#' }
#'
#' With the \code{estFUN}, \pkg{geex} computes the roots of the estimating equations
#' and/or the empirical sandwich variance estimator.
#'
#' The root finding algorithm defaults to \code{\link[rootSolve]{multiroot}} to
#' estimate roots though the solver algorithm can be specified in the \code{rootFUN}
#' argument. Starting values for \code{\link[rootSolve]{multiroot}} are passed via the
#' \code{root_control} argument. See \code{vignette("v03_root_solvers", package = "geex")}
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
#' correction API vignette: \code{vignette("v05_finite_sample_corrections", package = "geex")}
#'
#' @section Writing an estFUN:
#'
#' \subsection{Description}{
#' An \code{estFUN} is a function representing \eqn{\psi}{\psi}. \pkg{geex} works
#' by breaking \eqn{\psi}{\psi} into two parts:
#'
#' \itemize{
#' \item the "outer" part of the \code{estFUN} which manipulates \code{data} and
#' \code{outer_args} and returns an
#' \item "inner" function of \code{theta} and \code{inner_args}. Internally, this
#' "inner" function is this \code{psiFUN}
#' }
#'
#' In pseudo-code this looks like:
#' \code{
#' function(data, <<outer_args>>){
#'   O <- manipulate(data, <<outer_args>>)
#'   function(theta, <<inner_args>>){
#'     map(O, to = theta, and = <<inner_args>>)
#'   }
#' }}
#'
#' See the examples below or the package vignettes to see this in action.
#' }
#'
#' \subsection{Additional arguments}{
#' Additional arguments may be passed to both the inner and outer function of the
#' \code{estFUN}. Elements in an \code{outer_args} list are passed to the outer
#' function; any elements of the \code{inner_args} list are passed to the inner
#' function. For an example, see the finite sample correction vignette [\code{
#' vignette("v05_finite_sample_corrections", package = "geex")}].
#' }
#'
#' @section Setting up root_control:
#'
#' To estimate roots of the estimating functions, \pkg{geex} uses the \pkg{rootSolve}
#' \code{\link[rootSolve]{mutiroot}} function by default, which requires starting
#' values. The \code{root_control} argument expects a \code{\linkS4class{root_control}}
#' object, which the utility function \code{\link{setup_root_control}} aids in
#' creating. For example, \code{setup_root_control(start = 4)} creates a
#' \code{\linkS4class{root_control}} setting the starting value to 4. In general,
#' the dimension of \code{start} must the same as \code{theta} in the inner
#' \code{estFUN}
#'
#' @section Using weights:
#'
#' In some situations, use of weights can massively speed computations. Refer
#' to \code{vignette("v04_weights", package = "geex")} for an example.
#'
#' @return a \code{\linkS4class{geex}} object
#'
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
#' m_estimate(
#'  estFUN = ex_eeFUN,
#'  data  = geexex,
#'  root_control = setup_root_control(start = c(1,1)))
#'
#' # compare to the mean() and variance() functions
#' mean(geexex$Y1)
#' n <- nrow(geexex)
#' var(geexex$Y1) * (n - 1)/n
#'
#' # A simple linear model for regressing X1 and X2 on Y4
#' lm_eefun <- function(data){
#'  X <- cbind(1, data$X1, data$X2)
#'  Y <- data$Y4
#'  function(theta){
#'     t(X) %*% (Y - X %*% theta)
#'    }
#'  }
#'
#' m_estimate(
#'  estFUN = lm_eefun,
#'  data  = geexex,
#'  root_control = setup_root_control(start = c(0, 0, 0)))
#'
#' # Compare to lm() results
#' summary(lm(Y4 ~ X1 + X2, data = geexex))
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

  control <- methods::new('geex_control')
  if(!missing(deriv_control)){
    set_control(control, 'deriv') <- deriv_control
  }
  if(!missing(root_control)){
    set_control(control, 'root') <- root_control
  }
  if(!missing(approx_control)){
    set_control(control, 'approx') <- approx_control
  }

  basis <- methods::new("m_estimation_basis",
               .estFUN     = estFUN,
               .data       = data,
               .units      = units,
               .weights    = weights,
               .outer_args = outer_args,
               .inner_args = inner_args,
               .control    = control)

  out <- methods::new('geex',
             basis           = basis)
  ## Checks/Warnings ##
  if(is.null(roots) & !compute_roots){
    stop('If findroots = FALSE, estimates for the roots must be specified in the roots argument.')
  }

  ## Compute estimating equation roots ##
  if(compute_roots == TRUE){
    eesolved <- estimate_GFUN_roots(basis)
    out@rootFUN_results <- eesolved
    theta_hat <- eesolved[[root_control@.object_name]]
  } else {
    theta_hat <- roots
  }

  out@estimates <- theta_hat

  ## Compute component matrices ##
  if(compute_vcov == TRUE){
    mats <- estimate_sandwich_matrices(.basis = basis, .theta = theta_hat)
    ## Compute corrections ##
    if(!missing(corrections)){
      out@corrections <- make_corrections(mats, corrections)
    }

    out@sandwich_components <- mats
    ## Compute covariance estimate(s) ##
    out@vcov <- compute_sigma(A = grab_bread(mats), B = grab_meat(mats))
  } else {
    mats <- methods::new('sandwich_components')
  }

  out
}

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

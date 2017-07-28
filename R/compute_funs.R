#------------------------------------------------------------------------------#
# compute_** description:
# Functions that take in stuff and compute numeric values. The workhorses of
# performing M-estimation
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Compute roots for a set of estimating equations
#'
#' Using the \code{rootFUN} specified by the user (defaults to \code{\link[rootSolve]{multiroot}}),
#' this function estimates the roots of the equations: \deqn{G_m = sum_i psi(O_i, \hat{\theta}) = 0}{G_m = sum_i psi(O_i, theta) = 0}
#'
#' This is primilary an internal function used within \code{\link{estimate_equations}},
#' but it is exported for use in debugging and development.
#'
#' For an example of how to use a different \code{rootFUN},
#' see the root solver vignette, \code{vignette('geex_root_solvers', package = 'geex')}.
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
#' @return the output of the \code{rootFUN} function
#' @export
#'
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
                       inner_eeargs = geex_list$inner_eeargs,
                       weights      = geex_list$weights)

  # Find roots of psi
  rargs <- append(rootFUN_control, list(f = GmFUN))
  do.call(rootFUN, args = rargs)
}

#------------------------------------------------------------------------------#
#' Compute component matrices of the empirical sandwich covariance estimator
#'
#' For a given set of estimating equations computes the 'meat' (\eqn{B_m}{B_m} in Stefanski and Boos notation)
#' and 'bread' (\eqn{A_m}{A_m} in Stefanski and Boos notation) matrices necessary to compute the covariance matrix.
#'
#' @param geex_list a list containing \code{splitdt} (a \code{data.frame} that
#' has been \code{\link[base]{split}} by the grouping variable) and \code{eeFUN}
#' (see \code{\link{estimate_equations}})
#' @param theta vector of parameter estimates (i.e. estimated roots) passed to \code{geex_list$eeFUN}.
#' @param derivFUN the function used to take derivatives of the estimating equation functions.
#' Defaults to \code{\link[numDeriv]{jacobian}}.
#' @param derivFUN_control a list of options passed to \code{\link[numDeriv]{jacobian}}
#' (or the \code{derivFUN} function).
#' @inheritParams create_psi
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

  A_i_pre   <- if(!is.null(geex_list$weights)){ Map(`*`, A_i, geex_list$weights) } else A_i
  A_i_array <- check_array(simplify2array(A_i_pre))
  A   <- apply(A_i_array, 1:2, sum)

  # Compute outer product of observed estimating eqns
  B_i <- lapply(psi_i, function(ee) {
    ee_val <- do.call(ee, args = append(list(theta = theta), geex_list$inner_eeargs))
    ee_val %*% t(ee_val)
  })
  B_i_pre   <- if(!is.null(geex_list$weights)){ Map(`*`, B_i, geex_list$weights) } else B_i
  B_i_array <- check_array(simplify2array(B_i_pre))
  B   <- apply(B_i_array, 1:2, sum)

  list(A = A, A_i = A_i, B = B, B_i = B_i)
}

#------------------------------------------------------------------------------#
#' Compute empirical sandwich covariate estimator
#'
#' Computes \eqn{\Sigma = A^{-1} B (A^{-1})^T }{\Sigma = A^{-1} B (A^{-1})^T} with
#' provided \eqn{A}{A} and \eqn{B}{B} matrices.
#'
#' @param A a matrix, generally the \eqn{A}{A} matrix returned in the list of matrices from
#'   \code{\link{compute_matrices}}
#' @param B a matrix, generally the \eqn{B}{B} matrix returned in the list of matrices from
#'   \code{\link{compute_matrices}}
#' @export
#------------------------------------------------------------------------------#

compute_sigma <- function(A, B){
  Ainv <- solve(A)
  Ainv %*% B %*% t(Ainv)
}

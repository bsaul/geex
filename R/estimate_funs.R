#------------------------------------------------------------------------------#
# estimate_** description:
# Functions that estimate quantities from an m_estimation_basis *or* are higher
# level wrappers for estimate_** functions (e.g. m_estimate())
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

estimate_Gm_root <- function(basis,
                             rootFUN           = rootSolve::multiroot,
                             rootFUN_control   = NULL,
                             approxFUN         = NULL,
                             approxFUN_control = NULL){

  rootFUN <- match.fun(rootFUN)

  # Create estimating equation functions per group
  psi_i <- create_psi(.splitdt      = geex_list$splitdt,
                      .estFUN       = grab_estFUN(basis),
                      .outer_eeargs = geex_list$outer_eeargs)

  # Create psi function that sums over all ee funs
  # G_m = sum_i psi(O_i, theta) in SB notation]
  GmFUN <- create_GFUN(psi_list     = psi_i,
                       inner_eeargs = geex_list$inner_eeargs,
                       weights      = geex_list$weights)

  # Find roots of psi
  rargs <- append(rootFUN_control, list(f = GmFUN))
  do.call(rootFUN, args = rargs)
}

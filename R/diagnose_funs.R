#------------------------------------------------------------------------------#
# diagnose_** description:
# Functions that diagnose the output of geex, e.g. diagnose_roots(geex_object).
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Diagnose roots of estimating equations
#'
#' Computes the value of \deqn{G_m = sum_i psi(O_i, \hat{\theta})}{G_m = sum_i psi(O_i, theta)}, i.e., the estimating
#' equations at \code{theta}. Used to verify that \eqn{G_m = 0}{G_m = 0} (or close to 0).
#'
#' @param GFUN a function of theta
#' @param theta parameter estimates to use in evaluating the estimating equations.
#' @return a numeric vector
#' @export
#'
#------------------------------------------------------------------------------#

diagnose_roots <- function(GFUN, theta){
  GFUN(theta)
}



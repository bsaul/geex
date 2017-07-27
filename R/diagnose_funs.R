#------------------------------------------------------------------------------#
# diagnose_** description:
# Functions that diagnose the output of geex, e.g. diagnose_roots(geex_object).
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Evaluate estimating eqations.
#'
#' Computes the value of \deqn{G_m = sum_i psi(O_i, \hat{\theta})}{G_m = sum_i psi(O_i, theta)}, i.e., the estimating
#' equations at \code{theta}. Used to verify that \eqn{G_m = 0}{G_m = 0} (or close to 0).
#'
#' @param geex_list a list with at least \code{splitdt}, \code{eeFUN}, and \code{outer_eeargs} (which can be \code{NULL})
#' @param theta parameter estimates to use in evaluating the estimating equations.
#' @return a numeric vector
#' @export
#'
#------------------------------------------------------------------------------#

diagnose_roots <- function(geex_list, theta){
  psi_i <- create_psi(splitdt      = geex_list$splitdt,
                      eeFUN        = geex_list$eeFUN,
                      outer_eeargs = geex_list$outer_eeargs)
  GFUN <-  create_GFUN(psi_list = psi_i, inner_eeargs = geex_list$inner_eeargs)

  GFUN(theta)
}



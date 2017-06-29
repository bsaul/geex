#------------------------------------------------------------------------------#
#' Check value of estimating eqations
#'
#' Checks that the value of Gm = sum_i psi(O_i, theta) is close to zero when
#' evaluated at estimated parameters.
#' @param geex_list a list with at least \code{splitdt}, \code{eeFUN}, and \code{outer_eeargs} (which can be \code{NULL})
#' @param theta parameter estimates to use in evaluating the estimating equations.
#' @export
#'
#------------------------------------------------------------------------------#

check_GFUN <- function(geex_list, theta){
  psi_i <- create_psi(splitdt      = geex_list$splitdt,
                      eeFUN        = geex_list$eeFUN,
                      outer_eeargs = geex_list$outer_eeargs)
  GFUN <-  create_GFUN(psi_list = psi_i, inner_eeargs = geex_list$inner_eeargs)

  GFUN(theta)
}



#------------------------------------------------------------------------------#
#' Check value of estimating eqations
#'
#' Checks that the value of Gm = sum_i psi(O_i, theta) is close to zero when
#' evaluated at estimated parameters
#'
#------------------------------------------------------------------------------#

check_GFUN <- function(geex_list, theta, ...){
  psi_i <- create_psi(splitdt = geex_list$splitdt, eeFUN = geex_list$eeFUN, ...)
  GFUN <-  create_GFUN(psi_list = psi_i, ee_args = geex_list$ee_args)

  GFUN(theta)
}



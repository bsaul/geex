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
#' @examples
#' myee <- function(data) {
#'   function(theta) {
#'     c(
#'       data$Y1 - theta[1],
#'       (data$Y1 - theta[1])^2 - theta[2]
#'     )
#'   }
#' }
#'
#' mest <- m_estimate(
#'   estFUN = myee,
#'   data = geexex,
#'   root_control = setup_root_control(start = c(1, 1))
#' )
#'
#' f <- grab_GFUN(mest@basis)
#' # Should be close to zero
#' diagnose_roots(GFUN = f, theta = roots(mest))
#------------------------------------------------------------------------------#
diagnose_roots <- function(GFUN, theta) {
  GFUN(theta)
}

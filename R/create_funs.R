#------------------------------------------------------------------------------#
# create_** description:
# Functions that take in stuff and create a function or list of functions.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Creates list of psi functions
#'
#' Creates the estimating function (\eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)})
#' for each unit. That is, this function evaluates the outer function in \code{eeFUN}
#' for each independent unit and a returns the inner function in \code{eeFUN}.
#'
#' @param splitdt list of dataframes with data per unit
#' @param eeFUN the estimating equation function
#' @param approxFUN a function that approximates the inner function of \code{eeFUN}.
#' (EXPERIMENTAL).
#' @param approxFUN_control arguments passed to \code{approxFUN}
#' @param outer_eeargs a list of arguments passed to \code{eeFUN}
#' @return a list of functions, each function corresponding to a single unit
#' @export
#'
#------------------------------------------------------------------------------#

create_psi <- function(splitdt,
                       eeFUN,
                       approxFUN = NULL,
                       approxFUN_control = NULL,
                       outer_eeargs = NULL){
  out <- lapply(splitdt, function(data_i){
    do.call(eeFUN, args = append(list(data = data_i), outer_eeargs))
  })

  # if user specifies an approximation function, apply the function to each
  # evaluation of psi
  if(!is.null(approxFUN)){
    lapply(out, function(f){
      do.call(approxFUN, args = append(list(psi = f), approxFUN_control))
    }) -> out
  }

  out
}

#------------------------------------------------------------------------------#
#' Creates a function that sums over psi functions
#'
#' From a list of \eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)} for i = 1, ..., m,
#' creates \eqn{G_m = \sum_i \psi(O_i, \theta)}{G_m = \sum_i \psi(O_i, \theta)}. Here,
#' \eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)} is the *inner* part of an \code{eeFUN},
#' in that the data is fixed and \eqn{G_m}{G_m} is a function of \eqn{\theta)}{\theta}.
#'
#' @param psi_list list of psi functions
#' @param inner_eeargs list of arguments passed to psi
#' @export
#'
#------------------------------------------------------------------------------#

create_GFUN <- function(psi_list, inner_eeargs = NULL, weights = NULL){
  function(theta){
    psii <- lapply(psi_list, function(f) {
      do.call(f, args = append(list(theta = theta), inner_eeargs))
    })

    # If weights are provided, then multiply each psi function by its
    # respective weight
    if(is.null(weights)){
      psii_array <- simplify2array(psii)
    } else {
      psii_array <- simplify2array(Map(`*`, psii, weights))
    }
    # sum over unit-wise contributions to the estimating equations
    apply(check_array(psii_array), 1, sum)
  }
}

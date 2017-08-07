#------------------------------------------------------------------------------#
# create_** description:
# Functions that take in stuff and create a function or list of functions.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Creates list of psi functions
#'
#' Creates the estimating function (\eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)})
#' for each unit. That is, this function evaluates the outer function in
#' \code{estFUN} for each independent unit and a returns the inner function in
#' \code{estFUN}.
#'
#' @param .basis an object of class \code{\linkS4class{m_estimation_basis}}
#' @param .approx_control an object of class \code{\linkS4class{approx_control}}
#' or \code{NULL}
#' @return the \code{.basis} with the \code{.psiFUN_list} slot populated.
#' @export
#'
#------------------------------------------------------------------------------#

create_psiFUN_list <- function(.basis,
                               .approx_control){

  out <- lapply(.basis@.split_data, function(data_i){
    do.call(grab_estFUN(.basis),
            args = append(list(data = data_i), .basis@.outer_args))
  })

  # if user specifies an approximation function, apply the function to each
  # evaluation of psi

  # Use approx_control defaults if no options passed
  if(missing(.approx_control)){
    .approx_control <- new('approx_control')
  }
  approxFUN <- FUN(.approx_control)
  if(!(is.null(body(approxFUN)))){
    lapply(out, function(f){
      do.call(approxFUN, args = append(list(psi = f), options(.approx_control)))
    }) -> out
  }

  set_psiFUN_list(.basis) <- out
  .basis
}

#------------------------------------------------------------------------------#
#' Creates a function that sums over psi functions
#'
#' From a list of \eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)} for i = 1, ..., m,
#' creates \eqn{G_m = \sum_i \psi(O_i, \theta)}{G_m = \sum_i \psi(O_i, \theta)},
#' called \code{GFUN}. Here, \eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)} is the
#' *inner* part of an \code{estFUN}, in that the data is fixed and \eqn{G_m}{G_m}
#' is a function of \eqn{\theta)}{\theta}.
#'
#' @inheritParams create_psiFUN_list
#' @export
#'
#------------------------------------------------------------------------------#

create_GFUN <- function(.basis){
  psi_list <- get_psiFUN_list(.basis)
  function(theta){
    psii <- lapply(psi_list, function(psi) {
      do.call(psi, args = append(list(theta = theta), .basis@.inner_args))
    })

    # If weights are provided, then multiply each psi function by its
    # respective weight
    if(length(.basis@.weights) == 0){
      psii_array <- simplify2array(psii)
    } else {
      psii_array <- simplify2array(Map(`*`, psii, .weights))
    }
    # sum over unit-wise contributions to the estimating equations
    apply(check_array(psii_array), 1, sum)
  }
}

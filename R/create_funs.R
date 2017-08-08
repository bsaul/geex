#------------------------------------------------------------------------------#
# create_** description:
# Functions that take in stuff and create a function or list of functions.
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Creates an m_estimation_basis object
#'
#' @inheritParams m_estimate
#' @details Either \code{data} or \code{split_data} must be provided
#' @return a \code{\linkS4class{m_estimation_basis}}
#' @export
#------------------------------------------------------------------------------#

create_basis <- function(estFUN, data, units, outer_args, inner_args){
  new(Class="m_estimation_basis",
      .estFUN = estFUN,
      .data   = data,
      .units  =  if(!missing(units)) units else character(),
      .outer_args = if(!missing(outer_args)) outer_args else list(),
      .inner_args = if(!missing(inner_args)) inner_args else list() )
}

#------------------------------------------------------------------------------#
#' Creates list of psi functions
#'
#' Creates the estimating function (\eqn{\psi(O_i, \theta)}{\psi(O_i, \theta)})
#' for each unit. That is, this function evaluates the outer function in
#' \code{estFUN} for each independent unit and a returns the inner function in
#' \code{estFUN}.
#'
#' @param .basis an object of class \code{\linkS4class{m_estimation_basis}}
#' @return the \code{.basis} with the \code{.psiFUN_list} slot populated.
#' @export
#'
#------------------------------------------------------------------------------#
setGeneric("create_psiFUN_list", function(object, ...) standardGeneric("create_psiFUN_list"))
setMethod(
  f = "create_psiFUN_list",
  signature = "m_estimation_basis",
  definition = function(object){

    # Split data frame into data frames for each independent unit
    # if units are not specified, split into one per observation
    dt <- grab_basis_data(object)
    ut <- if(length(object@.units) == 0) 1:nrow(dt) else dt[[object@.units]]
    split_data <- split(x = dt, f = ut)

    # Apply estFUN to each unit's data
    out <- lapply(split_data, function(data_i){
      do.call(grab_estFUN(object),
              args = append(list(data = data_i), object@.outer_args))
    })

    # if user specifies an approximation function, apply the function to each
    # evaluation of psi

    appFUN  <- grab_approxFUN(object@.control)
    appopts <- grab_approx_options(object@.control)
    if(!(is.null(body(appFUN)))){
      lapply(out, function(f){
        do.call(appFUN, args = append(list(psi = f), appopts))
      }) -> out
    }

    set_psiFUN_list(object) <- out
    object
  }
)

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
  psi_list <- grab_psiFUN_list(.basis)
  function(theta){
    psii <- lapply(psi_list, function(psi) {
      do.call(psi, args = append(list(theta = theta), .basis@.inner_args))
    })

    compute_sum_of_list(psii, .basis@.weights)
  }
}

#------------------------------------------------------------------------------#
#' Grab something from an object
#'
#' @param from an object
#' @param what what to grab
#' @param ... additional arguments passed to \code{grab_**} function
#' @export
#------------------------------------------------------------------------------#

grab <- function(from, what, ...){
  switch(what,
         "response"         = grab_response(data = from, ...),
         "design_matrix"    = grab_design_matrix(data = from, ...),
         "response_formula" = grab_response_formula(model = from),
         "fixed_formula"    = grab_fixed_formula(model = from))
}

#------------------------------------------------------------------------------#
#' Grab the RHS formula from a model object
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @export
#------------------------------------------------------------------------------#
grab_fixed_formula <- function(model){
  stats::formula(model, fixed.only = TRUE)[-2]
}

#------------------------------------------------------------------------------#
#' Grab the LHS formula from a model object
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @export
#------------------------------------------------------------------------------#
grab_response_formula <- function(model){
  stats::formula(model)[-3]
}

#------------------------------------------------------------------------------#
#' Grab a matrix of fixed effects from a model object
#'
#' @param rhs_formula the right hand side of a model formula
#' @param data the data from which to extract the matrix
#' @export
#------------------------------------------------------------------------------#
grab_design_matrix <- function(data, rhs_formula){
  stats::model.matrix(object = rhs_formula, data = data)
}

#------------------------------------------------------------------------------#
#' Grab a vector of responses from a model object
#'
#' @param formula model formula
#' @param data data.frame from which to extract the vector of responses
#' @export
#------------------------------------------------------------------------------#

grab_response <- function(data, formula){
  stopifnot(class(formula) == 'formula')
  stats::model.response(stats::model.frame(formula = formula, data = data))
}

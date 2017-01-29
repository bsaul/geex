#------------------------------------------------------------------------------#
#' Check an array object
#'
#' Handles the case where there is a single estimating equation. This function
#' assumes that the object
#'
#' @param object the object to check whether it is an array
#' @return an array - either the orginal object or the given object converted
#' to an array
#------------------------------------------------------------------------------#

check_array <- function(object){
  if(is.array(object)){
    object
  } else if(is.numeric(object)){
    array(object, dim = c(1, 1, length(object)))
  } else if(is.matrix(object)){
    array(object, dim = c(1, 1, length(object)))
  } else {
    stop('Object is not an array, matrix, or numeric')
  }
}

#------------------------------------------------------------------------------#
#' Get the RHS formula from a model object
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @export
#------------------------------------------------------------------------------#
get_fixed_formula <- function(model){
  formula(model, fixed.only = TRUE)[-2]
}

#------------------------------------------------------------------------------#
#' Get the LHS formula from a model object
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @export
#------------------------------------------------------------------------------#
get_response_formula <- function(model){
  formula(model)[-3]
}

#------------------------------------------------------------------------------#
#' Get a matrix of fixed effects from a model object
#'
#' @param rhs_formula the right hand side of a model formula
#' @param data the data from which to extract the matrix
#' @export
#------------------------------------------------------------------------------#
get_design_matrix <- function(rhs_formula, data){
  model.matrix(rhs_formula, data)
}

#------------------------------------------------------------------------------#
#' Get a vector of responses from a model object
#'
#' @param formula model formula
#' @param data data.frame from which to extract the vector of responses
#' @export
#------------------------------------------------------------------------------#

get_response <- function(formula, data){
  model.response(model.frame(formula, data = data))
}

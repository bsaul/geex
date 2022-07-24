#------------------------------------------------------------------------------#
# grab_** description:
# Function that takes some object and "grabs" "what" "from" the object.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Grab something from an object
#'
#' @param from an object
#' @param what what to grab one of 'response', 'design_matrix', 'response_formula',
#' 'fixed_formula', 'eeFUN'
#' @param ... additional arguments passed to \code{grab_**} function
#' @seealso \code{\link{grab_response}}, \code{\link{grab_design_matrix}},
#' \code{\link{grab_response_formula}}, \code{\link{grab_fixed_formula}},
#' \code{\link{grab_design_levels}}
#' @return the value returns depends on the argument \code{what}.
#' @export
#------------------------------------------------------------------------------#

grab <- function(from, what, ...){
  switch(what,
         "response"         = grab_response(data = from, ...),
         "design_matrix"    = grab_design_matrix(data = from, ...),
         "design_levels"    = grab_design_levels(model = from),
         "response_formula" = grab_response_formula(model = from),
         "fixed_formula"    = grab_fixed_formula(model = from),
         'psiFUN'           = grab_psiFUN(object = from, ...),
         stop("'what' you want to grab() is not defined"))
}

#------------------------------------------------------------------------------#
#' Grab the RHS formula from a model object
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @return the right-hand side of a model's \code{\link[stats]{formula}} object
#' @export
#' @examples
#' fit <- lm(Sepal.Width ~ Petal.Width, data = iris)
#' grab_fixed_formula(fit)
#------------------------------------------------------------------------------#
grab_fixed_formula <- function(model){
  stats::formula(model, fixed.only = TRUE)[-2]
}

#------------------------------------------------------------------------------#
#' Grab the LHS formula from a model object
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @return the left-hand side of a model's \code{\link[stats]{formula}} object
#' @export
#' @examples
#' fit <- lm(Sepal.Width ~ Petal.Width, data = iris)
#' grab_response_formula(fit)
#------------------------------------------------------------------------------#
grab_response_formula <- function(model){
  stats::formula(model)[-3]
}

#------------------------------------------------------------------------------#
#' Grab a matrix of fixed effects from a model object
#'
#' @param rhs_formula the right hand side of a model formula
#' @param data the data from which to extract the matrix
#' @param ... Can be used to pass \code{xlev} to \code{\link[stats]{model.frame}}
#' @return a \code{\link[stats]{model.matrix}}
#' @export
#' @examples
#' # Create a "desigm" matrix for the first ten rows of iris data
#' fit <- lm(Sepal.Width ~ Petal.Width, data = iris)
#' grab_design_matrix(
#'   data = iris[1:10, ],
#'   grab_fixed_formula(fit))
#------------------------------------------------------------------------------#
grab_design_matrix <- function(data, rhs_formula, ...){
  stats::model.matrix(object = rhs_formula, data = data, ...)
}

#------------------------------------------------------------------------------#
#' Grab a list of the levels of factor variables in a model.
#'
#' Useful when splitting data later, used with \code{\link{grab_design_matrix}}
#' or especially when calling \code{\link{grab_psiFUN}} from within an eeFun.
#'
#' @param model a model object such as \code{lm}, \code{glm}, \code{merMod}
#' @export
#'
#' @return A named list of character vectors that provides the fentire set of
#' levels that each factor predictor in \code{model} will take on. This is
#' hopefully identical to what the \code{xlev} argument to
#' \code{link[stats]{model.frame}} desires. When \code{model} has no factors
#' as predictors, then an empty list is returned.
#'
#' @examples
#' \dontrun{
#'   geex::grab_design_matrix(
#'     data = data,
#'     rhs_formula = geex::grab_fixed_formula(model),
#'     xlev = geex::grab_design_levels(model)
#'   )
#'   ## Below is helpful within an eeFun.
#'   geex::grab_psiFUN(
#'     data = data,## Especially when this is a subset of the data
#'     rhs_formula = geex::grab_fixed_formula(model),
#'     xlev = geex::grab_design_levels(model)
#'   )
#' }
#------------------------------------------------------------------------------#
grab_design_levels <- function(model){

  full_model_frame <- stats::model.frame(model)
  data_classes <- attr(attr(full_model_frame, "terms"), "dataClasses")
  var_names <- names(data_classes)

  x_levels_list <- list()
  ii = 1
  for (var_num in 1:length(data_classes)){
    if (data_classes[var_num]=="factor") {
      these_levels <- levels( full_model_frame[,var_names[var_num] ])
      x_levels_list[[ii]] <- these_levels
      names(x_levels_list)[[ii]] <- var_names[var_num]
      ii <- ii+1
    }
  }

  x_levels_list
}

#------------------------------------------------------------------------------#
#' Grab a vector of responses from a model object
#'
#' @param formula model formula
#' @param data data.frame from which to extract the vector of responses
#' @return a \code{\link[stats]{model.response}}
#' @export
#' @examples
#' # Grab vector of responses for the first ten rows of iris data
#' fit <- lm(Sepal.Width ~ Petal.Width, data = iris)
#' grab_response(
#'   data = iris[1:10, ],
#'   formula(fit))
#------------------------------------------------------------------------------#

grab_response <- function(data, formula){
  stopifnot(class(formula) == 'formula')
  stats::model.response(stats::model.frame(formula = formula, data = data))
}

#------------------------------------------------------------------------------#
# define S4 Classes #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' estimation_function S4 class
#' @slot .estFUN the estimating function.
#' @slot .outer_args a named \code{list} of arguments passed to the outer
#' function of \code{estFUN}. Should *not* include the \code{data} argument.
#' @slot .inner_args a named \code{list} of arguments passed to the inner
#' function of \code{estFUN}. Should *not* include the \code{theta} argument.
#'
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = 'estimating_function',
  slots = c(.estFUN = "function",
            .outer_args = "list",
            .inner_args = "list"),
  validity = function(object){

    outer_estFUN_args <- formalArgs(grab_estFUN(object))
    inner_estFUN_args <- formalArgs(eval(body(grab_estFUN(object))))

    outer_args_names <- names(object@.outer_args)
    inner_args_names <- names(object@.inner_args)

    # check that first argument of outer eeFUN is data
    if(outer_estFUN_args[1] != 'data'){
      "First argument of the outer estFUN must be 'data'"
    }

    # check that first argument of inner eeFUN is theta
    else if(inner_estFUN_args[1] != 'theta'){
      "First argument of the inner estFUN must be 'theta'"
    }

    else if(length(outer_args_names) > 0){
      if(is.null(outer_args_names)){
        "outer_args must be a *named* list with names matching arguments in the outer estFUN"
      }

      else if("data" %in% outer_args_names){
        "the argument data should not be included in the outer_args"
      }

      # check that outer estFUN uses outer_args
      else if(any(!(names(object@outer_args) %in% outer_estFUN_args))){
        "outer_args contains elements not used in the outer estFUN"
      }
    }

    else if(length(inner_args_names) > 0){
      if(is.null(inner_args_names)){
        "inner_args must be a *named* list with names matching arguments in the inner estFUN"
      }

      else if("theta" %in% inner_args_names){
        "the argument data should not be included in the inner_args"
      }

      # check that inner estFUN uses inner_args,
      else if(any(!(names(object@inner_args) %in% inner_estFUN_args))){
        "inner_args contains elements not used in the inner estFUN"
      }
    }

    else TRUE
  }
)

#------------------------------------------------------------------------------#
#' m_estimation_basis S4 class
#'
#' @slot .data the analysis data.frame
#' @slot .units an (optional) character string identifying the variable in
#' @slot .split_data
#' \code{data} which
#' splits the data into indepedent units
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = "m_estimation_basis",
  slots = c(.data  = "data.frame",
            .units = "character",
            .split_data = "list",
            .weights = "numeric"),
  contains = "estimating_function",
  validity = function(object){

    if(length(object@.units) > 1){
      "units should be a character string identifying the name of the variable in object@.data"
    }

    else if(length(names(object@.data)) > 0 & length(object@.units) > 0){
      if(length(object@.units) == 1 & !(object@.units %in% names(object@.data))){
        paste(object@.units, " is not a variable in the data")
      }
    }

    else TRUE
})

#------------------------------------------------------------------------------#
#' Initialize a m_estimation_basis object
#'
#' @export
#------------------------------------------------------------------------------#

setMethod("initialize", "m_estimation_basis", function(.Object, ...){
  .Object <- callNextMethod()

  # TODO: This set up allows for the case where the split_data could be set to
  # a different split than that defined by .units. This shouldn't be the case

  # TODO: an m_estimation_basis object will carry around both the .data
  # and the .split_data, but this seems redundant

  if(length(.Object@.split_data) == 0){
    dt <- grab_basis_data(.Object)
    ut <- if(length(.Object@.units) == 0) 1:nrow(dt) else dt[[.Object@.units]]
    # Split data frame into data frames for each independent unit
    split(x = dt,
          # if units are not specified, split into one per observation
          f = ut) ->
      .Object@.split_data
  }
  .Object
})

#------------------------------------------------------------------------------#
#' root_control S4 class
#'
#' @slot .FUN a root finding function
#' @slot .options a list of options passed to \code{.FUN}
#' @slot .object_name a character string identifying the object containing the
#' roots in the output of \code{.FUN}
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "root_control",
  slots = c(.FUN = 'function',
            .options = 'list',
            .object_name = 'character'),
  validity = function(object){

    FUN_arg_names <- formalArgs(object@.FUN)

    if(FUN_arg_names[1] != 'f'){
      "The first argument of FUN must be 'f', as in multiroot() or uniroot()"
    }

    else if(any(!(names(object@.options) %in% FUN_arg_names))){
      "The root_control options contains arguments not used in the FUN"
    }

    else TRUE
  },
  prototype = prototype(
    .FUN = rootSolve::multiroot,
    .object_name = 'root'
  )
)

#------------------------------------------------------------------------------#
#' grab_estFUN generic
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_estFUN", function(object, ...) standardGeneric("grab_estFUN"))
setMethod("grab_estFUN", "estimating_function", function(object) object@.estFUN)

#------------------------------------------------------------------------------#
#' grab_basis_data generic
#'
#' A function that grabs the data from an m_estimation_basis object
#'
#------------------------------------------------------------------------------#

setGeneric("grab_basis_data", function(object, ...) standardGeneric("grab_basis_data"))
setMethod("grab_basis_data", "m_estimation_basis", function(object) object@.data)


#------------------------------------------------------------------------------#
#' Show the m_estimation_basis
#'
#' @export
#------------------------------------------------------------------------------#

setMethod(
  "show",
  signature = "m_estimation_basis",
  definition = function(object) {
    cat("An object of class ", class(object), "\n", sep = "")
    cat("psi: \n")
    print(body(object@.estFUN))
    cat("Data:\n")
    print(head(object@.data))
    cat("Units: ", object@.units)

    invisible(NULL)
  })

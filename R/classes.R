#------------------------------------------------------------------------------#
# define S4 Classes and methods used within geex #
#------------------------------------------------------------------------------#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## estimating_function class ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#------------------------------------------------------------------------------#
#' estimating_function S4 class
#'
#' @slot .estFUN the estimating function.
#' @slot .outer_args a named \code{list} of arguments passed to the outer
#' function of \code{.estFUN}. Should *not* include the \code{data} argument.
#' @slot .inner_args a named \code{list} of arguments passed to the inner
#' function of \code{.estFUN}. Should *not* include the \code{theta} argument.
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
    # TODO: Checking the inner_estFUN_args doesn't work because body() grabs
    # more than just the  anynomous function returned by estFUN (duh).
    # inner_estFUN_args <- formalArgs(eval(body(grab_estFUN(object))))

    outer_args_names <- names(object@.outer_args)
    inner_args_names <- names(object@.inner_args)

    # check that first argument of outer eeFUN is data
    if(outer_estFUN_args[1] != 'data'){
      "First argument of the outer estFUN must be 'data'"
    }

    # check that first argument of inner eeFUN is theta
    # See TODO above
    # else if(inner_estFUN_args[1] != 'theta'){
    #   "First argument of the inner estFUN must be 'theta'"
    # }

    # Check outer_args
    else if(length(outer_args_names) > 0){
      if(is.null(outer_args_names)){
        "outer_args must be a *named* list with names matching arguments in the outer estFUN"
      }

      else if("data" %in% outer_args_names){
        "the argument data should not be included in the outer_args"
      }

      # check that outer estFUN uses outer_args
      else if(any(!(outer_args_names %in% outer_estFUN_args))){
        "outer_args contains elements not used in the outer estFUN"
      }
    }

    # Check inner_args
    else if(length(inner_args_names) > 0){
      if(is.null(inner_args_names)){
        "inner_args must be a *named* list with names matching arguments in the inner estFUN"
      }

      else if("theta" %in% inner_args_names){
        "the argument data should not be included in the inner_args"
      }

      # check that inner estFUN uses inner_args,
      # See TODO above
      # else if(any(!(inner_args_names %in% inner_estFUN_args))){
      #   "inner_args contains elements not used in the inner estFUN"
      # }
    }

    else TRUE
  }
)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## sandwich_components class and methods ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#------------------------------------------------------------------------------#
#' sandwich_components S4 class
#'
#' @slot .A the "bread" matrix
#' @slot .A_i a list of "bread" matrices per unit
#' @slot .B the "meat" matrix
#' @slot .B_i a list of "meat" matrices per unit
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "sandwich_components",
  slots = c(.A   = 'matrix',
            .A_i = 'list',
            .B   = 'matrix',
            .B_i = 'list'),
  validity = function(object){
    A_dims <- dim(object@.A)
    B_dims <- dim(object@.B)
    if(A_dims[1] != A_dims[2]){
      "The .A (bread) matrix must be square"
    }
    else if(B_dims[1] != B_dims[2]){
      "The .B (meat) matrix must be square"
    }
    else TRUE
  }
)

#------------------------------------------------------------------------------#
#' Shows the sandwich_components S4 class
#'
#' @param object a \code{\linkS4class{sandwich_components}} object
#' @export
#------------------------------------------------------------------------------#

setMethod(
  "show",
  signature = "sandwich_components",
  definition = function(object) {
    cat("An object of class ", class(object), "\n", sep = "")
    cat("A (bread matrix):\n")
    print(object@.A)
    cat("B (meat matrix):\n")
    print(object@.B)

    invisible(NULL)
  })
#------------------------------------------------------------------------------#
#' Grabs the .A (bread matrix) slot
#'
#' @param object a \code{\linkS4class{sandwich_components}} object
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_bread",function(object){standardGeneric ("grab_bread")})
setMethod(
  f = "grab_bread",
  signature = "sandwich_components",
  function(object){
    return(object@.A)
  })

#------------------------------------------------------------------------------#
#' Gets the .A_i (list of bread matrices) slot
#'
#' @param object a \code{\linkS4class{sandwich_components}} object
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_bread_list",function(object){standardGeneric ("grab_bread_list")})
setMethod(
  f = "grab_bread_list",
  signature = "sandwich_components",
  function(object){
    return(object@.A_i)
  })

#------------------------------------------------------------------------------#
#' Gets the .B (meat matrix) slot
#'
#' @param object a \code{\linkS4class{sandwich_components}} object
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_meat",function(object){standardGeneric ("grab_meat")})
setMethod(
  f = "grab_meat",
  signature = "sandwich_components",
  function(object){
    return(object@.B)
  })

#------------------------------------------------------------------------------#
#' Gets the .B_i (list of bread matrices) slot
#'
#' @param object a \code{\linkS4class{sandwich_components}} object
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_meat_list",function(object){standardGeneric ("grab_meat_list")})
setMethod(
  f = "grab_meat_list",
  signature = "sandwich_components",
  function(object){
    return(object@.B_i)
  })

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## control classes ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#------------------------------------------------------------------------------#
#' basic_control S4 class
#'
#' A general class for defining a \code{function}, and the options passed to the
#' function
#'
#' @slot .FUN a function
#' @slot .options a list of options passed to \code{.FUN}
#' @seealso root_control deriv_control approx_control
#'
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = "basic_control",
  slots = c(.FUN = 'function',
            .options = 'list')
)

#------------------------------------------------------------------------------#
#' root_control S4 class
#'
#' @slot .FUN a root finding function whose first argument must be named \code{f}.
#' @slot .options a list of options passed to \code{.FUN}.
#' @slot .object_name a character string identifying the object containing the
#' roots in the output of \code{.FUN}.
#'
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = "root_control",
  slots = c(.object_name = 'character'),
  contains = "basic_control",
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
#' Setup a root_control object
#'
#' @param FUN a root finding function whose first argument must be named \code{f}.
#' @param roots_name a character string identifying the object containing the
#' roots in the output of \code{FUN}.
#' @param ... arguments passed to \code{FUN}
#' @return a \code{\linkS4class{root_control}} object
#' @export
#------------------------------------------------------------------------------#

setup_root_solver <- function(FUN, roots_name, ...){
  dots <- list(...)
  hold <- call('new')
  hold[['Class']] <- "root_control"
  if(!missing(FUN)){
    hold[[".FUN"]] <- FUN
  }
  if(length(dots) > 0){
    hold[[".options"]] <- dots
  }
  if(!missing(roots_name)){
    hold[[".object_name"]] <- roots_name
  }

  eval(hold)
}

#------------------------------------------------------------------------------#
#' deriv_control S4 class
#'
#' @slot .FUN a function which computes a numerical derivation. This functions
#' first argument must the function on which the derivative is being compute.
#' Defaults to \code{\link[numDeriv]{jacobian}}.
#' @slot .options a list of options passed to \code{.FUN}. Defaults to
#' \code{list(method = 'Richardson')}
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "deriv_control",
  contains = "basic_control",
  prototype = prototype(
    .FUN = numDeriv::jacobian,
    .options = list(method = 'Richardson')
  )
)

#------------------------------------------------------------------------------#
#' approx_control S4 class
#'
#' EXPERIMENTAL. See example 7 in \code{vignette("01_additional_examples", package = "geex")}
#' for usage.
#'
#' @slot .FUN a function which approximates an \code{estFUN}.
#' @slot .options a list of options passed to \code{.FUN}.
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "approx_control",
  contains = "basic_control"
)

#------------------------------------------------------------------------------#
#' correct_control S4 class
#'
#' @slot .FUN a function which "corrects" a \code{\linkS4class{sandwich_components}}
#' object. Usually a small-sample correction
#' @slot .options a list of options passed to \code{.FUN}.
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "correct_control",
  contains = "basic_control",
  validity = function(object){
    args_names <- formalArgs(grab_FUN(object))
    option_names <- names(grab_options(object))
    if(!('components' == args_names[1])){
      "'components' must be the first argument of a correction function"
    } else if(any(!(option_names %in% args_names[-1]))) {
      "All correction options must be an argument for the correction function"
    } else
      TRUE
  }
)

#------------------------------------------------------------------------------#
#' Extracts the \code{.options} slot from a \code{\linkS4class{basic_control}} object
#'
#' @param object a \code{\linkS4class{basic_control}} object
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_options", function(object, ...) standardGeneric("grab_options"))
setMethod("grab_options", "basic_control", function(object) object@.options)

#------------------------------------------------------------------------------#
#' Extracts the \code{.FUN} slot from a \code{\linkS4class{basic_control}} object
#'
#' @param object a \code{\linkS4class{basic_control}} object
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_FUN", function(object, ...) standardGeneric("grab_FUN"))
setMethod("grab_FUN", "basic_control", function(object) object@.FUN)

#------------------------------------------------------------------------------#
#' geex_control S4 class
#'
#' An object which control all the \code{\linkS4class{basic_control}} objects
#' necessary to perform M-estimation
#'
#' @slot .approx an \code{\linkS4class{approx_control}} object
#' @slot .root  a \code{\linkS4class{root_control}} object
#' @slot .deriv  a \code{\linkS4class{deriv_control}} object
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "geex_control",
  slots = c(.approx = 'approx_control',
            .root   = 'root_control',
            .deriv  = 'deriv_control')
)

#------------------------------------------------------------------------------#
#' Extracts the \code{.approx@.FUN} slot from a \code{\linkS4class{geex_control}} object
#'
#' @param object a \code{\linkS4class{geex_control}} object
#' @param slot name of the slot from which to grab the function. One of "deriv",
#' "approx", or "root"
#' @export
#------------------------------------------------------------------------------#
setMethod("grab_FUN", "geex_control", function(object, slot) {
  switch(slot,
         "approx" = grab_FUN(object@.approx),
         "root"   = grab_FUN(object@.root),
         "deriv"  = grab_FUN(object@.deriv))
})

#------------------------------------------------------------------------------#
#' Extracts the \code{.approx@.options} slot from a \code{\linkS4class{geex_control}} object
#'
#' @param object a \code{\linkS4class{geex_control}} object
#' @param slot name of the slot from which to grab the function. One of "deriv",
#' "approx", or "root"
#'
#' @export
#------------------------------------------------------------------------------#
setMethod("grab_options", "geex_control", function(object, slot) {
  switch(slot,
         "approx" = grab_options(object@.approx),
         "root"   = grab_options(object@.root),
         "deriv"  = grab_options(object@.deriv))
})

#------------------------------------------------------------------------------#
#' Sets the slots from for \code{\linkS4class{geex_control}} object
#'
#' @param object a \code{\linkS4class{geex_control}} object
#' @param slot name of the slot from which to grab the function. One of "deriv",
#' "approx", or "root"
#' @param value a \code{\linkS4class{root_control}}, \code{\linkS4class{apprix_control}},
#' or \code{\linkS4class{deriv_control}} object
#'
#' @export
#------------------------------------------------------------------------------#
setGeneric("set_control<-", function(object, slot, value){
  standardGeneric("set_control<-")})

setReplaceMethod(
  f = "set_control",
  signature="geex_control",
  definition = function(object, slot, value){
    if(slot == 'deriv'){
      object@.deriv <- value
    }
    if(slot == 'root'){
      object@.root <- value
    }
    if(slot == 'approx'){
      object@.approx <- value
    }
    validObject(object)
    return (object)
  }
)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## m_estimation class and methods ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#------------------------------------------------------------------------------#
#' m_estimation_basis S4 class
#'
#' @slot .data the analysis data.frame
#' @slot .units an (optional) character string identifying the variable in
#' \code{.data} which splits the data into indepedent units
#' @slot .weights a numeric vector of weights used in weighting the estimating
#' functions
#' @slot .psiFUN_list a list of \code{psiFUN}s created by \code{\link{create_psiFUN_list}}
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = "m_estimation_basis",
  slots = c(.data        = "data.frame",
            .units       = "character",
            .weights     = "numeric",
            .psiFUN_list = "list",
            .GFUN        = "function",
            .control     = "geex_control"),
  contains = "estimating_function",
  validity = function(object){

    if(length(object@.units) > 1){
      "units should be a character string of the name of single variable in data"
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
#'# @export
#------------------------------------------------------------------------------#

setMethod("initialize", "m_estimation_basis", function(.Object, ...){
  .Object <- callNextMethod()

  .Object <- create_psiFUN_list(.Object)
  .Object <- create_GFUN(.Object)
  .Object
})

#------------------------------------------------------------------------------#
#' Sets the .psi_list slot in a m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#' @param value a \code{list} of psiFUNs created by \code{\link{create_psiFUN_list}}
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("set_psiFUN_list<-", function(object, value){
  standardGeneric("set_psiFUN_list<-")})

setReplaceMethod(
  f = "set_psiFUN_list",
  signature="m_estimation_basis",
  definition = function(object, value){
    object@.psiFUN_list <- value
    validObject(object)
    return (object)
  }
)

#------------------------------------------------------------------------------#
#' Gets the .psi_list slot in a m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_psiFUN_list",function(object){standardGeneric ("grab_psiFUN_list")})
setMethod(
  f = "grab_psiFUN_list",
  signature = "m_estimation_basis",
  function(object){
    return(object@.psiFUN_list)
  })

#------------------------------------------------------------------------------#
#' Sets the .psi_list slot in a m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#' @param value a \code{function} created \code{\link{create_GFUN}}
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("set_GFUN<-", function(object, value){
  standardGeneric("set_GFUN<-")})

setReplaceMethod(
  f = "set_GFUN",
  signature="m_estimation_basis",
  definition = function(object, value){
    object@.GFUN <- value
    validObject(object)
    return (object)
  }
)

#------------------------------------------------------------------------------#
#' Gets the .psi_list slot in a m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#'
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_GFUN",function(object){standardGeneric ("grab_GFUN")})
setMethod(
  f = "grab_GFUN",
  signature = "m_estimation_basis",
  function(object){
    return(object@.GFUN)
  })

#------------------------------------------------------------------------------#
#' Shows the m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
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

#------------------------------------------------------------------------------#
#' grab_basis_data generic
#'
#' Grabs the \code{.data} from an \code{\linkS4class{m_estimation_basis}} object
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#------------------------------------------------------------------------------#

setGeneric("grab_basis_data", function(object, ...) standardGeneric("grab_basis_data"))
setMethod("grab_basis_data", "m_estimation_basis", function(object) object@.data)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## geex class and methods ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#------------------------------------------------------------------------------#
#' geex S4 class
#'
#' @slot basis a \code{\linkS4class{m_estimation_basis}} object
#' @slot rootFUN_results the results of call to the root finding algorithm function
#' @slot sandwich_components a \code{\linkS4class{sandwich_components}} object
#' @slot GFUN the \code{function} of which the roots are computed.
#' @slot corrections a \code{list} of correction performed on \code{sandwich_components}
#' @slot estimates a \code{numeric} vector of parameter estimates
#' @slot vcov the empirical sandwich variance \code{matrix}
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = "geex",
  slots = c(basis           = "m_estimation_basis",
            rootFUN_results = "ANY",
            sandwich_components = "sandwich_components",
            GFUN            = "function",
            corrections     = "list",
            estimates       = "numeric",
            vcov            = "matrix"))


#------------------------------------------------------------------------------#
#' Shows the geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @export
#------------------------------------------------------------------------------#

setMethod(
  "show",
  signature = "geex",
  definition = function(object) {
    cat("M-estimation results based on the estimating function:\n", sep = "")
    print(body(object@basis@.estFUN))
    cat("\nEstimates:\n")
    print(object@estimates)
    cat("\nCovariance:\n")
    print(object@vcov)
    if(length(object@corrections) > 0){
      cat("The results include", length(object@corrections),
          "covariance corrections")
    }

    invisible(NULL)
  })


#------------------------------------------------------------------------------#
#' Gets the variance-covariance matrix from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#'
#' @export
#------------------------------------------------------------------------------#

setMethod(
  "vcov",
  signature = "geex",
  definition = function(object) {
    object@vcov
  })

#------------------------------------------------------------------------------#
#' Gets the parameter estimates from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @export
#------------------------------------------------------------------------------#

setMethod(
  "coef",
  signature = "geex",
  definition = function(object) {
    object@estimates
  })

#------------------------------------------------------------------------------#
#' Gets the parameter estimates matrix from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @export
#------------------------------------------------------------------------------#

setGeneric("roots", function(object, ...) standardGeneric("roots"))
setMethod(
  "roots",
  signature = "geex",
  definition = function(object) {
    object@estimates
  })

#------------------------------------------------------------------------------#
#' Gets the corrections from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @export
#------------------------------------------------------------------------------#

setGeneric("get_corrections", function(object, ...) standardGeneric("get_corrections"))
setMethod(
  "get_corrections",
  signature = "geex",
  definition = function(object) {
    if(length(object@corrections) == 0){
      "No corrections were performed on this object"
    } else {
      object@corrections
    }
  })

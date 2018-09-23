#------------------------------------------------------------------------------#
# define S4 Classes and methods used on these classes within geex #
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

#------------------------------------------------------------------------------#
#' Grab estimating functions from a model object
#'
#' @param object a \code{\linkS4class{estimating_function}} object
#' @docType methods
#' @rdname grab_estFUN-methods
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_estFUN",function(object){ standardGeneric("grab_estFUN") })

#' @rdname grab_estFUN-methods
#' @aliases grab_estFUN,estimating_function,estimating_function-method

setMethod("grab_estFUN", "estimating_function", function(object) object@.estFUN)

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
#' @slot .ee_i a list of observed estimating function values per unit
#'
#' @export
#------------------------------------------------------------------------------#
setClass(
  Class = "sandwich_components",
  slots = c(.A    = 'matrix',
            .A_i  = 'list',
            .B    = 'matrix',
            .B_i  = 'list',
            .ee_i = 'list'),
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
#' Show (print) the S4 geex classes
#'
#' @param object the object to print
#' @docType methods
#' @rdname show-methods
#' @export

setGeneric("show",function(object){standardGeneric("show")})

#' Shows the sandwich_components S4 class
#'
#' \code{\linkS4class{m_estimation_basis}}, or \code{\linkS4class{geex}} object
#' @rdname show-methods
#' @aliases show,sandwich_components,sandwich_components-method
#' @export

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
#' @docType methods
#' @rdname grab_bread-methods
#' @export
#' @examples
#' myee <- function(data){
#'  function(theta){
#'    c(data$Y1 - theta[1],
#'    (data$Y1 - theta[1])^2 - theta[2])
#'   }
#' }
#'
#' results <- m_estimate(
#'    estFUN = myee,
#'    data = geexex,
#'    root_control = setup_root_control(start = c(1,1)))
#'
#' grab_bread(results@sandwich_components)
#------------------------------------------------------------------------------#

setGeneric("grab_bread",function(object){standardGeneric("grab_bread")})

#' @rdname grab_bread-methods
#' @aliases grab_bread,sandwich_components,sandwich_components-method

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
#' @docType methods
#' @rdname grab_bread_list-methods
#' @export
#' @examples
#' myee <- function(data){
#'  function(theta){
#'    c(data$Y1 - theta[1],
#'    (data$Y1 - theta[1])^2 - theta[2])
#'   }
#' }
#'
#' results <- m_estimate(
#'    estFUN = myee,
#'    data = geexex,
#'    root_control = setup_root_control(start = c(1,1)))
#'
#' head(grab_bread_list(results@sandwich_components))
#------------------------------------------------------------------------------#

setGeneric("grab_bread_list",function(object){standardGeneric ("grab_bread_list")})

#' @rdname grab_bread_list-methods
#' @aliases grab_bread_list,sandwich_components,sandwich_components-method

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
#' @docType methods
#' @rdname grab_meat-methods
#' @export
#' @examples
#' myee <- function(data){
#'  function(theta){
#'    c(data$Y1 - theta[1],
#'    (data$Y1 - theta[1])^2 - theta[2])
#'   }
#' }
#'
#' results <- m_estimate(
#'    estFUN = myee,
#'    data = geexex,
#'    root_control = setup_root_control(start = c(1,1)))
#'
#' grab_meat_list(results@sandwich_components)
#------------------------------------------------------------------------------#

setGeneric("grab_meat",function(object){standardGeneric ("grab_meat")})

#' @rdname grab_meat-methods
#' @aliases grab_meat,sandwich_components,sandwich_components-method

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
#' @docType methods
#' @rdname grab_meat_list-methods
#' @export
#' @examples
#' myee <- function(data){
#'  function(theta){
#'    c(data$Y1 - theta[1],
#'    (data$Y1 - theta[1])^2 - theta[2])
#'   }
#' }
#'
#' results <- m_estimate(
#'    estFUN = myee,
#'    data = geexex,
#'    root_control = setup_root_control(start = c(1,1)))
#'
#' head(grab_meat_list(results@sandwich_components))
#------------------------------------------------------------------------------#

setGeneric("grab_meat_list",function(object){standardGeneric ("grab_meat_list")})

#' @rdname grab_meat_list-methods
#' @aliases grab_meat_list,sandwich_components,sandwich_components-method

setMethod(
  f = "grab_meat_list",
  signature = "sandwich_components",
  function(object){
    return(object@.B_i)
  })


#------------------------------------------------------------------------------#
#' Gets the .ee_i (observed estimating function) slot
#'
#' @param object a \code{\linkS4class{sandwich_components}} object
#' @docType methods
#' @rdname grab_ee-methods
#' @export
#' @examples
#' myee <- function(data){
#'  function(theta){
#'    c(data$Y1 - theta[1],
#'    (data$Y1 - theta[1])^2 - theta[2])
#'   }
#' }
#'
#' results <- m_estimate(
#'    estFUN = myee,
#'    data = geexex,
#'    root_control = setup_root_control(start = c(1,1)))
#'
#' grab_ee_list(results@sandwich_components)
#------------------------------------------------------------------------------#


setGeneric("grab_ee_list",function(object){standardGeneric ("grab_ee_list")})

#' @rdname grab_meat_list-methods
#' @aliases grab_meat_list,sandwich_components,sandwich_components-method

setMethod(
  f = "grab_ee_list",
  signature = "sandwich_components",
  function(object){
    return(object@.ee_i)
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
#' @seealso \code{\link{root_control-class}}, \code{\link{deriv_control-class}}
#' \code{\link{approx_control-class}}
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
#' Setup a basic_control object
#'
#' @param type one of \code{c("deriv", "approx", "root")}
#' @param FUN a function
#' @param ... arguments passed to \code{FUN}
#' @return a \code{\linkS4class{basic_control}} object
#' @seealso \code{\link{setup_root_control}}, \code{\link{setup_deriv_control}},
#' \code{\link{setup_approx_control}}
#' @export
#------------------------------------------------------------------------------#

setup_control <- function(type, FUN, ...){
  if(!(type %in% c("deriv", "approx", "root"))){
    stop("type must be one of deriv, approx, root")
  }
  dots <- list(...)
  hold <- call('new')
  hold[['Class']] <- paste0(type, "_control")
  if(!missing(FUN)){
    hold[[".FUN"]] <- FUN
  }
  if(length(dots) > 0){
    hold[[".options"]] <- dots
  }

  eval(hold)
}

#------------------------------------------------------------------------------#
#' Setup a root_control object
#'
#' @param roots_name a character string identifying the object containing the
#' @inheritParams setup_control
#' @return a \code{\linkS4class{root_control}} object
#' @examples
#' # Setup the default
#' setup_root_control(start = c(3, 5, 6))
#'
#' # Also setup the default
#' setup_root_control(FUN = rootSolve::multiroot,
#'                    start = c(3, 5, 6))
#'
#' # Or use uniroot()
#' setup_root_control(FUN = stats::uniroot,
#'                    interval = c(0, 1))
#' @export
#------------------------------------------------------------------------------#

setup_root_control <- function(FUN, roots_name, ...){
  hold <- setup_control(
    type = "root",
    FUN  = FUN,
    ...  = ...)

  if(!missing(roots_name)){
    hold@.object_name <- roots_name
  }

  hold
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
#' Setup a deriv_control object
#'
#' @inheritParams setup_control
#' @return a \code{\linkS4class{deriv_control}} object
#' @export
#' @examples
#' setup_deriv_control() # default
#' setup_deriv_control(method = "simple") # will speed up computations
#------------------------------------------------------------------------------#

setup_deriv_control <- function(FUN, ...){
  setup_control(type = "deriv", FUN  = FUN, ...  = ...)
}

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
#' Setup an approx_control object
#'
#' @inheritParams setup_control
#' @return a \code{\linkS4class{approx_control}} object
#' @export
#' @examples
#' # For usage, see example 7 in
#' \dontrun{vignette("01_additional_examples", package = "geex")}
#------------------------------------------------------------------------------#

setup_approx_control <- function(FUN, ...){
  setup_control(type = "approx", FUN  = FUN, ... = ...)
}

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
# Extracts the \code{.options} slot from a \code{\linkS4class{basic_control}} object
#
# @param object a \code{\linkS4class{basic_control}} object
# @param ... additional arguments passed to other methods
# @docType methods
# @rdname grab_options-methods
#------------------------------------------------------------------------------#

setGeneric("grab_options", function(object, ...) standardGeneric("grab_options"))

# @rdname grab_options-methods
# @aliases grab_options,basic_control,basic_control-method

setMethod("grab_options", "basic_control", function(object) object@.options)

#------------------------------------------------------------------------------#
# Extracts the \code{.FUN} slot from a \code{\linkS4class{basic_control}} object
#
# @param object a \code{\linkS4class{basic_control}} object
# @param ... additional arguments passed to other methods
# @docType methods
# @rdname grab_FUN-methods
#------------------------------------------------------------------------------#

setGeneric("grab_FUN", function(object, ...) standardGeneric("grab_FUN"))

# @rdname grab_FUN-methods
# @aliases grab_FUN,basic_control,basic_control-method

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
# Extracts functions from a \code{\linkS4class{geex_control}} object
#
# @param slot name of the slot from which to grab the function. One of "deriv",
# "approx", or "root"
#
# @rdname grab_FUN-methods
# @aliases grab_FUN,geex_control,geex_control-method
#
# @export
#------------------------------------------------------------------------------#
setMethod("grab_FUN", "geex_control", function(object, slot) {
  switch(slot,
         "approx" = grab_FUN(object@.approx),
         "root"   = grab_FUN(object@.root),
         "deriv"  = grab_FUN(object@.deriv))
})

#------------------------------------------------------------------------------#
# Extracts options from a \code{\linkS4class{geex_control}} object
#
# @param slot name of the slot from which to grab the function. One of "deriv",
# "approx", or "root"
# @rdname grab_options-methods
# @aliases grab_options,geex_control,geex_control-method
#
# @export
#------------------------------------------------------------------------------#
setMethod("grab_options", "geex_control", function(object, slot) {
  switch(slot,
         "approx" = grab_options(object@.approx),
         "root"   = grab_options(object@.root),
         "deriv"  = grab_options(object@.deriv))
})

#------------------------------------------------------------------------------#
# Sets the slots from for \code{\linkS4class{geex_control}} object
#
# @param object a \code{\linkS4class{geex_control}} object
# @param slot name of the slot from which to grab the function. One of "deriv",
# "approx", or "root"
# @param value a \code{\linkS4class{root_control}}, \code{\linkS4class{approx_control}},
# or \code{\linkS4class{deriv_control}} object
#
# @export
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
    methods::validObject(object)
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
#' @slot .GFUN a function created by \code{\link{create_GFUN}}
#' @slot .control a \code{\linkS4class{geex_control}} object
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
# Initialize a m_estimation_basis object
#
# @export
#------------------------------------------------------------------------------#

setMethod("initialize", "m_estimation_basis", function(.Object, ...){
  .Object <- methods::callNextMethod()

  .Object <- create_psiFUN_list(.Object)
  .Object <- create_GFUN(.Object)
  .Object
})

#------------------------------------------------------------------------------#
# Sets the .psi_list slot in a m_estimation_basis
#
# @param object a \code{\linkS4class{m_estimation_basis}} object
# @param value a \code{list} of psiFUNs created by \code{\link{create_psiFUN_list}}
# @export
#------------------------------------------------------------------------------#

setGeneric("set_psiFUN_list<-", function(object, value){
  standardGeneric("set_psiFUN_list<-")})

setReplaceMethod(
  f = "set_psiFUN_list",
  signature="m_estimation_basis",
  definition = function(object, value){
    object@.psiFUN_list <- value
    methods::validObject(object)
    return (object)
  }
)

#------------------------------------------------------------------------------#
#' Gets the .psi_list slot in a m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#' @docType methods
#' @rdname grab_psiFUN_list-methods
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_psiFUN_list",function(object){standardGeneric("grab_psiFUN_list")})

#' @rdname grab_psiFUN_list-methods
#' @aliases grab_psiFUN_list,m_estimation_basis,m_estimation_basis-method

setMethod(
  f = "grab_psiFUN_list",
  signature = "m_estimation_basis",
  function(object){
    return(object@.psiFUN_list)
  })

#------------------------------------------------------------------------------#
# Sets the .psi_list slot in a m_estimation_basis
#
# @param object a \code{\linkS4class{m_estimation_basis}} object
# @param value a \code{function} created \code{\link{create_GFUN}}
#
# @export
#------------------------------------------------------------------------------#

setGeneric("set_GFUN<-", function(object, value){
  standardGeneric("set_GFUN<-")})

setReplaceMethod(
  f = "set_GFUN",
  signature="m_estimation_basis",
  definition = function(object, value){
    object@.GFUN <- value
    methods::validObject(object)
    return (object)
  }
)

#------------------------------------------------------------------------------#
#' Gets the .psi_list slot in a m_estimation_basis
#'
#' @param object a \code{\linkS4class{m_estimation_basis}} object
#' @docType methods
#' @rdname grab_GFUN-methods
#' @export
#------------------------------------------------------------------------------#

setGeneric("grab_GFUN",function(object){standardGeneric ("grab_GFUN")})

#' @rdname grab_GFUN-methods
#' @aliases grab_GFUN,m_estimation_basis,m_estimation_basis-method

setMethod(
  f = "grab_GFUN",
  signature = "m_estimation_basis",
  function(object){
    return(object@.GFUN)
  })

#------------------------------------------------------------------------------#
#' Shows the m_estimation_basis
#'
#' @rdname show-methods
#' @aliases show,m_estimation_basis,m_estimation_basis-method
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
    print(utils::head(object@.data))
    cat("Units: ", object@.units)

    invisible(NULL)
  })

#------------------------------------------------------------------------------#
# grab_basis_data generic
#
# Grabs the \code{.data} from an \code{\linkS4class{m_estimation_basis}} object
#
# @param object a \code{\linkS4class{m_estimation_basis}} object
#------------------------------------------------------------------------------#

setGeneric("grab_basis_data", function(object, ...) standardGeneric("grab_basis_data"))
setMethod("grab_basis_data", "m_estimation_basis", function(object) object@.data)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## geex class and methods ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#------------------------------------------------------------------------------#
#' geex S4 class
#'
#' @slot call the \code{m_estimate} call
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
  slots = c(call            = "call",
            basis           = "m_estimation_basis",
            rootFUN_results = "ANY",
            sandwich_components = "sandwich_components",
            GFUN            = "function",
            corrections     = "list",
            estimates       = "numeric",
            vcov            = "matrix"))


#------------------------------------------------------------------------------#
#' Shows the geex object
#'
#' @rdname show-methods
#' @aliases show,geex,geex-method
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
#' geex summary object
#'
#' @slot estFUN a \code{estimating-function}
#' @slot outer_args the \code{list} arguments passed to the \code{m_estimate} call
#' @slot inner_args the \code{list} arguments passed to the \code{m_estimate} call
#' @slot data the \code{data.frame} passed to the \code{m_estimate} call
#' @slot weights the \code{weights} passed to the \code{m_estimate} call
#' @slot nobs the number of observational units used to compute the M-estimator
#' @slot units the name of the variable identifying the observational units
#' @slot corrections a \code{list} of correction performed on \code{sandwich_components}
#' @slot estimates a \code{numeric} vector of parameter estimates
#' @slot vcov the empirical sandwich variance \code{matrix}
#' @export
#------------------------------------------------------------------------------#

setClass(
  Class = "geex_summary",
  slots = c(estFUN          = "function",
            outer_args      = "list",
            inner_args      = "list",
            data            = "data.frame",
            weights         = "numeric",
            nobs            = "numeric",
            units           = "character",
            corrections     = "list",
            estimates       = "numeric",
            vcov            = "matrix"))

#------------------------------------------------------------------------------#
#' Object Summaries
#'
#' @param object a \code{\linkS4class{geex}} object
#' @param keep_data keep the original data or not
#' @param keep_args keep the \code{outer_args} and \code{inner_args} passed to \code{estFUN} or not
#' @rdname summary-methods
#' @export
#' @examples
#' \dontrun{
#' library(geepack)
#' data('ohio')
#' glmfit  <- glm(resp ~ age, data = ohio,
#'               family = binomial(link = "logit"))
#' example_ee <- function(data, model){
#'   f <- grab_psiFUN(model, data)
#'   function(theta){
#'     f(theta)
#'   }
#' }
#' z  <- m_estimate(
#' estFUN = example_ee,
#' data = ohio,
#' compute_roots = FALSE,
#' units = 'id',
#' roots = coef(glmfit),
#' outer_args = list(model = glmfit))
#'
#' object.size(z)
#' object.size(summary(z))
#' object.size(summary(z, keep_data = FALSE))
#' object.size(summary(z, keep_data = FALSE, keep_args = FALSE))
#' }
#------------------------------------------------------------------------------#

setMethod(
  f         = "summary",
  signature = "geex",
  function(object, keep_data = TRUE, keep_args = TRUE){
    methods::new("geex_summary",
        estFUN      = object@basis@.estFUN,
        outer_args  = if(keep_args) object@basis@.outer_args else list(),
        inner_args  = if(keep_args) object@basis@.inner_args else list(),
        data        = if(keep_data) object@basis@.data else data.frame(),
        units       = object@basis@.units,
        weights     = object@basis@.weights,
        nobs        = nobs(object),
        corrections = object@corrections,
        estimates   = object@estimates,
        vcov        = object@vcov)
  }
)

#------------------------------------------------------------------------------#
#' Shows the geex_summary object
#'
#' @rdname show-methods
#' @aliases show,geex_summary,geex_summary-method
#' @export
#------------------------------------------------------------------------------#

setMethod(
  "show",
  signature = "geex_summary",
  definition = function(object) {
    cat("M-estimation results based on the estimating function:\n", sep = "")
    print(body(object@estFUN))
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
#' @rdname vcov-methods
#' @aliases vcov,geex,geex-method
#' @export
#' @examples
#' ex_eeFUN <- function(data){
#'  function(theta){
#'    with(data,
#'     c(Y1 - theta[1],
#'      (Y1 - theta[1])^2 - theta[2] ))
#' }}
#'
#' results <- m_estimate(
#'  estFUN = ex_eeFUN,
#'  data  = geexex,
#'  root_control = setup_root_control(start = c(1,1)))
#'
#' vcov(results)

setMethod(
  "vcov",
  signature = "geex",
  definition = function(object) {
    object@vcov
  })

#' @rdname vcov-methods
#' @aliases vcov,geex_summary,geex_summary-method
#' @export

setMethod(
  "vcov",
  signature = "geex_summary",
  definition = function(object) {
    object@vcov
  })

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Gets the parameter estimates from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @rdname coef-methods
#' @aliases coef,geex,geex-method
#' @export
#' @examples
#' ex_eeFUN <- function(data){
#'  function(theta){
#'    with(data,
#'     c(Y1 - theta[1],
#'      (Y1 - theta[1])^2 - theta[2] ))
#' }}
#'
#' results <- m_estimate(
#'  estFUN = ex_eeFUN,
#'  data  = geexex,
#'  root_control = setup_root_control(start = c(1,1)))
#'
#' coef(results)

setMethod(
  "coef",
  signature = "geex",
  definition = function(object) {
    object@estimates
  })

#' @rdname coef-methods
#' @aliases coef,geex_summary-method
#' @export

setMethod(
  "coef",
  signature = "geex_summary",
  definition = function(object) {
    object@estimates
  })
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#' Gets the parameter estimates matrix from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @param ... arguments passed to other methods
#' @docType methods
#' @rdname roots-methods
#' @export
#' @examples
#' ex_eeFUN <- function(data){
#'  function(theta){
#'    with(data,
#'     c(Y1 - theta[1],
#'      (Y1 - theta[1])^2 - theta[2] ))
#' }}
#'
#' results <- m_estimate(
#'  estFUN = ex_eeFUN,
#'  data  = geexex,
#'  root_control = setup_root_control(start = c(1,1)))
#'
#' roots(results)
#------------------------------------------------------------------------------#

setGeneric("roots", function(object, ...) standardGeneric("roots"))

#' @rdname roots-methods
#' @aliases roots,geex,geex-method

setMethod(
  "roots",
  signature = "geex",
  definition = function(object) {
    object@estimates
  })

#' @rdname roots-methods
#' @aliases roots,geex,geex-method

setMethod(
  "roots",
  signature = "geex_summary",
  definition = function(object) {
    object@estimates
  })

#------------------------------------------------------------------------------#
#' Gets the corrections from a geex object
#'
#' @param object a \code{\linkS4class{geex}} object
#' @param ... arguments passed to other methods
#' @docType methods
#' @rdname get_corrections-methods
#' @export
#' @examples
#' myee <- function(data){
#'  function(theta){
#'    c(data$Y1 - theta[1],
#'    (data$Y1 - theta[1])^2 - theta[2])
#'   }
#' }
#'
#' results <- m_estimate(
#'    estFUN = myee,
#'    data = geexex,
#'    root_control = setup_root_control(start = c(1,1)),
#'    corrections  = list(
#'      bias_correction_.1 = correction(fay_bias_correction, b = .1),
#'      bias_correction_.3 = correction(fay_bias_correction, b = .3))
#'    )
#'
#' get_corrections(results)
#------------------------------------------------------------------------------#

setGeneric("get_corrections", function(object, ...) standardGeneric("get_corrections"))

#' @rdname get_corrections-methods
#' @aliases get_corrections,m_estimation_basis,m_estimation_basis-method
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

#' @rdname get_corrections-methods
#' @export

setMethod(
  "get_corrections",
  signature = "geex_summary",
  definition = function(object) {
    if(length(object@corrections) == 0){
      "No corrections were performed on this object"
    } else {
      object@corrections
    }
  })

#' @rdname grab_psiFUN_list-methods
#' @aliases grab_psiFUN_list,geex,geex-method

setMethod(
  f = "grab_psiFUN_list",
  signature = "geex",
  function(object){
    return(object@sandwich_components@.psiFUN_list)
  })

#' @rdname grab_GFUN-methods
#' @aliases grab_GFUN,geex,geex-method

setMethod(
  f = "grab_GFUN",
  signature = "geex",
  function(object){
    return(object@basis@.GFUN)
  })

#------------------------------------------------------------------------------#
#' Extract the number observations
#'
#' @param object a \code{\linkS4class{geex}} object
#' @rdname nobs-methods
#' @aliases nobs,geex,geex-method
#' @export
#' @examples
#' \dontrun{
#' library(geepack)
#' data('ohio')
#'
#' glmfit  <- glm(resp ~ age, data = ohio,
#'               family = binomial(link = "logit"))
#' z  <- m_estimate(
#'   estFUN = example_ee,
#'   data = ohio,
#'   compute_roots = FALSE,
#'   units = 'id',
#'   roots = coef(glmfit),
#'   outer_args = list(model = glmfit))
#'
#' nobs(z)
#' }

setMethod(
  f         = "nobs",
  signature = "geex",
  function(object){
    if(length(object@basis@.units) == 0){
      nrow(object@basis@.data)
    } else {
      length(unique(object@basis@.data[[object@basis@.units]]))
    }
  }
)

#' @rdname nobs-methods
#' @aliases nobs,geex_summary,geex_summary-method
#' @export

setMethod(
  f         = "nobs",
  signature = "geex_summary",
  function(object){
    object@nobs
  }
)

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Extract Model weights
#'
#' @param object a \code{\linkS4class{geex}} object
#' @rdname weights-methods
#' @aliases weights,geex,geex-method
#' @export

setMethod(
  f         = "weights",
  signature = "geex",
  function(object){
    object@basis@.weights
  }
)

#' @rdname weights-methods
#' @aliases weights,geex_summary,geex_summary-method
#' @export

setMethod(
  f         = "weights",
  signature = "geex_summary",
  function(object){
    object@weights
  }
)

#------------------------------------------------------------------------------#


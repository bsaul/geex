#------------------------------------------------------------------------------#
# check_** description:
# Functions for use within geex computations that check that functions or other
# objects conform to what downstream functions expect.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Check and convert object to an array
#
# Used in, e.g., \code{compute_matrices} to either (a) confirm that an object is
# an array or (b) if the object is numeric or a matrix, then convert to an array,
# or (c) if the object is not an array, numeric, or matrix, then return an error.
#
# @param object the object to check whether it is an array
# @return an array - either the orginal object or the given object converted
# to an array
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
# Check that a list of corrections conforms to geex standards
#
# @param corrections a list of small sample corrections to perform
# @export
#------------------------------------------------------------------------------#
check_corrections <- function(corrections){
  corrs <- names(corrections)

  out <- lapply(seq_along(corrs), function(i){
    items <- names(corrections[[i]])
    if(!('correctFUN' %in% items)){
      stop(paste0('correctFUN is not specified in ', corrs[i]))
    }
    if(any(!(items %in% c('correctFUN', 'correctFUN_control')))){
      warning(paste0('Additional list items in the ', corrs[i], ' correction will be ignored. Pass all additional arguments to correctFUN using the correctFUN_control item.'))
    }
  })
  invisible(out)
}

#------------------------------------------------------------------------------#
# Check an eeFUN object
#
# Checks that eeFUN returns a function. May be developed to perform additional
# checks in the future.
#
# @param geex_list a list of \code{eeFUN}, \code{splitdt}, \code{inner_eeargs}, and
# \code{outer_eeargs}
# @export
#------------------------------------------------------------------------------#

check_eeFUN <- function(geex_list){
  # Check: eeFUN returns a function
  f <- do.call(geex_list$eeFUN, args = append(list(data = geex_list$splitdt[[1]]),
                                              geex_list$outer_eeargs))
  if(!is.function(f)){
    stop('eeFUN does not return a function')
  }
}




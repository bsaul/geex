#------------------------------------------------------------------------------#
# check_** description:
# Functions for use within geex computations that check that functions or other
# objects conform to what downstream functions expect.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Check and convert object to an array
#
# Used in, e.g., \code{process_matrix_list} to either (a) confirm that an object is
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

#------------------------------------------------------------------------------#
# compute_** description:
# Functions that take in numeric values (e.g. numeric(), matrices, arrays)
# and compute numeric values.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Compute empirical sandwich covariate estimator
#'
#' Computes \eqn{\Sigma = A^{-1} B (A^{-1})^T }{\Sigma = A^{-1} B (A^{-1})^T} with
#' provided \eqn{A}{A} and \eqn{B}{B} matrices.
#'
#' @param A a matrix, generally the \eqn{A}{A} matrix returned in the list of matrices from
#'   \code{\link{compute_matrices}}
#' @param B a matrix, generally the \eqn{B}{B} matrix returned in the list of matrices from
#'   \code{\link{compute_matrices}}
#' @export
#------------------------------------------------------------------------------#

compute_sigma <- function(A, B){
  Ainv <- solve(A)
  Ainv %*% B %*% t(Ainv)
}


#------------------------------------------------------------------------------#
# Compute the sum of  a list of matrices to sum
#
# @param l a list of matrices
# @param w a numeric vector of weights
#------------------------------------------------------------------------------#

compute_sum_of_matrix_list <- function(.l, .w = numeric(0)){
  M_i_pre   <- if(length(.w) > 0){ Map(`*`, .l, .w) } else .l
  Reduce(`+`, M_i_pre)
}

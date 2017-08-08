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
#' @param A a matrix, generally the \code{.A} slot in a
#' \code{\linkS4class{sandwich_components}} object created in
#' \code{\link{estimate_sandwich_components}}
#' @param B a matrix, generally the \code{.B} slot in a
#' \code{\linkS4class{sandwich_components}} object created in
#' \code{\link{estimate_sandwich_components}}
#'
#' @return the \code{matrix} \code{Ainv \%*\% B \%*\% t(Ainv)}
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

compute_sum_of_list <- function(.l, .w = numeric(0)){
  M_i_pre   <- if(length(.w) > 0){ Map(`*`, .l, .w) } else .l
  Reduce(`+`, M_i_pre)
}

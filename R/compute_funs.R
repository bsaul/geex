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
#' \code{\link{estimate_sandwich_matrices}}
#' @param B a matrix, generally the \code{.B} slot in a
#' \code{\linkS4class{sandwich_components}} object created in
#' \code{\link{estimate_sandwich_matrices}}
#'
#' @return the \code{matrix} \code{Ainv \%*\% B \%*\% t(Ainv)}
#' @export
#' @examples
#' A <- diag(2, nrow = 2, ncol = 2)
#' B <- matrix(4, nrow = 2, ncol = 2)
#' compute_sigma(A = A, B = B)
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
  dimw <- length(.w)

  M_i_pre   <- if(dimw > 0){
    stopifnot(dimw == length(.l))
    Map(`*`, .l, .w)
  } else .l
  Reduce(`+`, M_i_pre)
}

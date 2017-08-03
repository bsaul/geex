#------------------------------------------------------------------------------#
# compute_** description:
# Functions that take in stuff and compute numeric values. The workhorses of
# performing M-estimation
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

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
#' @param solver the function used to compute the inverse of \code{A}, Defaults
#' to \code{\link{solve}}
#'
#' @return the \code{matrix} \code{Ainv \%*\% B \%*\% t(Ainv)}
#' @export
#' @examples
#' A <- diag(2, nrow = 2, ncol = 2)
#' B <- matrix(4, nrow = 2, ncol = 2)
#' compute_sigma(A = A, B = B)
#------------------------------------------------------------------------------#
compute_sigma <- function(A, B, solver = solve) {
  Ainv <- solver(A)
  Ainv %*% B %*% t(Ainv)
}

#------------------------------------------------------------------------------#
#' Compute the sum of  a list of matrices to sum
#'
#' @param .l a list of matrices
#' @param .w a numeric vector of weights
#' @keywords internal
#------------------------------------------------------------------------------#

compute_sum_of_list <- function(.l, .w = numeric(0)) {
  dimw <- length(.w)

  M_i_pre <- if (dimw > 0) {
    stopifnot(dimw == length(.l))
    Map(`*`, .l, .w)
  } else {
    .l
  }
  Reduce(`+`, M_i_pre)
}

#------------------------------------------------------------------------------#
#' Compute the sum of  a list of matrices to sum
#'
#' @param .l a list of matrices
#' @param .w a numeric vector of weights
#' @param .wFUN a function of \code{i}, \code{j}, and (optionally) additional
#' arguments
#' @param ... additional arguments passed to \code{.wFUN}
#'
#' Either \code{.w} or \code{.wFUN} must be specified but not both.
#'
#' @importFrom methods formalArgs
#' @keywords internal
#------------------------------------------------------------------------------#

compute_pairwise_sum_of_list <- function(.l, .w = NULL, .wFUN = NULL, ...) {
  use_w <- !missing(.w)
  use_wFUN <- !missing(.wFUN)

  if ((use_w && use_wFUN) || (!use_w && !use_wFUN)) {
    stop("Either a vector of weights (.w) or a function (.wFUN) must be specified and not both.")
  }

  if (use_w && (length(.l) != length(.w))) {
    stop("The length of the weight vector must equal the length of the list.")
  }

  if (use_wFUN && all(methods::formalArgs(.wFUN)[1:2] != c("i", "j"))) {
    stop("The first two arguments of .wFUN must be i and j.")
  }

  lapply(seq_along(.l), function(i) {
    lapply(seq_along(.l), function(j) {
      if (use_wFUN) {
        .w <- .wFUN(i, j, ...)
      }
      .w * tcrossprod(.l[[i]], .l[[j]])
      # w * .l[[i]] %*% t(.l[[j]])
    }) -> hold1
    compute_sum_of_list(hold1)
  }) -> hold2
  compute_sum_of_list(hold2)
}

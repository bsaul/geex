#------------------------------------------------------------------------------#
# correct_by_** description:
# Functions that correct the output of compute_matrices *by* some manipulation.
# Also includes necessary (non-exported) internal functions.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Make corrections to basic sandwich estimators
#
# @param mats list of matrices resulting from \code{\link{compute_matrices}}
# @param corrections a list with one sublist for each correction to make. The
# sublist must contain at least \code{correctFUN} and optionally additional
# arguments in \code{correctFUN_control}.
# @return a list with corrected values. The list retains the names of the input
#  \code{corrections}
# @export
#------------------------------------------------------------------------------#
make_corrections <- function(components, corrections){
  lapply(corrections, function(correction){
     correct_by(.components = components, .correct_control = correction)
  })
}

###############################################################################
# correct_by functions provided by geex ####
###############################################################################

#------------------------------------------------------------------------------#
#' Correct sandwich components
#'
#' @param .components an object of class \code{\linkS4class{sandwich_components}}
#' @param .correct_control an object of class \code{\linkS4class{correct_control}}
#'
#' @return the result of \code{.FUN} in \code{.correct_control}.
#' @export
#------------------------------------------------------------------------------#

correct_by <- function(.components, .correct_control){
  f <- grab_FUN(.correct_control)
  opts <- grab_options(.correct_control)
  do.call(f, args = append(list(components = .components), opts))
}


#------------------------------------------------------------------------------#
#' Creates a correct_control object
#'
#' @param FUN a correction to perform. \code{components} must be the
#' first argument
#' @param ,,, additional arguments passed to \code{FUN}
#' @return a \code{\linkS4class{correct_control}} object
#' @export
#------------------------------------------------------------------------------#

correction <- function(FUN, ...){
  dots <- list(...)
  new(Class="correct_control",
      .FUN = FUN,
      .options   = dots)
}

#------------------------------------------------------------------------------#
#' Correct sandwich variance estimator byFay's bias correction
#'
#' Computes the bias corrected sandwich covariance matrix described in Fay and
#' Graubard (2001). See \code{vignette("05_finite_sample_corrections", package = "geex")}
#' for further information.
#'
#' @param components an object of class \code{\linkS4class{sandwich_components}}
#' @param b a numeric value < 1. Defaults to 0.75 as in Fay.
#' @return a corrected covariance matrix
#' @references Fay, M. P., & Graubard, B. I. (2001). Small-Sample adjustments for
#' Wald-type tests using sandwich estimators. Biometrics, 57(4), 1198-1206
#' @export
#------------------------------------------------------------------------------#

fay_bias_correction <- function(components, b = 0.75){
  fay_bias_correction_partial(components, b = b) ->  corrected_matrices

  compute_sigma(A = grab_bread(components), B = corrected_matrices$Bbc)
}

#------------------------------------------------------------------------------#
#' Correct sandwich variance inference by Fay's degrees of freedom correction
#'
#' Computes the degrees of freedom correction described in Fay and
#' Graubard (2001). See \code{vignette("05_finite_sample_corrections", package = "geex")}
#' for further information.
#'
#' @inheritParams fay_bias_correction
#' @param L a k x p matrix where p is the dimension of theta
#' @param version either 1 or 2, corresponding to hat(d) or tilde(d), respectively
#' @references Fay, M. P., & Graubard, B. I. (2001). Small-Sample adjustments for
#' Wald-type tests using sandwich estimators. Biometrics, 57(4), 1198-1206
#' @return a scalar corresponding to the estimated degrees of freedom
#' @export
#------------------------------------------------------------------------------#

fay_df_correction <- function(components, b = .75, L, version){

  ## Prepare necessary matrices ##
  bias_mats <- fay_bias_correction_partial(components, b = b)
  A   <- grab_bread(components)
  A_i <- grab_bread_list(components)
  df_prep   <- df_correction_prep(L = L, A = A, A_i = A_i, H_i = bias_mats$H_i)

  ## Compute DF corrections ##
  if(version == 1){
    out <- df_correction_1(df_prep$A_d, df_prep$C)
  } else if(version == 2){
    out <- df_correction_2(A = A, A_i = A_i, C = df_prep$C, L = L, Bbc = bias_mats$Bbc)
  }
  out
}

###############################################################################
# functions internal to provided correct_by ####
###############################################################################


#------------------------------------------------------------------------------#
# Compute the matrices necessary for Fay's bias correction
#
# @inheritParams fay_bias_correction
# @return a list with
# \itemize{
# \item H_i - the H_i matrix
# \item Bbc - the 'bias corrected' meat matrix
# }
# @export
#------------------------------------------------------------------------------#

fay_bias_correction_partial <- function(components, b){

  A   <- grab_bread(components)
  A_i <- grab_bread_list(components)
  B_i <- grab_meat_list(components)

  Ainv <- solve(A)
  H_i <- lapply(A_i, function(Ai){
    diag( (1 - pmin(b, diag(Ai %*% Ainv) ))^(-0.5) )
  })

  stopifnot(length(B_i) == length(H_i))

  Bbc_i <- lapply(seq_along(B_i), function(i){
    H_i[[i]] %*% B_i[[i]] %*% H_i[[i]]
  })
  Bbc <- compute_sum_of_list(Bbc_i) #apply(simplify2array(Bbc_i), 1:2, sum)

  list(H_i = H_i, Bbc = Bbc)
}

#------------------------------------------------------------------------------#
# Preparations for Fay's degrees of freedom corrections
#
# @param H_i a list of bias adjusted matrices
# @inheritParams fay_bias_correction
# @inheritParams fay_df_correction
# @export
#------------------------------------------------------------------------------#
df_correction_prep <- function(L, A, A_i, H_i){
  p <- ncol(A)
  m <- length(A_i)

  Ainv <- solve(A)
  II   <- diag(1, p*m)
  AA   <- do.call(rbind, A_i)
  calI <- do.call(cbind, args = lapply(1:m, function(i) diag(1, p) ))
  G    <- II - (AA %*% solve(A) %*% calI)

  M_i  <- lapply(H_i, function(mat){
    mat %*% Ainv %*% L %*% t(L) %*% t(Ainv) %*% mat
  })
  M    <- Matrix::bdiag(M_i)

  C    <- t(G) %*% M %*% G

  A_diag  <- Matrix::bdiag(A_i)

  list(A_d = A_diag, C = C)
}

#------------------------------------------------------------------------------#
# Estimate Fay's degrees of freedom corrections 1
#
# @param A_d see \code{\link{df_correction_prep}}
# @param C see \code{\link{df_correction_prep}}
#
# @export
#------------------------------------------------------------------------------#
df_correction_1 <- function(A_d, C){
  estimate_df(A = A_d, C = C)
}

#------------------------------------------------------------------------------#
# Estimate Fay's degrees of freedom correction 2
#
# @inheritParams fay_bias_correction
# @inheritParams fay_df_correction
# @param C same as the C matrix in the Fay (2001) notation Section 2.3
# @param Bbc the bias corrected "B" matrix
#
# @export
#------------------------------------------------------------------------------#
df_correction_2 <- function(A, A_i, C, L, Bbc){
  w_i  <- lapply(seq_along(A_i), function(i) {
    # exclude the ith element
    Oi <- apply(simplify2array(A_i[-i]), 1:2, sum)
    t(L) %*% (solve(Oi) - solve(A) ) %*% L
  })
  wbar <- sum(unlist(w_i))

  Abc_i <- lapply(w_i, function(w){
    as.numeric(w/wbar) * Bbc
  })
  Abc  <- Matrix::bdiag(Abc_i)

  estimate_df(A = Abc, C = C)
}

#------------------------------------------------------------------------------#
# Estimate Fay's degrees of freedom corrections
#
# @inheritParams df_correction_2
# @inheritParams fay_df_correction
# @export
#------------------------------------------------------------------------------#

estimate_df <- function(A, C){
  AC <- A %*% C
  (sum(Matrix::diag(AC)))^2 / sum(Matrix::diag(AC %*% AC))
}


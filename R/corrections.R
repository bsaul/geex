#------------------------------------------------------------------------------#
#' Estimate Fay's bias correction
#'
#' @param m the number of groups
#' @param A the bread matrix
#' @param Ai a list of bread matrices per group
#' @param Bi a list of meat matrices per group
#' @param b a numeric value < 1. Defaults to 0.75 as in Fay
#' @return a list with
#' \itemize{
#' \item H_i - the H_i matrix
#' \item Bbc - the 'bias corrected' meat matrix
#' }
# @export
#------------------------------------------------------------------------------#

bias_correction <- function(m, A, Ai, Bi, b = 0.75){

  H_i <- lapply(Ai, function(m){
    diag( (1 - pmin(b, diag(m %*% solve(A)) ) )^(-0.5) )
  })

  Bbc_i <- lapply(1:m, function(i){
    H_i[[i]] %*% Bi[[i]] %*% H_i[[i]]
  })
  Bbc   <- apply(simplify2array(Bbc_i), 1:2, sum)

  list(H_i = H_i, Bbc = Bbc)
}

#------------------------------------------------------------------------------#
#' Preparations for Fay's degrees of freedom corrections
#'
# @export
#------------------------------------------------------------------------------#
df_correction_prep <- function(m, L, A, A_i, H_i){
  p <- ncol(A)

  II   <- diag(1, p*m)
  AA   <- do.call(rbind, A_i)
  calI <- do.call(cbind, args = lapply(1:m, function(i) diag(1, p) ))
  G    <- II - (AA %*% solve(A) %*% calI)

  M_i  <- lapply(H_i, function(mat){
    mat %*% solve(A) %*% L %*% t(L) %*% t(solve(A)) %*% mat
  })
  M    <- Matrix::bdiag(M_i)

  C    <- t(G) %*% M %*% G

  A_d  <- Matrix::bdiag(A_i)

  list(A_d = A_d, C = C)
}

#------------------------------------------------------------------------------#
#' Estimate Fay's degrees of freedom corrections 1
#'
#' @param A_d see \code{\link{df_correction_prep}}
#' @param C see \code{\link{df_correction_prep}}
#'
# @export
#------------------------------------------------------------------------------#
df_correction_1 <- function(A_d, C){
  estimate_df(A = A_d, C = C)
}

#------------------------------------------------------------------------------#
#' Estimate Fay's degrees of freedom correction 2
#'
# @export
#------------------------------------------------------------------------------#
df_correction_2 <- function(m, A, A_i, C, L,  Bbc){
  w_i  <- lapply(1:m, function(i) {
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
#' Estimate Fay's degrees of freedom corrections
#'
# @export
#------------------------------------------------------------------------------#

estimate_df <- function(A, C){
  AC <- A %*% C
  (sum(Matrix::diag(AC)))^2 / sum(Matrix::diag(AC %*% AC))
}

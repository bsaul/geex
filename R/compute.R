#------------------------------------------------------------------------------#
#' Handles the case where there a single estimating equation
#' @export
#------------------------------------------------------------------------------#

check_array <- function(vals){
  if(!is.array(vals)){
    array(vals, dim = c(1, 1, length(vals)))
  } else {
    vals
  }
}


#------------------------------------------------------------------------------#
#' Compute roots for a set of estimating equations
#' @export
#------------------------------------------------------------------------------#

eeroot <- function(obj, start, root_options = NULL, ...){

  psi_i <- lapply(obj$splitdt, function(data_i){
    obj$eeFUN(data = data_i, ...)
  })

  psi <- function(theta){
    psii <- lapply(psi_i, function(f) f(theta))
    apply(check_array(simplify2array(psii)), 1, sum)
  }

  rargs <- append(root_options, list(f = psi, start = start))

  do.call(rootSolve::multiroot, args = rargs)
}

#------------------------------------------------------------------------------#
#' Compute component matrices for covariance matrix for a set of
#' estimating equations
#'
#' @export
#------------------------------------------------------------------------------#

compute_matrices <- function(obj,
                             contrast = NULL,
                             theta,
                             corrections = NULL,
                             correction_options = list(),
                             numDeriv_options = list(method = 'Richardson'),
                             ...){
  # Warnings
  if('bias' %in% corrections & is.null(correction_options$b)){
    stop('b argument must be present if using bias correction')
  }

  correct_bias <- any(c('bias') %in% corrections)
  correct_df   <- any(c('df') %in% corrections)

  with(obj, {
    m <- length(splitdt)

    # Create list of estimating eqn functions per unit
    psi_i <- lapply(splitdt, function(data_i){
      eeFUN(data = data_i, ...)
    })

    # Compute the negative of the derivative matrix of estimating eqn functions
    # (the information matrix)
    A_i <- lapply(psi_i, function(ee){
      args <- append(list(fun = ee, x = theta), numDeriv_options)
      val  <- do.call(numDeriv::jacobian, args = args)
      -val
    })
    A_i <- check_array(simplify2array(A_i))
    A   <- apply(A_i, 1:2, sum)

    # Compute outer product of observed estimating eqns
    B_i <- lapply(psi_i, function(ee) ee(theta) %*% t(ee(theta)) )
    B   <- apply(check_array(simplify2array(B_i)), 1:2, sum)

    out <-  list(A = A, A_i = A_i, B = B, B_i = B_i)

    #### Compute corrections ####

    # Bias correction #
    if(correct_bias|correct_df){
      bias_try <- try(bias_correction(m = m, A = A, Ai = A_i, Bi = B_i,
                                      b = correction_options$b),
                       silent = TRUE)

      if(is(bias_try, 'try-error')){
        bias_fail <- TRUE
      } else {
        bias_fail <- FALSE
        bias_mats <- bias_try
        Bbc <- bias_mats$Bbc
        H_i <- bias_mats$H_i
        out$Bbc <- Bbc
      }

      out$bias_fail <- bias_fail
    }

    # Degrees of Freedom corrections #
    if(correct_df){
      if(is.null(contrast)){
        stop('contrast must be specified for df correction')
      }

      if(!bias_fail){

        df_prep <- df_correction_prep(m = m, L = contrast, A = A, A_i = A_i, H_i = H_i)

        # DF correction 1 #
        df1 <- df_correction_1(A_d = df_prep$A_d, C = df_prep$C)

        # DF correction 2 #
        df2_try <- try(df_correction_2(m = m, A = A, A_i = A_i, C = df_prep$C,
                                       L = contrast, Bbc = Bbc), silent = TRUE)
        if(is(df2_try, 'try-error')){
          df2 <- NA_real_
        } else {
          df2 <- df2_try
        }

        out$df1 <- df1
        out$df2 <- df2
      }
    }

    out
  })
}

#------------------------------------------------------------------------------#
#' Estimate Fay's bias correction
#'
#' @export
#------------------------------------------------------------------------------#

bias_correction <- function(m, A, Ai, Bi, b){

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
#' @export
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
#' @export
#------------------------------------------------------------------------------#
df_correction_1 <- function(A_d, C){
  estimate_df(A = A_d, C = C)
}


#------------------------------------------------------------------------------#
#' Estimate Fay's degrees of freedom correction 2
#'
#' @export
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
#' @export
#------------------------------------------------------------------------------#

estimate_df <- function(A, C){
  AC <- A %*% C
  (sum(Matrix::diag(AC)))^2 / sum(Matrix::diag(AC %*% AC))
}

#------------------------------------------------------------------------------#
#' Compute covariance matrix for set of estimating equations
#'
#' @export
#------------------------------------------------------------------------------#

compute_sigma <- function(matrices, corrections = NULL){
  with(matrices,
       if(any('bias' %in% corrections)){
         solve(A) %*% Bbc %*% t(solve(A))
       } else {
         solve(A) %*% B %*% t(solve(A))
       })
}

#------------------------------------------------------------------------------#
#' ee
#'
#' @export
#------------------------------------------------------------------------------#

estimate_equations <- function(eeFUN,
                               data,
                               units,
                               corrections = NULL,
                               correction_options = list(),
                               numDeriv_options = list(method = 'Richardson'),
                               rootsolver = rootSolve::multiroot,
                               rootsolver_options = NULL,
                               findroots  = TRUE,
                               roots = NULL,
                               ...){

  ## Warnings ##
  if(missing(roots) & findroots){
    stop('If findroots = TRUE, then starting values for the rootsolver must be specified in roots argument.')
  }

  split_data <- split(x = data, f = data[[units]] )
  eeobj      <- list(eeFUN = eeFUN, splitdt = split_data)

  ## Compute estimating equation roots ##
  if(findroots){
    eesolved <- eeroot(eeobj, start = roots,
                       root_options = rootsolver_options,
                       ...)
    theta_hat <- eesolved$root
  } else {
    theta_hat <- roots
  }

  ## Compute variance estimates ##
  mats <- compute_matrices(obj   = eeobj,
                           theta = theta_hat,
                           corrections = corrections,
                           numDeriv_options = numDeriv_options,
                           correction_options = correction_options,
                           ...)

  Sigma_hat <- compute_sigma(mats)

  list(parameters = theta_hat, vcov = Sigma_hat)
}


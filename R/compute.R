
compute_matrices <- function(obj,
                             contrast = NULL,
                             corrections = NULL,
                             correction_options = list(),
                             numDeriv_options = list(method = 'simple'),
                             ...){
  # Warnings
  if('bias' %in% corrections & is.null(correction_options$b)){
    stop('b argument must be present if using bias correction')
  }

  with(obj, {
    m <- length(splitdt)

    # Create list of estimating eqn functions per unit
    psi_i <- lapply(splitdt, function(data_i){
      psi_i <- eeFUN(data = data_i, ...)
    })

    # Compute the negative of the derivative matrix of estimating eqn functions
    A_i <- lapply(psi_i, function(ee){
      args <- append(list(fun = ee, x = theta), numDeriv_options)
      val  <- do.call(numDeriv::jacobian, args = args)
      -val
    })
    A   <- apply(simplify2array(A_i), 1:2, sum)

    # Compute outer product of observed estimating eqns
    B_i <- lapply(psi_i, function(ee) ee(theta) %*% t(ee(theta)) )
    B   <- apply(simplify2array(B_i), 1:2, sum)

    out <-  list(A = A, A_i = A_i, B = B, B_i = B_i)

    #### Compute corrections ####

    # Bias correction
    if(any(c('bias', 'df') %in% corrections)){
      b   <- correction_options$b
      H_i <- lapply(A_i, function(m){
        diag( (1 - pmin(b, diag(m %*% solve(A)) ) )^(-0.5) )
      })

      Bbc_i <- lapply(1:m, function(i){
        H_i[[i]] %*% B_i[[i]] %*% H_i[[i]]
      })
      Bbc   <- apply(simplify2array(Bbc_i), 1:2, sum)

      out$Bbc <- Bbc
    }

    # Degrees of Freedom corrections
    if('df' %in% corrections){
      if(is.null(contrast)){
        stop('contrast must be specified for df correction')
      }

      L <- contrast
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

      out$df1  <- estimate_df(A = A_d, C = C)
      out$df2  <- estimate_df(A = Abc, C = C)
    }
    out
  })
}

estimate_df <- function(A, C){
  AC <- A %*% C
  (sum(Matrix::diag(AC)))^2 / sum(Matrix::diag(AC %*% AC))
}

compute_sigma <- function(matrices, ...){
  with(matrices, solve(A) %*% B %*% t(solve(A)))
}

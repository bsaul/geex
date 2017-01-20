#------------------------------------------------------------------------------#
#' Compute roots for a set of estimating equations
#'
#' @param geex_list a list containing \code{splitdt} (a \code{data.frame} that
#' has been \code{\link[base]{split}} by the grouping variable) and \code{eeFUN}
#' (see \code{\link{estimate_equations}})
#' @param start vector with length of the number of parameters to find. Passed to
#' \code{\link[rootSolve]{multiroot}} if not NULL. Defaults to NULL.
#' @param rootsolver the function used to find roots of the estimating equations.
#' Defaults to \code{\link[rootSolve]{multiroot}}.
#' @param root_options a list of options to be passed to the \code{rootsolver}
#' function
#' @param ... additional arguments passed to \code{geex_list$eeFUN}
#' @return the output of the \code{rootsolver} function
#' @export
#------------------------------------------------------------------------------#

eeroot <- function(geex_list,
                   start        = NULL,
                   rootsolver   = rootSolve::multiroot,
                   root_options = NULL,
                   ...){

  # Create estimating equation functions per group
  psi_i <- lapply(geex_list$splitdt, function(data_i){
    geex_list$eeFUN(data = data_i, ...)
  })

  # Create psi function that sums over all ee funs
  psi <- function(theta){
    psii <- lapply(psi_i, function(f) f(theta))
    apply(check_array(simplify2array(psii)), 1, sum)
  }

  # Find roots of psi
  rargs <- append(root_options, list(f = psi, start = start))
  do.call(rootsolver, args = rargs)
}

#------------------------------------------------------------------------------#
#' Compute component matrices for covariance matrix
#'
#' For a given set of estimating equations computes the 'meat' and 'bread' matrices
#' necessary to compute the covariance matrix.
#'
#' @param geex_list a list containing \code{splitdt} (a \code{data.frame} that
#' has been \code{\link[base]{split}} by the grouping variable) and \code{eeFUN}
#' (see \code{\link{estimate_equations}})
#' @param theta vector of parameters passed to \code{geex_list$eeFUN}.
#' @param numDeriv_options a list of options passed to \code{\link[numDeriv]{jacobian}}.
#' @param ... additional arguments passed to \code{geex_list$eeFUN}.
#' @return a list with
#' \itemize{
#' \item A - the 'bread' matrix
#' \item B - the 'meat' matrix
#' \item A_i - the 'bread' matrix for each group
#' \item B_i - the 'meat' matrix for each group
#' }
#'
#' @export
#------------------------------------------------------------------------------#

compute_matrices <- function(geex_list,
                             theta,
                             numDeriv_options = list(method = 'Richardson'),
                             silent = TRUE,
                             ...){

  with(geex_list, {
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
    A_i_array <- check_array(simplify2array(A_i))
    A   <- apply(A_i_array, 1:2, sum)

    # Compute outer product of observed estimating eqns
    B_i <- lapply(psi_i, function(ee) ee(theta) %*% t(ee(theta)) )
    B   <- apply(check_array(simplify2array(B_i)), 1:2, sum)

    list(A = A, A_i = A_i, B = B, B_i = B_i)
  })
}



#------------------------------------------------------------------------------#
#' Compute covariance matrix for set of estimating equations
#'
#' @param matrices the list of matrices returned from \code{\link{compute_matrices}}
#' @param corrections character vector of corrections to perform
#'
#' @export
#------------------------------------------------------------------------------#

compute_sigma <- function(A, B){
  solve(A) %*% B %*% t(solve(A))
}

#------------------------------------------------------------------------------#
#' Estimate parameters and their covariance from a set of estimating equations
#'
#' @param eeFUN a function that takes in group-level data and returns a function
#' that takes parameters as its first argument
#' @param data a data.frame
#' @param units a string identifying the grouping variable in \code{data}
#' @param corrections character vector of small sample corrections to apply
#' @param correction_options a list of options for the small sample corrections
#' @param numDeriv_options a list of options for \code{\link[numDeriv]{jacobian}}
#' @param rootsolver the function that will find roots of the estimating equations
#' when \code{findroots = TRUE}.
#' @param rootsolver_options a list of options for the \code{rootsolver} function
#' @param findroots whether or not to find the roots of the estimating equations.
#' Defaults to \code{TRUE}.
#' @param roots a numeric vector containing either starting values for the roots when using
#' the default \code{rootsolver} or roots that have been estimated elsewhere
#' @param ... additional arguments passed to the \code{eeFUN}.
#' @return a list with the following
#' \itemize{
#' \item \code{parameters} - a vector of estimated parameters
#' \item \code{vcov} - the variance-covariance matrix for the parameters
#' }
#' @export
#------------------------------------------------------------------------------#

estimate_equations <- function(eeFUN,
                               data,
                               units,
                               corrections_list = NULL,
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
  geex_list  <- list(eeFUN = eeFUN, splitdt = split_data)

  ## Compute estimating equation roots ##
  if(findroots){
    eesolved <- eeroot(geex_list,
                       start        = roots,
                       rootsolver   = rootsolver,
                       root_options = rootsolver_options,
                       ...)
    theta_hat <- eesolved$root
  } else {
    theta_hat <- roots
  }

  ## Compute core matrices ##
  mats <- compute_matrices(geex_list   = geex_list,
                           theta       = theta_hat,
                           numDeriv_options = numDeriv_options,
                           ...)

  ## Compute corrections ##
  if(!is.null(corrections_list)){
    corrections <- make_corrections(mats, corrections_list)
  }

  ## Compute covariance estimate(s) ##
  Sigma_hat <- compute_sigma(A = mats$A, B = mats$B)

  list(parameters = theta_hat, vcov = Sigma_hat, corrections = corrections)
}


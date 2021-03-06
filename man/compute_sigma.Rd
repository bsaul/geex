% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_funs.R
\name{compute_sigma}
\alias{compute_sigma}
\title{Compute empirical sandwich covariate estimator}
\usage{
compute_sigma(A, B, solver = solve)
}
\arguments{
\item{A}{a matrix, generally the \code{.A} slot in a
\code{\linkS4class{sandwich_components}} object created in
\code{\link{estimate_sandwich_matrices}}}

\item{B}{a matrix, generally the \code{.B} slot in a
\code{\linkS4class{sandwich_components}} object created in
\code{\link{estimate_sandwich_matrices}}}

\item{solver}{the function used to compute the inverse of \code{A}, Defaults
to \code{\link{solve}}}
}
\value{
the \code{matrix} \code{Ainv \%*\% B \%*\% t(Ainv)}
}
\description{
Computes \eqn{\Sigma = A^{-1} B (A^{-1})^T }{\Sigma = A^{-1} B (A^{-1})^T} with
provided \eqn{A}{A} and \eqn{B}{B} matrices.
}
\examples{
A <- diag(2, nrow = 2, ncol = 2)
B <- matrix(4, nrow = 2, ncol = 2)
compute_sigma(A = A, B = B)
}

#' geex: M-estimation API
#'
#' geex provides an extensible API for estimating parameters and their covariance
#' from a set of estimating functions (M-estimation). M-estimation theory has a
#' long history [see reference in the M-estimation bibliography:
#' \url{https://bsaul.github.io/geex/articles/articles/mestimation_bib.html}.
#' For an excellent introduction, see the primer by L.A. Stefanski and D.D. Boos,
#' "The Calculus of M-estimation" (The American Statistician (2002), 56(1), 29-38)
#' (\url{http://www.jstor.org/stable/3087324}).
#'
#' M-estimation encompasses a broad swath of statistical estimators and ideas including:
#'
#' \itemize{
#' \item the empirical "sandwich" variance estimator
#' \item generalized estimating equations (GEE)
#' \item many maximum likelihood estimators
#' \item robust regression
#' \item and many more}
#'
#' geex can implement all of these using a user-defined estimating function.
#'
#' To learn more about geex, see the package vignettes: \code{browseVignettes(package = 'geex')}.
#'
#' @section Goals:
#' If you can specify a set of unbiased estimating equations, geex does the rest.
#' The goals of geex are simply:
#'
#' \itemize{
#' \item To minimize the translational distance between a set of estimating
#' functions and R code;
#' \item To return numerically accurate point and covariance estimates from
#' a set of unbiased estimating functions.
#' }
#'
#' geex does not, by itself, necessarily aim to be fast nor precise. Such goals
#' are left to the user to implement or confirm.
#'
"_PACKAGE"

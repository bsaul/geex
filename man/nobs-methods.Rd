% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{nobs,geex-method}
\alias{nobs,geex-method}
\alias{nobs,geex,geex-method}
\alias{nobs,geex_summary-method}
\alias{nobs,geex_summary,geex_summary-method}
\title{Extract the number observations}
\usage{
\S4method{nobs}{geex}(object)

\S4method{nobs}{geex_summary}(object)
}
\arguments{
\item{object}{a \code{\linkS4class{geex}} object}
}
\description{
Extract the number observations
}
\examples{
\dontrun{
library(geepack)
data('ohio')

glmfit  <- glm(resp ~ age, data = ohio,
              family = binomial(link = "logit"))
z  <- m_estimate(
  estFUN = example_ee,
  data = ohio,
  compute_roots = FALSE,
  units = 'id',
  roots = coef(glmfit),
  outer_args = list(model = glmfit))

nobs(z)
}
}

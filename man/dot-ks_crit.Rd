% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/halting_tests.R
\name{.ks_crit}
\alias{.ks_crit}
\title{Criterion function for ks_halt.}
\usage{
.ks_crit(covariate, condition)
}
\arguments{
\item{covariate}{A vector containing a covariate to match the conditions on.}

\item{condition}{A factor vector containing condition labels.}
}
\value{
The p-value.
}
\description{
Warnings such as "cannot compute exact p-value with ties" are suppressed.
}
\keyword{internal}

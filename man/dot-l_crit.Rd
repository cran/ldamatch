% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/halting_tests.R
\name{.l_crit}
\alias{.l_crit}
\title{Criterion function for l_halt.}
\usage{
.l_crit(covariate, condition)
}
\arguments{
\item{covariate}{A vector containing a covariate to match the conditions on.}

\item{condition}{A factor vector containing condition labels.}
}
\value{
The p-value.
}
\description{
Warnings such as "ANOVA F-tests on an essentially perfect fit are unreliable"
are suppressed.
}
\keyword{internal}

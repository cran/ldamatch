% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.calc_p_thresh_ratio}
\alias{.calc_p_thresh_ratio}
\title{Calculates p-value-threshold ratio.}
\usage{
.calc_p_thresh_ratio(
  condition,
  covariates,
  halting_test,
  thresh,
  silent = TRUE
)
}
\arguments{
\item{condition}{A factor vector containing condition labels.}

\item{covariates}{A columnwise matrix containing
covariates to match the conditions on.}

\item{halting_test}{A function to apply to `covariates` (in matrix form)
which is TRUE iff the conditions are matched.
Signature: halting_test(condition, covariates, thresh).
The following halting tests are part of this package:
\code{\link{t_halt}}, \code{\link{U_halt}},
\code{\link{l_halt}}, \code{\link{ad_halt}},
\code{\link{ks_halt}}, \code{\link{wilks_halt}},
\code{\link{f_halt}}.
You can create the intersection of two or more halting
tests using \code{\link{create_halting_test}}.}

\item{thresh}{The return value of halting_test has to be greater than
or equal to thresh for the matched groups.}

\item{silent}{If FALSE, prints warning when the test statistic cannot be
calculated; if TRUE (the default) they are not printed.}
}
\value{
The p-value-threshold ratio, or NA if the p-value could not be
        calculated.

The p-value / thresh ratio.
}
\description{
Calculates p-value-threshold ratio.
}
\keyword{internal}

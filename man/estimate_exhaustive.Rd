% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/estimate_exhaustive.R
\name{estimate_exhaustive}
\alias{estimate_exhaustive}
\title{Estimates the maximum number of cases to be checked during exhaustive search.}
\usage{
estimate_exhaustive(min_preserved, condition, cases_per_second = 100,
  print_info = get("PRINT_INFO", .ldamatch_globals))
}
\arguments{
\item{min_preserved}{Assumes that at least a total of this many subjects
                         will be preserved.}

\item{condition}{A factor vector containing condition labels.}

\item{cases_per_second}{Assumes that this number of cases are checked out
                         per second, for estimating the time it takes to run
                         the exhaustive search; default: 100.}

\item{print_info}{If TRUE, prints partial calculations as well for
                         the number of cases and estimated time when removing
                         1, 2, ... subjects.}
}
\value{
The maximum number of cases.
}
\description{
Estimates the maximum number of cases to be checked during exhaustive search.
}
\examples{
estimate_exhaustive(58, as.factor(c(rep('ALN', 25), rep('TD', 44))))
estimate_exhaustive(84, as.factor(c(rep('ASD', 51), rep('TD', 44))))
}

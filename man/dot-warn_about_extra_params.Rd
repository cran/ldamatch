% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.warn_about_extra_params}
\alias{.warn_about_extra_params}
\title{Warns about extra (i.e. unused) parameters.}
\usage{
.warn_about_extra_params(given_args = NULL, ...)
}
\arguments{
\item{given_args}{The names of arguments given to the search function.}

\item{...}{Consumes extra parameters that are not used by the
search algorithm at hand; this function gives a warning
about the ones whose value is not NULL that their value
is not used.}
}
\description{
Warns about extra (i.e. unused) parameters.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_exhaustive.R
\name{.create_Cartesian_iterable}
\alias{.create_Cartesian_iterable}
\title{Creates Cartesian product of iterators.}
\usage{
.create_Cartesian_iterable(initializers, get_next, sspace)
}
\arguments{
\item{initializers}{A list of initializer functions (with no arguments)
for iterators.}

\item{get_next}{A function for retrieving next item for an iterator
argument; it assumes that the iterator returns NULL
when finished.}

\item{sspace}{elements to be used (a list of vectors)}
}
\value{
A function that returns list of values, and stops with
       "StopIteration" message when finished, so that it can be used with
       the iterators::iter() function to create an iterator that works with
       foreach.
}
\description{
Creates Cartesian product of iterators.
}
\keyword{internal}

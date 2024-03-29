% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_heuristic3.R
\name{search_heuristic3}
\alias{search_heuristic3}
\title{Finds matching using depth-first search, looking ahead n steps.}
\usage{
search_heuristic3(
  condition,
  covariates,
  halting_test,
  thresh,
  props,
  max_removed_per_cond,
  tiebreaker = NULL,
  min_preserved = length(levels(condition)),
  lookahead = 2,
  prefer_test = TRUE,
  print_info = TRUE,
  max_removed_per_step = 1,
  max_removed_percent_per_step = 0.5,
  ratio_for_slowdown = 0.5,
  given_args = NULL,
  ...
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

\item{props}{Either the desired proportions (percentage) of the
sample for each condition as a named vector,
or the names of the conditions
for which we prefer to preserve the subjects,
in decreasing order of preference. If not specified, the
(full) sample proportions are used.
This is preferred among configurations with the same
taken into account by the other methods to some extent.
For example, c(A = 0.4, B = 0.4, C = 0.2) means that
we would like the number of subjects in groups A, B, and
C to be around 40\%, 40\%, and 20\% of the total number of
subjects, respectively. Whereas c("A", "B", "C") means
that if possible, we would like to keep all subjects
in group A, and prefer keeping subjects in B, even if
it results in losing more subjects from C.}

\item{max_removed_per_cond}{The maximum number of subjects that can be removed from
each group. It must have a valid number for each group.}

\item{tiebreaker}{NULL, or a function similar to halting_test, used to
decide between cases for which halting_test yields
equal values.}

\item{min_preserved}{The minimum number of preserved subjects.
It can be used to ensure that the search will not take
forever to run, but instead fail when a solution is not
found when preserving this number of subjects.}

\item{lookahead}{The lookahead to use: a positive integer.
It is used by the heuristic3 and heuristic4 algorithms,
with a default of 2.
The running time is O(N ^ lookahead), wheren N is the
number of subjects.}

\item{prefer_test}{If TRUE, prefers higher test statistic more than
the expected group size proportion; default is TRUE.
Used by all algorithms except exhaustive, which always}

\item{print_info}{If TRUE, prints summary information on the input and the
results, as well as progress information for the
exhaustive search and random algorithms. Default: TRUE;
can be changed using
\code{\link{set_param}("PRINT_INFO", FALSE)}.}

\item{max_removed_per_step}{The number of equivalent subjects
that can be removed in each step. (The actual allowed
number may be less depending on the p-value / theshold ratio.)
This parameters is used by the heuristic3 and heuristic4
algorithms, with a default value of 1.}

\item{max_removed_percent_per_step}{The percentage of remaining subjects
that can be removed in each step.
Used when max_removed_per_step > 1,
with a default value of 0.5.}

\item{ratio_for_slowdown}{The p-value / threshold ratio at which
it starts removing subjects one by one.
Used when max_removed_per_step > 1,
with a default value of 0.5.}

\item{given_args}{The names of arguments given to the search function.}

\item{...}{Consumes extra parameters that are not used by the
search algorithm at hand; this function gives a warning
about the ones whose value is not NULL that their value
is not used.}
}
\value{
All results found by search method in a list. It raises a
        "Convergence failure" error if it cannot find a matched set.
}
\description{
In each step, it removes one subject from the set of subjects with the
smallest associated p-value after "lookahead" steps.
}
\details{
Note that this algorithm is not deterministic, as it chooses one possible
path randomly when there are multiple apparently equivalent ones. In practice
this means that it may return different results on different runs (including
the case that it fails to converge to a solution in one run, but converges in
another run). If print_info = TRUE (the default), you will see a message
about "Random choices" if the algorithm needed to make random path choices.
}

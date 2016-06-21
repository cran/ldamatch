#' Returns human readable format for number of seconds.
#'
#' @param seconds       The number of seconds to convert to human-readable form.
#'
#' @param num_decimals  The number of decimals to print in the output.
#'
#' @return A string containing "<number> seconds/minutes/hours/days/years".
#'
.get_human_readable <- function(seconds, num_decimals = 1) {
    value <- seconds
    measurement <- "second"
    for (convert in list(
        list(60, "minute"),
        list(60, "hour"),
        list(24, "day"),
        list(7, "week"),
        list(365 / (12 * 7), "month"),
        list(12, "year")
    )) {
        if (round(value, num_decimals) >= convert[[1]]) {
            value <- value / convert[[1]]
            measurement <- convert[[2]]
        } else {
            break
        }
    }
    value <- round(value, num_decimals)
    paste(value, if (value == 1)
        measurement
        else
            paste0(measurement, "s"))
}


#' Estimates the maximum number of cases to be checked during exhaustive search.
#'
#' @param min_preserved     Assumes that at least a total of this many subjects
#'                          will be preserved.
#'
#' @param cases_per_second  Assumes that this number of cases are checked out
#'                          per second, for estimating the time it takes to run
#'                          the exhaustive search; default: 100.
#'
#' @param print_info        If TRUE, prints partial calculations as well for
#'                          the number of cases and estimated time when removing
#'                          1, 2, ... subjects.
#'
#' @param group_sizes       A particular set of group sizes that we know a
#'                          matched solution for; min_preserved need not be
#'                          specified if this one is.
#'
#' @param props             The desired proportions (percentage) of the sample
#'                          for each condition; if this and group_sizes are both
#'                          specified, the maximum number of cases to considered
#'                          by the exhaustive search can be calculated more
#'                          precisely.
#'
#' @inheritParams match_groups
#'
#' @return The maximum number of cases: an integer if not greater than the
#' maximum integer size (.Machine$integer.max), otherwise a Big Integer
#' (see the gmp package).
#'
#' @examples
#' estimate_exhaustive(58, as.factor(c(rep("ALN", 25), rep("TD", 44))))
#' estimate_exhaustive(84, as.factor(c(rep("ASD", 51), rep("TD", 44))))
#'
#' @import data.table
#' @import gmp
#' @importFrom iterpc iterpc getlength
#'
#' @export
estimate_exhaustive <-
    function(min_preserved = sum(group_sizes),
             condition,
             cases_per_second = 100,
             print_info = TRUE,
             max_removed = NULL,
             group_sizes = NULL,
             props = NULL) {
        divergence <- NULL  # to suppress codetools warnings
        stopifnot(is.factor(condition), length(cases_per_second) == 1)
        min_preserved <- max(min(min_preserved, length(condition)),
                             length(levels(condition)))
        sspace <- split(seq_along(condition), condition)
        max_removed <- .normalize_max_removed(max_removed, condition)
        minpergrp <- vapply(sspace, length, 0) - max_removed
        grpnames <- names(sspace)
        num_cases <- 0
        grpsizes <- data.table::data.table(t(vapply(sspace, length, 0)))
        if (!is.null(group_sizes)) {
            RUnit::checkTrue(
                setequal(names(group_sizes), names(grpsizes)),
                "The names of group_sizes and the conditions must be the same"
            )
            group_sizes <- group_sizes[names(grpsizes)]
            if (!is.null(props))
                props <- .normalize_props(props, condition)
        }
        while (min_preserved < sum(grpsizes[1, grpnames, with = FALSE])) {
            grpsizes <- .decrease_group_sizes(grpsizes, grpnames, minpergrp)
            if (nrow(grpsizes) == 0)
                break
            if (!is.null(props)) {
                grpsizes <- .sort_group_sizes(grpsizes, grpnames, props)
                grpsizes[, divergence := NULL]
            }
            for (grpsizes_row in 1:nrow(grpsizes)) {
                inc(num_cases) <- do.call(prod,
                                          lapply(grpnames, function(cond)
                                              iterpc::getlength(iterpc::iterpc(length(
                                                  sspace[[cond]]
                                              ), grpsizes[[grpsizes_row, cond]]),
                                              bigz = TRUE)))
                if (!is.null(props) &&
                    isTRUE(all.equal(unlist(grpsizes[grpsizes_row]), group_sizes)))
                    break
            }
            if (print_info) {
                cat(
                    "If",
                    sum(grpsizes[1, grpnames, with = FALSE]),
                    "of",
                    length(condition),
                    "kept: at most",
                    as.character(num_cases),
                    "cases. "
                )
                cat(
                    "If ",
                    cases_per_second,
                    " cases per second evaluated: ",
                    .get_human_readable(as.double(num_cases / cases_per_second)),
                    ".\n",
                    sep = ""
                )
            }
        }
        if (num_cases > .Machine$integer.max)
            num_cases
        else
            as.integer(num_cases)
    }

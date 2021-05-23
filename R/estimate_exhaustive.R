#' Returns human readable format for number of seconds.
#'
#' @param seconds       The number of seconds to convert to human-readable form.
#'
#' @param num_decimals  The number of decimals to print in the output.
#'
#' @return A string containing "<number> seconds/minutes/hours/days/years".
#'
#' @keywords internal
.get_human_readable <- function(seconds, num_decimals = 1) {
    value <- seconds
    measurement <- "second"
    for (convert in list(
        # multiplier from previous measurement; name [and possibly plural name]; max. decimals
        list(60, "minute"),
        list(60, "hour"),
        list(24, "day"),
        list(7, "week"),
        list(365 / (12 * 7), "month"),
        list(12, "year"),
        list(100, c("century", "centuries"), 0)
    )) {
        num_decimals_for_value <-
            min(num_decimals, if (length(convert) >= 3)
                convert[[3]]
                else
                    NULL)
        if (round(value, num_decimals_for_value) >= convert[[1]]) {
            value <- value / convert[[1]]
            measurement <- convert[[2]]
        } else {
            break
        }
    }
    rounded_value <- round(value, num_decimals_for_value)
    paste(
        rounded_value,
        if (rounded_value == 1)
            measurement[[1]]
        else if (length(measurement) >= 2)
            measurement[[2]]
        else
            paste0(measurement[[1]], "s")
    )
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
#' @param props             Either the desired proportions (percentage) of the
#'                          sample for each condition as a named vector,
#'                          or the names of the conditions
#'                          for which we prefer to preserve the subjects,
#'                          in decreasing order of preference. If not specified, the
#'                          (full) sample proportions are used.
#'
#' @param max_cases         Once it is certain that the number of cases is
#'                          definitely above this number, calculation stops. In this case,
#'                          the returned number is guaranteed to be larger than max_cases,
#'                          but it is not the exact number of exhaustive cases.
#'                          Default is infimum, i.e. the exact number of cases is calculated.
#'
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
             max_removed_per_cond = NULL,
             group_sizes = NULL,
             props = prop.table(table(condition)),
             max_cases = Inf) {
        divergence <- NULL  # to suppress codetools warnings
        stopifnot(is.factor(condition), length(cases_per_second) == 1)
        stopifnot(is.logical(print_info) && length(print_info) == 1)
        if (!is.null(props))
            props <- .normalize_props(props, condition, keep_last_item = TRUE)
        condition_names <- intersect(levels(condition), unique(c(
            as.character(condition),
            names(max_removed_per_cond),
            names(group_sizes[group_sizes != 0]),
            names(props[props != 0.0])
        )))
        condition <- factor(as.character(condition),
                    intersect(levels(condition), condition_names))
        props <- props[condition_names]
        min_preserved <- max(min(min_preserved, length(condition)),
                             length(levels(condition)))
        sspace <- split(seq_along(condition), condition)
        max_removed_per_cond <-
            .normalize_max_removed_per_cond(max_removed_per_cond, condition)
        minpergrp <-
            vapply(sspace, length, 0) - max_removed_per_cond
        grpnames <- names(sspace)
        grpsizes <-
            data.table::data.table(t(vapply(sspace, length, 0)))
        if (!is.null(group_sizes)) {
            RUnit::checkTrue(
                all(group_sizes[setdiff(names(group_sizes), condition_names)] == 0, na.rm = TRUE),
                sprintf("The names of group_sizes and the conditions must be the same: %s; %s",
                        paste(names(group_sizes), collapse = ', '),
                        paste(names(grpsizes), collapse = ', '))
            )
            group_sizes <- group_sizes[condition_names]
        }

        num_cases <- 0
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
                    isTRUE(all.equal(
                        unlist(grpsizes[grpsizes_row]),
                        group_sizes,
                        tolerance = .tolerance
                    )))
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
            if (num_cases > max_cases) {
                cat(
                    "Number of cases is definitely more than max_cases argument (",
                    max_cases,
                    "), stopping calculation.\n",
                    sep = ""
                )
                break
            }
        }
        if (num_cases > .Machine$integer.max)
            num_cases
        else
            as.integer(num_cases)
    }

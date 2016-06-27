#' Calculates multipliers used in search_random.
#'
#' Derives multiplier for rcounts (the number of subjects that can be removed)
#' such that the proportion of the expected sizes of groups will be props.
#' The returned multipliers will be in the range of 0 to 1.
#'
#' @param counts   The number of subjects for each group.
#'
#' @param rcounts  The number of subjects that can be removed for each group.
#'
#' @param props The expected proportion of subjects for each group.
.calc_multipliers <- function(counts, rcounts, props) {
    m <-
        (counts - props / props[[1]] * (counts[[1]] - rcounts[[1]])) / rcounts
    m <- vapply(m / max(m), max, numeric(1), 0.0)
    m
}


#' Searches by randomly selecting subspaces with decreasing expected size.
#'
#' @param max_removed   The maximum number of subjects that can be removed from
#'                      each group. It must have a valid number for each group,
#'                      and the groups must be in the same order as in sspace.
#'
#' @inheritParams match_groups
#' @inheritParams .warn_about_extra_params
#'
#' @return All results found by search method in a list. It raises a
#          "Convergence failure" error if it cannot find a matched set.
#'
#' @importFrom RUnit checkTrue
#' @importFrom stats rbinom
search_random <- function(condition,
                          covariates,
                          halting_test,
                          thresh,
                          props,
                          max_removed,
                          tiebreaker = NULL,
                          replicates,
                          print_info = TRUE,
                          ...) {
    .warn_about_extra_params(...)
    # Checks replicates argument.
    if (is.null(replicates)) {
        replicates <- get("RND_DEFAULT_REPLICATES", .ldamatch_globals)
    } else {
        RUnit::checkTrue(
            length(replicates) == 1 && replicates %% 1 == 0,
            "The replicates parameter must be one integer number"
        )
    }
    # Searches subject space.
    sspace <- split(seq_along(condition), condition)
    counts <- table(condition)  # total number of subjects
    rcounts <-
        vapply(sspace, length, 0)  # number of subjects for removal
    rcounts <-
        rcounts - (rcounts == counts)  # do not remove all subjects
    multipliers <- .calc_multipliers(counts, rcounts, props)
    start <-
        max(1, floor(min((replicates - 1) / (rcounts * multipliers - 1)
        )))
    end <- replicates + start - 1
    best <- NULL
    best_num <- 0
    for (i in start:end) {
        # number of subjects to remove per group
        if (!best_num)
            pos <- i
        nrs <-
            stats::rbinom(length(sspace), rcounts, multipliers * pos / end)
        nrs <- pmin(nrs, max_removed)
        is.in <- rep(TRUE, length(condition))
        mapply(function(s, nr, len)
            is.in[sample(s, nr)] <<- FALSE, sspace, nrs)
        if (sum(is.in) < best_num)
            next
        # Tests and returns binary vector if anything is found.
        ratio <-
            halting_test(condition[is.in], covariates[is.in, , drop = FALSE], thresh)
        if (!ratio) {
            next
        } else if (!best_num) {
            cmp <- 1
            best <- list(metric = ratio, sets = list(is.in))
            best_num <- sum(is.in)
        } else {
            cmp <- compare_ldamatch_outputs(
                is.in,
                best$sets[[1]],
                condition,
                covariates,
                halting_test,
                props,
                tiebreaker = tiebreaker
            )
            if (cmp > 0) {
                best <- list(metric = ratio, sets = list(is.in))
                best_num <- sum(is.in)
            } else if (cmp == 0 &&
                       !any(sapply(best$sets, identical, is.in))) {
                best$sets <- c(best$sets, list(is.in))
            }
        }
        if (cmp >= 0 && print_info)
            cat(
                "Found matching: ",
                paste(
                    levels(condition),
                    table(condition[is.in]),
                    sep = ": ",
                    collapse = "; "
                ),
                " (total: ",
                best_num,
                ")\n",
                sep = ""
            )
    }
    if (!best_num)
        stop("Convergence failure")
    best$sets
}

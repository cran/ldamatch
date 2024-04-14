#' Chooses rows with best test statistic.
#'
#' @param dat  A data.table with an ind column with indices for items to
#'             consider dropping.
#'
#' @inheritParams match_groups
#'
#' @return A data.table containing only the rows of dat with for best test
#'         statistic values (decided primarily by halting_test, then by
#'         tiebreaker).
#'
#' @keywords internal
.choose_best_test_statistic <- function(dat,
                                        condition,
                                        covariates,
                                        halting_test,
                                        thresh,
                                        tiebreaker) {
    ratio <-
        tiebreaker_ratio <- ind <- NULL  # to suppress codetools warnings
    dat[, ratio :=  vapply(ind, function(ind)
        .calc_p_thresh_ratio(condition[-ind], covariates[-ind, , drop = FALSE],
                             halting_test, thresh), 0.0)]
    dat <- suppressWarnings(dat[ratio == max(ratio, na.rm = TRUE)])
    if (nrow(dat) == 0)
        stop("Convergence failure")
    if (nrow(dat) > 1 && is.function(tiebreaker)) {
        dat[, tiebreaker_ratio :=  vapply(ind, function(ind)
            .calc_p_thresh_ratio(condition[-ind], covariates[-ind, , drop = FALSE],
                                 tiebreaker, thresh), 0.0)]
        dat_tiebreaker <-
            suppressWarnings(dat[tiebreaker_ratio == max(tiebreaker_ratio, na.rm = TRUE)])
        if (nrow(dat_tiebreaker) >= 1)
            dat <- dat_tiebreaker
    }
    dat
}


#' OBSOLETE: Finds matching using depth-first search recursively.
#'
#' Please use the heuristic3 search algorithm with lookahead=1 instead
#' for nearly equivalent results. Note that heuristic3 is parallelized,
#' more memory efficient, and chooses subject to remove randomly from among
#' equivalent choices instead of choosing the first one deterministically.
#' This function is implemented recursively, so may run out of memory when
#' applied to many subjects.
#'
#' In each step, it removes one subject from the set of subjects with
#' the smallest p-value recursively.
#'
#' @param max_removed_per_cond   The maximum number of subjects that can be removed from
#'                      each group. It must have a valid number for each group.
#'
#' @inheritParams match_groups
#' @inheritParams .warn_about_extra_params
#'
#' @return All results found by search method in a list. It raises a
#'         "Convergence failure" error if it cannot find a matched set.
#'
#' @import data.table
search_heuristic2 <- function(condition,
                              covariates,
                              halting_test,
                              thresh,
                              props,
                              max_removed_per_cond,
                              tiebreaker = NULL,
                              prefer_test = TRUE,
                              print_info = TRUE,
                              given_args = NULL,
                              ...) {
    .warn_about_extra_params(given_args, ...)
    if (print_info) {
        cat(sprintf("Number of subjects: %d\n", length(condition)))
    }
    # Ends recursion when matched.
    if (halting_test(condition, covariates, thresh))
        return(list(rep(TRUE, length(condition))))
    # Finds best direction.
    dat <-
        data.table::data.table(ind = which(condition %in% names(which(
            max_removed_per_cond > 0
        ))))
    if (prefer_test)
        dat <- .choose_best_test_statistic(dat,
                                           condition,
                                           covariates,
                                           halting_test,
                                           thresh,
                                           tiebreaker)
    divergence <- NULL
    if (nrow(dat) > 1) {
        dat[, divergence :=  vapply(ind, function(ind)
            .calc_subject_balance_divergence(table(condition[-ind]), props), 0.0)]
        dat <-
            suppressWarnings(dat[divergence == max(divergence, na.rm = TRUE)])
        if (nrow(dat) == 0)
            stop("Convergence failure")
    }
    if (!prefer_test && nrow(dat) > 1)
        dat <- .choose_best_test_statistic(dat,
                                           condition,
                                           covariates,
                                           halting_test,
                                           thresh,
                                           tiebreaker)
    # Removes each subject for best p-value / thresh ratio in turn, and finds
    # balance recursively.
    for (ind in dat$ind) {
        max_removed_per_cond.for_subset <- max_removed_per_cond
        dec(max_removed_per_cond.for_subset[condition[ind]]) <- 1
        lis.in.for_subset  <- try(search_heuristic2(
            condition[-ind],
            covariates[-ind, , drop = FALSE],
            halting_test,
            thresh,
            props,
            max_removed_per_cond.for_subset,
            tiebreaker
        ),
        silent = TRUE)
        if (inherits(lis.in.for_subset, "try-error")) {
            is.in <- logical(length(condition))
            is.in[-ind] <- lis.in.for_subset[[1]]
            is.in[[ind]] <- FALSE
            return(list(is.in))
        }
        break  # we do not want it to take a lot of time
    }
    stop("Convergence failure")
}

#' Returns the subjects whose inclusion status changed in the highest number of
#' cases.
#'
#' @param best_sets   A list of logical vectors showing which subjects stay in.
#'                    (E.g., the return value from .choose_best_subjects().)
#'
#' @param is.in       A logical vector showing which items are preserved
#'                    currently; versions resulting by changing indices for
#'                    each candidate are then compared.
#'
#' @return A list of subject indices corresponding to is.in.
.get_most_frequently_excluded_subjects <-
    function(best_sets, is.in) {
        excluded_subject_counts <- table(unlist(lapply(
            best_sets, .get_difference_inds, is.in
        )))
        best_subjects <-
            as.integer(names(excluded_subject_counts)[excluded_subject_counts == max(excluded_subject_counts)])
        best_subjects
    }

#' Finds matching using depth-first search, looking ahead n steps.
#'
#' In each step, it removes one subject from the set of subjects that were
#' removed on most paths after "lookahead" steps, preferring one with the
#' smallest associate p-value.
#'
#' Note that this algorithm is not deterministic, as it chooses one possible
#' subject for removal randomly when there are multiple apparently equivalent ones.
#' In practice it means that it may return different results on different runs
#' (including the case that it fails to converge to a solution in one run,
#' but converges in another run). If print_info = TRUE (the default), you will
#' see a message about "Random choices" if the algorithm needed to make such
#' random decisions.
#'
#' @param max_removed   The maximum number of subjects that can be removed from
#'                      each group. It must have a valid number for each group.
#'
#' @inheritParams match_groups
#' @inheritParams .warn_about_extra_params
#'
#' @return All results found by search method in a list. It raises a
#'         "Convergence failure" error if it cannot find a matched set.
#'
#' @import foreach
#' @importFrom iterpc iterpc iter_wrapper
search_heuristic4 <- function(condition,
                              covariates,
                              halting_test,
                              thresh,
                              props,
                              max_removed,
                              tiebreaker = NULL,
                              min_preserved = NULL,
                              lookahead = NULL,
                              print_info = TRUE,
                              ...) {
    .warn_about_extra_params(...)
    if (is.null(lookahead))
        lookahead <- 2
    else
        stopifnot(lookahead >= 1)
    if (is.null(min_preserved))
        min_preserved <- length(levels(condition))
    is.in <- rep(TRUE, length(condition))
    look <-
        1  # level we are looking at (the number of subjects dropped)
    random_choices <- list()
    step <- 0
    while (sum(is.in) > min_preserved) {
        if (print_info)
            cat(sprintf("Number of subjects: %d\n", sum(is.in)))
        can_be_removed <-
            condition %in% names(which(max_removed > 0))
        repeat {
            if (print_info)
                cat("Lookahead", look)
            is.in_candidates <- is.in & can_be_removed
            candidates <- iterpc::iter_wrapper(iterpc::iterpc(
                sum(is.in_candidates),
                look,
                which(is.in_candidates)
            ))
            best_sets <- .choose_best_subjects(
                candidates,
                is.in,
                condition,
                covariates,
                halting_test,
                thresh,
                tiebreaker,
                props
            )
            if (length(best_sets) == 0)
                stop("Convergence failure")
            if (halting_test(condition[best_sets[[1]]],
                             covariates[best_sets[[1]], , drop = FALSE], thresh)) {
                if (print_info) {
                    cat(
                        sprintf(
                            "\nFound solution at %d: %s\n",
                            look,
                            .vector_list_to_string(
                                lapply(best_sets, .get_difference_inds, is.in)
                            )
                        )
                    )
                    if (length(random_choices) > 0)
                        cat(
                            "Random choices:",
                            .vector_list_to_string(random_choices, " "),
                            "\n"
                        )
                }
                return(best_sets)  # returns list of best results
            }
            if (print_info)
                cat(", best sets:",
                    .vector_list_to_string(lapply(
                        best_sets, .get_difference_inds, is.in
                    )),
                    "\n")
            if (look >= lookahead)
                break
            look <- look + 1
        }
        # choose one subject for removal that was among the most frequently
        # removed ones; among those, prefer the one with the lowest p-value,
        # then choose randomly among the remaining ones
        inc(step) <- 1
        best_subjects <-
            .get_most_frequently_excluded_subjects(best_sets, is.in)
        if (length(best_subjects) > 1)
            best_subjects <- .get_most_frequently_excluded_subjects(
                .choose_best_subjects(
                    best_subjects,
                    is.in,
                    condition,
                    covariates,
                    halting_test,
                    thresh,
                    tiebreaker,
                    props
                ),
                is.in
            )
        if (length(best_subjects) > 1) {
            chosen_ind <- sample(best_subjects, 1)
            random_choices <- c(random_choices, list(
                list(
                    step = step,
                    num_choices = length(best_subjects),
                    chosen_ind = chosen_ind
                )
            ))
        } else {
            chosen_ind <- best_subjects
        }
        is.in <- .flip_ind(is.in, chosen_ind)
        dec(max_removed[condition[chosen_ind]]) <- 1
    }
    if (print_info && length(random_choices) > 0)
        cat("Random choices:",
            .vector_list_to_string(random_choices, "; "),
            "\n")
    stop("Convergence failure")
}

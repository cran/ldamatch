#' Chooses best one(s) of a set of subjects having the best p-value(s).
#'
#' Used as first parameter of .search_heuristic_with_lookahead.
#' It traces the best parents for a set of subjects to be removed back to the first level.
#' Note that the subject indices may not be unique: ones that occur in more configurations may be listed multiple times.
#'
#' @return The table of counts for the chosen indices within is.in.
#'
#' @importFrom utils combn
#'
#' @keywords internal
.choose_subject_with_best_p_value_from_subject_tuples <-
    function(is.in,
             best_sets,
             look,
             condition,
             covariates,
             halting_test,
             thresh,
             tiebreaker,
             props,
             prefer_test,
             max_removed_per_cond,
             max_removed_in_next_step,
             ratio_for_slowdown,
             remove_best_only,
             print_info) {
        step_look <- look
        repeat {
            best_subject_counts <- table(unlist(best_sets))
            if (step_look == 1 ||
                length(best_subject_counts) <= max_removed_in_next_step)
                return(best_subject_counts)
            dec(step_look) <- 1
            parent_inds <-
                .unique_list(unlist(lapply(best_sets, function(ind)
                    utils::combn(ind, step_look, simplify = FALSE)), recursive = FALSE))
            if (print_info)
                cat(
                    sprintf(
                        "Looking for parent at %d, parent inds: %s\n",
                        step_look,
                        .vector_list_to_string(parent_inds)
                    )
                )
            best_sets <- .choose_best_subjects(
                parent_inds,
                is.in,
                condition,
                covariates,
                halting_test,
                thresh,
                tiebreaker,
                props,
                prefer_test,
                max_removed_per_cond,
                max_removed_in_next_step,
                ratio_for_slowdown,
                remove_best_only
            )
            if (print_info) {
                cat(
                    sprintf(
                        "Looking for parent at %d, best parent: %s\n",
                        step_look,
                        .vector_list_to_string(best_sets)
                    )
                )
            }
        }
    }


#' Finds matching using depth-first search, looking ahead n steps.
#'
#' In each step, it removes one subject from the set of subjects with the
#' smallest associated p-value after "lookahead" steps.
#'
#' Note that this algorithm is not deterministic, as it chooses one possible
#' path randomly when there are multiple apparently equivalent ones. In practice
#' this means that it may return different results on different runs (including
#' the case that it fails to converge to a solution in one run, but converges in
#' another run). If print_info = TRUE (the default), you will see a message
#' about "Random choices" if the algorithm needed to make random path choices.
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
#' @importFrom utils combn
search_heuristic3 <- function(condition,
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
                              ...) {
    .search_heuristic_with_lookahead(
        choose_from_subject_tuples = .choose_subject_with_best_p_value_from_subject_tuples,
        condition = condition,
        covariates = covariates,
        halting_test = halting_test,
        thresh = thresh,
        props = props,
        max_removed_per_cond = max_removed_per_cond,
        tiebreaker = tiebreaker,
        min_preserved = min_preserved,
        lookahead = lookahead,
        prefer_test = prefer_test,
        print_info = print_info,
        max_removed_per_step = max_removed_per_step,
        max_removed_percent_per_step = max_removed_percent_per_step,
        ratio_for_slowdown = ratio_for_slowdown,
        given_args = given_args,
        ...
    )
}

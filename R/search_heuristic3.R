#' Flips logical vector at specified indices
#'
#' @param is.in        A logical vector showing which items are preserved.
#'
#' @param ind  Integer indices for the is.in logical vector.
#'
#' @return A logical vector identical to is.in except for indices in ind where
#'         it is is.in negated.
.flip_ind <- function(is.in, ind) {
    v <- as.vector(ind)
    is.in[v] <- !is.in[v]
    is.in
}


#' Creates string from list of vectors.
#'
#' @param lv   A list of vectors.
#'
#' @param sep  A string to be inserted between the name of a vector item and
#' its value.
#'
#' @return A character string.
.vector_list_to_string <- function(lv, sep = "") {
    paste(sapply(lv, function(v)
        paste(
            names(v), v, sep = sep, collapse = ", "
        )),
        collapse = "; ")
}


#' Returns indices of differences from reference vector.
#'
#' @param v     A logical vector.
#' @param vref  The reference logical vector to which v is compared.
#'
#' @return An integer vector containing the indices where v differs from vref.
.get_difference_inds <- function(v, vref) {
    which(v != vref)
}


#' Chooses best set of subjects in a set.
#'
#' @param candidates  An iterator returning (or a list containing)
#'                    indices for the is.in logical vector whose in / out status
#'                    is to be changed.
#'
#' @param is.in       A logical vector showing which items are preserved
#'                    currently; versions resulting by changing indices for
#'                    each candidate are then compared.
#'
#' @inheritParams match_groups
#' @inheritParams .warn_about_extra_params
#'
#' @return A list containing the best is.in vectors resulting from changing
#'         is.in at the indices for each candidate.
#'
#' @import foreach
.choose_best_subjects <- function(candidates,
                                  is.in,
                                  condition,
                                  covariates,
                                  halting_test,
                                  thresh,
                                  tiebreaker,
                                  props) {
    ind <- NULL  # to suppress codetools warnings
    best <- foreach::foreach(
        ind = candidates,
        .init = list(metric = 0, sets = list()),
        .combine = .combine_sets,
        .inorder = FALSE
    ) %dopar% {
        cis.in <- .flip_ind(is.in, ind)
        ratio <- .calc_p_thresh_ratio(condition[cis.in], covariates[cis.in, , drop = FALSE],
                                      halting_test, thresh)
        list(metric = min(ratio, 1.0), set = cis.in)
    }
    if (length(best$sets) > 1 && is.function(tiebreaker)) {
        best <- foreach::foreach(
            cis.in = best$sets,
            .init = list(metric = 0, sets = list()),
            .combine = .combine_sets,
            .inorder = FALSE
        ) %dopar% {
            ratio <- .calc_p_thresh_ratio(condition[cis.in], covariates[cis.in, , drop = FALSE],
                                          tiebreaker, thresh)
            list(metric = min(ratio, 1.0), set = cis.in)
        }
    }
    if (length(best$sets) > 1) {
        best <- foreach::foreach(
            cis.in = best$sets,
            .init = list(metric = -Inf, sets = list()),
            .combine = .combine_sets,
            .inorder = FALSE
        ) %dopar% {
            divergence <- .calc_subject_balance_divergence(table(condition[cis.in]), props)
            list(metric = -1 - divergence, set = cis.in)
        }
    }
    best$sets
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
#' @param max_removed   The maximum number of subjects that can be removed from
#'                      each group. It must have a valid number for each group.
#'
#' @inheritParams match_groups
#' @inheritParams .warn_about_extra_params
#'
#' @return All results found by search method in a list. It raises a
#'         "Convergence failure" error if it cannot find a matched set.
#'
#' @importFrom iterpc iterpc iter_wrapper
#' @importFrom utils combn
search_heuristic3 <- function(condition,
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
        # tracing best parents back to current level
        step_look <- look
        while (step_look > 1) {
            step_look <- step_look - 1
            inds <- sapply(best_sets,
                           function(set)
                               which(is.in != set), simplify = FALSE)
            parent_inds <- .unique_list(unlist(
                sapply(inds, function(ind)
                    utils::combn(ind, step_look, simplify = FALSE),
                    simplify = FALSE),
                recursive = FALSE
            ))  # unique list items
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
                props
            )
            if (print_info) {
                cat(
                    sprintf(
                        "Looking for parent at %d, best parent: %s\n",
                        step_look,
                        .vector_list_to_string(
                            lapply(best_sets, .get_difference_inds, is.in)
                        )
                    )
                )
            }
        }
        # committing to best direction
        inc(step) <- 1
        num_choices <- length(best_sets)
        chosen_set <- sample(seq_len(num_choices), 1)
        if (num_choices > 1)
            random_choices <- c(random_choices, list(
                list(
                    step = step,
                    num_choices = num_choices,
                    chosen_set = chosen_set
                )
            ))
        old_is.in <- is.in
        is.in <- best_sets[[chosen_set]]
        ind <- which(old_is.in != is.in)
        dec(max_removed[condition[ind]]) <- 1
    }
    if (print_info && length(random_choices) > 0)
        cat("Random choices:",
            .vector_list_to_string(random_choices, "; "),
            "\n")
    stop("Convergence failure")
}

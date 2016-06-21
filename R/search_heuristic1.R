#' Finds matching using heuristic based on linear discriminant analysis.
#'
#' At each vertex of the search graph, this takes a step which moves the
#' proportions of conditions in the subspace closer to the desired (or sample)
#' proportions, so the expected proportions are enforced.
#'
#' @param max_removed   The maximum number of subjects that can be removed from
#'   each group. It must have a valid number for each group.
#'
#' @inheritParams match_groups
#' @inheritParams .warn_about_extra_params
#'
#' @return All results found by search method in a list. It raises a
#'   "Convergence failure" error if it cannot find a matched set.
search_heuristic1 <- function(condition,
                              covariates,
                              halting_test,
                              thresh,
                              props,
                              max_removed,
                              print_info = FALSE,
                              ...) {
    .warn_about_extra_params(...)
    sspace <-
        .create_subject_subspace_using_LDA(condition, covariates)
    is.in <- rep(TRUE, length(condition))
    # Computes sample proportions to use if not specified.
    count <- table(condition)
    # Walks the search space.
    depth <- count
    depth[] <- 1
    sspace <- sspace[order(names(sspace))]
    max_removed <- max_removed[names(sspace)]
    limit <- pmin(vapply(sspace, length, 0), max_removed)
    repeat {
        candidates <- names(which(depth < limit))
        if (length(candidates) == 0)
            stop("Convergence failure")
        divergence <- sapply(candidates, function(candidate) {
            count_without_candidate <- count
            dec(count_without_candidate[[candidate]]) <- 1
            .calc_subject_balance_divergence(count_without_candidate, props)
        })
        excess <- names(which.min(divergence))
        is.in[sspace[[excess]][depth[[excess]]]] <- FALSE
        inc(depth[[excess]]) <- 1
        dec(count[[excess]]) <- 1
        if (halting_test(condition[is.in], covariates[is.in, , drop = FALSE], thresh))
            return(list(is.in))
    }
}

#' Finds matching using depth-first search, looking ahead n steps.
#'
#' In each step, it removes one subject from the set of subjects with the
#' smallest associated p-value after "lookahead" steps.
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
#' @keywords internal
.search_heuristic_with_lookahead <-
    function(choose_from_subject_tuples,
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
             ...) {
        .warn_about_extra_params(given_args, ...)
        stopifnot(lookahead >= 1)
        stopifnot(max_removed_percent_per_step <= 1.0)
        stopifnot(ratio_for_slowdown <= 1.0)
        remove_best_only <-
            if (max_removed_per_step <= 1) {
                # either negative or 1
                max_removed_per_step <- abs(max_removed_per_step)
                TRUE
            } else {
                FALSE
            }
        is.in <- rep(TRUE, length(condition))
        look <- 1  # number of subjects removed for checking p-value
        random_choices <- list()
        max_removed_in_next_step <- 1
        step <- 0
        repeat {
            subject_count <- sum(is.in)
            is.in_candidates <-
                is.in &
                (condition %in% names(which(max_removed_per_cond > 0)))
            if (subject_count <= min_preserved ||
                !any(is.in_candidates))
                break
            if (max_removed_per_step > 1) {
                max_removed_in_next_step <- max(1, floor(
                    min(
                        max_removed_per_step,
                        subject_count * max_removed_percent_per_step,
                        subject_count - min_preserved,
                        sum(max_removed_per_cond)
                    )
                ))
                remove_best_only <- (max_removed_in_next_step == 1)
            }
            if (print_info) {
                count <- table(condition[is.in])
                cat(
                    sprintf(
                        "Number of subjects: %s; p/thresh ratio: %f",
                        paste(
                            names(count),
                            count,
                            sep = ':',
                            collapse = ', '
                        ),
                        .calc_p_thresh_ratio(condition[is.in], covariates[is.in, , drop = FALSE], halting_test, thresh)
                    )
                )
                if (max_removed_per_step > 1)
                    cat("; maximum removed in next step:",
                        max_removed_in_next_step)
                cat('\n')
            }
            repeat {
                inc(step) <- 1
                candidates <- iterpc::iterpc(sum(is.in_candidates),
                                             look,
                                             which(is.in_candidates))
                best_sets <- .choose_best_subjects(
                    candidates,
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
                if (print_info)
                    cat("Lookahead ",
                        look,
                        ", number of best sets: ",
                        length(best_sets),
                        "                    \n",
                        sep = '')
                if (length(best_sets) == 0)
                    break
                # if acceptable solution found, returns list of best results
                new_is.in <- .flip_ind(best_sets[[1]], is.in)
                new_ratio <-
                    .calc_p_thresh_ratio(condition[new_is.in], covariates[new_is.in, , drop = FALSE], halting_test, thresh)
                if (new_ratio >= 1.0) {
                    if (print_info) {
                        cat(sprintf(
                            "\nFound %d solution(s) in %d steps\n\n",
                            length(best_sets),
                            step
                        ))
                        if (length(random_choices) > 0)
                            cat(
                                "Random choices:",
                                .vector_list_to_string(random_choices, ": "),
                                "\n"
                            )
                    }
                    return(lapply(best_sets, .flip_ind, is.in))
                }
                if (look >= lookahead)
                    break
                inc(look) <- 1
            }
            if (length(best_sets) == 0)
                break
            # choose best subjects
            best_subject_counts <-
                choose_from_subject_tuples(
                    is.in,
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
                    print_info
                )
            # choose at most max_removed_in_next_step subjects for removal,
            # but not more per condition than given in max_removed_per_cond
            best_subjects <- as.integer(names(best_subject_counts))
            chosen_inds <-
                if (length(best_subjects) <= max_removed_in_next_step) {
                    best_subjects
                } else {
                    # find number of subjects that can be removed from each group: gsc
                    # limit overall and per group, prefering group sizes closest
                    # to what we expect in props
                    can_be_removed_per_cond <-
                        pmin(table(condition[best_subjects]), max_removed_per_cond)
                    to_remove <-
                        min(max_removed_in_next_step,
                            sum(can_be_removed_per_cond))
                    if (sum(can_be_removed_per_cond > 0) == 1) {
                        gsc <- pmin(can_be_removed_per_cond, to_remove)
                    } else if (is.integer(props)) {
                        gsc <- can_be_removed_per_cond  # expected group size changes
                        gsc[] <- 0
                        conds_to_decrease <-
                            sample(setdiff(
                                names(can_be_removed_per_cond),
                                names(props)
                            ))
                        for (cond in c(conds_to_decrease, rev(names(props)))) {
                            gsc[[cond]] <- min(can_be_removed_per_cond[[cond]], to_remove)
                            to_remove <- to_remove - gsc[[cond]]
                        }
                    } else if (is.numeric(props)) {
                        ugs <-
                            table(condition[is.in])  # upper bound
                        lgs <-
                            ugs - can_be_removed_per_cond  # lower bound
                        egs <- pmax(lgs, pmin(ugs, floor((subject_count - to_remove) * props)))  # expected group sizes
                        while (subject_count - sum(egs) > to_remove) {
                            i <- sample(seq_along(egs), 1, prob = ugs - egs)
                            egs[[i]] <- egs[[i]] + 1
                        }
                        repeat {
                            cgs <- pmax(lgs, pmin(ugs, egs))
                            fr <- sum(cgs - egs)  # for removal
                            gtc <- (cgs > lgs)  # groups to change
                            if (fr <= 0 || !any(gtc))
                                break
                            cgs[gtc] <-
                                cgs[gtc] - fr * props[gtc] / sum(props[gtc])  # approximation
                            egs <- cgs
                        }
                        gsc <- ugs - cgs
                    } else {
                        stop("The props parameter has unknown type ", typeof(props))
                    }
                    stopifnot(any(gsc > 0))
                    unlist(lapply(names(gsc)[gsc > 0], function(cond) {
                        candidates <- intersect(best_subjects, which(condition == cond))
                        candidate_count <-
                            min(length(candidates), gsc[[cond]])
                        if (candidate_count > 0)
                            sample(candidates, candidate_count, prob = best_subject_counts[best_subjects %in% candidates])
                        else
                            c()
                    }))
                }
            # if multiple subjects are to be removed
            new_is.in <- .flip_ind(chosen_inds, is.in)
            if (max_removed_per_step > 1) {
                new_ratio <-
                    .calc_p_thresh_ratio(condition[new_is.in], covariates[new_is.in, , drop = FALSE], halting_test, thresh)
                # if we have exceeded ratio for slowdown, do not commit it
                if (new_ratio >= ratio_for_slowdown) {
                    cat(
                        "Slowing down search, since removing subjects would result in p / threshold ratio =",
                        new_ratio,
                        ">=",
                        ratio_for_slowdown,
                        ".\n"
                    )
                    max_removed_per_step <- 1
                    max_removed_in_next_step <- 1
                    remove_best_only <- TRUE
                    next
                } else {
                    cat("Removing",
                        length(chosen_inds),
                        "subjects.\n")
                }
            }
            # committing to best direction
            is.in <- new_is.in
            max_removed_per_cond <- max_removed_per_cond -
                table(condition[chosen_inds])
            if (length(best_subjects) > max_removed_in_next_step)
                random_choices <-
                c(random_choices, list(list(
                    step = step,
                    num_choices = length(best_subjects)
                    # chosen_set = chosen_inds
                )))
        }
        if (print_info)
            cat("Search failed in", step, "steps.\n")
        if (print_info && length(random_choices) > 0)
            cat("Random choices:",
                .vector_list_to_string(random_choices, ": "),
                "\n")
        stop("Convergence failure")
    }

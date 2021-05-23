## Helpers

`inc<-` <- function(x, value)
    x + value
`dec<-` <- function(x, value)
    x - value

#' An infinitesimally small amount, used to check if values are
#' approximately the same.
#'
#' @keywords internal
.tolerance <- sqrt(.Machine$double.eps)


#' Calculates p-value-threshold ratio.
#'
#'
#' @param silent  If FALSE, prints warning when the test statistic cannot be
#'                calculated; if TRUE (the default) they are not printed.
#'
#' @return The p-value-threshold ratio, or NA if the p-value could not be
#'         calculated.
#'
#' @inheritParams match_groups
#'
#' @return The p-value / thresh ratio.
#'
#' @keywords internal
.calc_p_thresh_ratio <- function(condition,
                                 covariates,
                                 halting_test,
                                 thresh,
                                 silent = TRUE) {
    ratio <-
        try(-halting_test(condition, covariates, -thresh), silent = silent)
    if (class(ratio) == "try-error")
        ratio <- NA
    ratio
}


#' Calculates p-value using specified halting test.
#'
#' @inheritParams match_groups
#'
#' @return The p-value.
#'
#' @export
calc_p_value <- function(condition, covariates, halting_test) {
    .calc_p_thresh_ratio(
        condition,
        as.matrix(covariates),
        halting_test,
        thresh = 1,
        silent = FALSE
    )
}


#' Characterizes closeness of actual group sizes to what is expected.
#'
#' @param table_condition  The number of different condition values,
#'                         usually created by calling table(condition).
#'
#' @inheritParams match_groups
#'
#' @return KL divergence of the actual group size proportions from the expected ones.
#'
#' @seealso \code{\link{match_groups}} for meaning of condition parameter.
#'
#' @importFrom entropy KL.plugin
#'
#' @keywords internal
.calc_subject_balance_divergence <-
    function(table_condition, props) {
        if (is.integer(props)) {
            # calculates divergence from preference for keeping all subjects
            divergence <- 0
            multiplier <- 1
            for (cond in names(props)) {
                N <- props[[cond]]
                divergence <-
                    divergence * multiplier + (N - table_condition[[cond]])
                multiplier <- N
            }
            divergence
        } else if (is.numeric(props)) {
            # calculates divergence from group size proportions
            entropy::KL.plugin(prop.table(table_condition), props)
        } else {
            stop("The props parameter has unknown type ", typeof(props))
        }
    }


#' Combines current best and candidate sets, keeping the highest metric value.
#'
#' @param best       A list(metric, sets); metric is a number, set is a list of
#'                   vectors.
#' @param candidate  A list(metric, set); metric is a number, set is a vector.
#'                   Candidate is only considered if metric is not zero.
#'
#' @return A list containing the highest metric and a list of set values (sets).
#'
#' @keywords internal
.combine_sets <- function(best, candidate) {
    if (is.finite(candidate$metric) && candidate$metric &&
        candidate$metric >= best$metric - .tolerance) {
        if (candidate$metric > best$metric + .tolerance) {
            best$metric <- candidate$metric
            best$sets <- list()
        } else {
            best$metric <- min(candidate$metric, best$metric)
        }
        best$sets <- c(best$sets, list(candidate$set))
    }
    best
}


#' Warns about extra (i.e. unused) parameters.
#'
#' @param given_args    The names of arguments given to the search function.
#' @param ...           Consumes extra parameters that are not used by the
#'                      search algorithm at hand; this function gives a warning
#'                      about the ones whose value is not NULL that their value
#'                      is not used.
#'
#'
#' @keywords internal
.warn_about_extra_params <- function(given_args = NULL, ...) {
    ignored_args <- intersect(names(list(...)), given_args)
    if (length(ignored_args) > 0)
        warning("argument(s) ignored: ",
                paste(ignored_args, collapse = ", "),
                "\n")
}


#' Uniquifies a list.
#'
#' @param l  A list.
#'
#' @return The unique list items.
#'
#' @keywords internal
.unique_list <- function(l) {
    j <- 2
    while (j <= length(l)) {
        e <- l[[j]]
        for (i in seq_len(j - 1)) {
            if (isTRUE(all.equal(l[[i]], e, tolerance = .tolerance))) {
                l[[j]] <- NULL
                j <- j - 1
                break
            }
        }
        j <- j + 1
    }
    l
}


#' Flips logical vector at specified indices
#'
#' @param is.in        A logical vector showing which items are preserved.
#'
#' @param ind  Integer indices for the is.in logical vector.
#'
#' @return A logical vector identical to is.in except for indices in ind where
#'         it is is.in negated.
#'
#' @keywords internal
.flip_ind <- function(ind, is.in) {
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
#'
#' @keywords internal
.vector_list_to_string <- function(lv, sep = "") {
    paste(sapply(lv, function(v)
        paste(
            names(v), v, sep = sep, collapse = ", "
        )),
        collapse = "; ")
}


#' Wrapper to foreach::foreach called from .choose_best_subjects.
#'
#' @param input  An iterator created using either the iterpc or the iterators
#' package, or anything else foreach::foreach can interpret (esp. a list).
#'
#' @param operation  The operation to be performed for each item in input
#' (possibly after preprocessing it; see preprocess_input).
#'
#' @param preprocess_input Processes each value retrieved from the input iterator.
#'
#' @param .init,.combine  The same as the parameters of foreach::foreach
#' with identical names.
#'
#' @param max_chunk_size  The maximum number of items to be retrieved from input
#' if it is an iterator.
#'
#' @param print_progress   If TRUE, prints messages about the progress.
#'
#' Used to use iterpc::iter_wrapper() on iterpc iterators, but realized that
#' foreach doesn't handle iterators in a nice way (converts it to a list, which
#' may be huge, instead of gradually retrieving the contensts), so feeding
#' segments of the iterators to foreach instead.
#'
#' @import foreach
#'
#' @keywords internal
.foreach <- function(input,
                     operation,
                     preprocess_input,
                     .init,
                     .combine,
                     max_chunk_size = get("PROCESSED_CHUNK_SIZE", .ldamatch_globals),
                     print_progress = get("PRINT_PROGRESS", .ldamatch_globals)) {
    process <- function(current_input, current_init) {
        foreach::foreach(
            current_input = current_input,
            .init = current_init,
            .combine = .combine,
            .inorder = FALSE
        ) %dopar% operation(preprocess_input(current_input))
    }
    if ('comb' %in% class(input)) {
        get_next <- iterpc::getnext
        input_count <- iterpc::getlength(input)
    } else if ('iter' %in% class(input)) {
        get_next <- function(input) {
            value <- try(iterators::nextElem(input), silent = TRUE)
            if (class(value) == 'try-error')
                NULL
            else
                value
        }
        input_count <- NA
    } else {
        get_next <- NULL
    }
    if (!is.null(get_next)) {
        result <- .init
        input_pos <- 0
        running <- TRUE
        while (running) {
            current_input <-
                replicate(max_chunk_size,
                          get_next(input),
                          simplify = FALSE)
            failure <- sapply(current_input, is.null)
            if (any(failure)) {
                running <- FALSE
                current_input <-
                    current_input[seq_len(which(failure)[[1]] - 1)]
            }
            current_chunk_size <- length(current_input)
            remaining_chunk_count <-
                floor((input_count - (
                    input_pos + current_chunk_size
                )) / max_chunk_size)
            if (print_progress &&
                (input_pos > 0 ||
                 is.na(remaining_chunk_count) ||
                 remaining_chunk_count > 0)) {
                cat(sprintf("\r%d to %d",
                            input_pos + 1,
                            input_pos + current_chunk_size))
                if (!is.na(input_count)) {
                    cat(sprintf(" of %d (%d chunks remaining)",
                                input_count, remaining_chunk_count))
                }
                input_pos <- input_pos + current_chunk_size
            }
            result <- process(current_input, current_init = result)
        }
    } else {
        result <- process(input, .init)
    }
    if (print_progress)
        cat("\r")
    result
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
#' @return list(inds): A list containing the best index vectors indicating the
#' positions to flip in is.in.
#'
#' @keywords internal
.choose_best_subjects <- function(candidates,
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
                                  remove_best_only) {
    inds_as_chars <- as.character(seq_along(condition))

    #' Calculates counts per condition that are within allowed limits.
    #'
    #' @keywords internal
    .relevant_count <- function(counts) {
        sum(pmin(max_removed_per_cond, table(condition[as.integer(names(counts)[counts > 0])])))
    }

    #' Combines candidate lists keeping just enough to exceed maximum values.
    #'
    #' Included within the function body to make parameters of
    #' .choose_best_subjects available as closures instead of having to pass
    #' them as its parameters.
    #'
    #' @param best       list(counts, lsets = list(list(metric, sets)))
    #' @param candidate  list(metric, set)
    #'
    #' @return list(counts, lsets = list(list(metric, sets)))
    #'
    #'
    #' @keywords internal
    .combine_enough_sets <- function(best, candidate) {
        candidate_stored <- FALSE
        lsets_index <- 1
        while (lsets_index <= length(best$lsets)) {
            lset <- best$lsets[[lsets_index]]
            if (candidate$metric >= lset$metric - .tolerance) {
                if (candidate$metric > lset$metric + .tolerance) {
                    best$lsets <-
                        append(best$lsets, list(
                            list(
                                metric = candidate$metric,
                                sets = list(candidate$set)
                            )
                        ), after = lsets_index - 1)
                } else {
                    best$lsets[[lsets_index]] <-
                        list(
                            metric = min(candidate$metric, lset$metric),
                            sets = c(lset$sets, list(candidate$set))
                        )
                }
                counts_for_candidate <-
                    table(factor(unlist(candidate$set), inds_as_chars))
                best$counts <- best$counts + counts_for_candidate
                candidate_stored <- TRUE
                break
            }
            lsets_index <- lsets_index + 1
        }
        if (!candidate_stored) {
            counts_for_candidate <-
                table(factor(unlist(candidate$set), inds_as_chars))
            if (min(max_removed_in_next_step,
                    .relevant_count(best$counts)) < min(
                        max_removed_in_next_step,
                        .relevant_count(best$counts + counts_for_candidate)
                    )) {
                best$lsets <-
                    c(best$lsets, list(list(
                        metric = candidate$metric,
                        sets = list(candidate$set)
                    )))
                best$counts <- best$counts + counts_for_candidate
            }
        } else if (length(best$lsets) > 1) {
            counts_for_last_one <-
                table(factor(unlist(best$lsets[[length(best$lsets)]]), inds_as_chars))
            if (min(max_removed_in_next_step,
                    .relevant_count(best$counts)) == min(
                        max_removed_in_next_step,
                        .relevant_count(best$counts - counts_for_last_one)
                    )) {
                best$lsets <- best$lsets[-length(best$lsets)]
                best$counts <- best$counts - counts_for_last_one
            }
        }
        best
    }

    # filter down subject candidates in at most three steps
    preprocess_input <- as.vector
    if (remove_best_only) {
        combine <- .combine_sets
        init <- list(metric = -Inf, sets = list())
    } else {
        combine <- .combine_enough_sets
        init <-
            list(counts = table(factor(c(), inds_as_chars)), lsets = list())
    }
    for (i in 1:3) {
        if ((i == 1 && prefer_test) || (i == 2 && !prefer_test)) {
            best <- .foreach(
                input = candidates,
                preprocess_input = preprocess_input,
                .init = init,
                .combine = combine,
                operation = function(ind) {
                    cis.in <- .flip_ind(ind, is.in)
                    ratio <-
                        .calc_p_thresh_ratio(condition[cis.in],
                                             covariates[cis.in, , drop = FALSE],
                                             halting_test,
                                             thresh)
                    list(metric = min(ratio, 1.0),
                         set = ind)
                }
            )
            if (!remove_best_only &&
                best$lsets[[1]]$metric >= ratio_for_slowdown) {
                max_removed_in_next_step <- 1
                best$lsets <- best$lsets[1]
            }
        } else if ((i == 2 &&
                    prefer_test) || (i == 1 && !prefer_test)) {
            best <- .foreach(
                input = candidates,
                preprocess_input = preprocess_input,
                .init = init,
                .combine = combine,
                operation = function(ind) {
                    cis.in <- .flip_ind(ind, is.in)
                    divergence <-
                        .calc_subject_balance_divergence(table(condition[cis.in]), props)
                    list(metric = -1 - divergence,
                         set = ind)
                }
            )
        } else if (i == 3 && is.function(tiebreaker)) {
            best <- .foreach(
                input = candidates,
                preprocess_input = preprocess_input,
                .init = init,
                .combine = combine,
                operation = function(ind) {
                    cis.in <- .flip_ind(ind, is.in)
                    ratio <-
                        .calc_p_thresh_ratio(condition[cis.in],
                                             covariates[cis.in, , drop = FALSE],
                                             tiebreaker,
                                             thresh)
                    list(metric = min(ratio, 1.0),
                         set = ind)
                }
            )
        } else {
            break
        }
        preprocess_input <- identity
        if (remove_best_only) {
            candidates <- best$sets
        } else {
            candidates <-
                unlist(lapply(best$lsets, `[[`, 'sets'), recursive = FALSE)
            if (i == 1 && !prefer_test)
                next
        }
        if (.relevant_count(table(unlist(candidates))) <= max_removed_in_next_step)
            break
    }
    candidates
}

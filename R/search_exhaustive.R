#' Creates all group sizes by reducing one group in all rows of grpsizes.
#'
#' Used for generating all group size combinations for one specific total size
#' iteratively, starting from grpsizes with one row containing original group
#' sizes.
#'
#' @param grpsizes  A data.table with the columns containing the group names,
#'                  and the rows containing a particular setup of group sizes.
#'                  All rows are expected to have the same sum (not checked).
#'
#' @param grpnames  The group names (specified because the table can have other
#'                  columns as well).
#'
#' @param minpergrp The minimum number of subjects to be preserved per group.
#'
#' @return A data.table with the same format as grpsizes, containing all
#' possible group setups totaling to one less than the total in grpsizes.
#'
#' @import data.table
#'
#' @keywords internal
.decrease_group_sizes <- function(grpsizes, grpnames, minpergrp) {
    # Table for all distribution of group sizes for one less than previous total.
    d <- data.table::rbindlist(lapply(grpnames, function(col) {
        dg <- data.table::copy(grpsizes)
        dg[, (col) := get(col) - 1]
        dg <- dg[eval(as.name(col)) >= minpergrp[[col]]]
    }))
    data.table::setkeyv(d, grpnames)
    unique(d)
}


#' Orders rows by similarity to expected group size proportions.
#'
#' @inheritParams .decrease_group_sizes
#' @inheritParams match_groups
#'
#' @import data.table
#'
#' @keywords internal
.sort_group_sizes <- function(grpsizes, grpnames, props) {
    divergence <- NULL  # to suppress codetools warnings
    grpsizes[, divergence := vapply(seq_len(nrow(.SD)), function(row)
        .calc_subject_balance_divergence(.SD[row,], props), 0.0),
        .SDcols = grpnames]
    data.table::setorder(grpsizes, divergence)
    grpsizes
}


#' Creates Cartesian product of iterators.
#'
#' @param initializers  A list of initializer functions (with no arguments)
#'                      for iterators.
#'
#' @param get_next      A function for retrieving next item for an iterator
#'                      argument; it assumes that the iterator returns NULL
#'                      when finished.
#'
#' @param sspace        elements to be used (a list of vectors)
#'
#' @return A function that returns list of values, and stops with
#'        "StopIteration" message when finished, so that it can be used with
#'        the iterators::iter() function to create an iterator that works with
#'        foreach.
#'
#' @keywords internal
.create_Cartesian_iterable <-
    function(initializers, get_next, sspace) {
        values <- NULL
        iterators <-
            sapply(initializers, function(fn)
                fn(), simplify = FALSE)
        function() {
            if (is.null(values)) {
                values <<- lapply(seq_along(iterators), function(pos)
                    sspace[[pos]][get_next(iterators[[pos]])])
            } else {
                pos <- length(iterators)
                repeat {
                    index <- get_next(iterators[[pos]])
                    if (!is.null(index)) {
                        values[[pos]] <<- sspace[[pos]][index]
                        inc(pos) <- 1
                        if (pos > length(iterators))
                            break
                    } else {
                        iterators[[pos]] <<- initializers[[pos]]()
                        dec(pos) <- 1
                        if (!pos) {
                            values <<- NULL
                            stop("StopIteration")
                        }
                    }
                }
            }
            values
        }
    }


#' Searches over all possible subspaces for specified group size setup.
#'
#' Results are optimized for the following, in decreasing order of preference:
#' number of subjects; proportion of group sizes close to props;
#' p-value as large as possible.
#'
#' @param best       The best matched groups so far together with its
#'                   p-value / thresh ratio; a list containing ratio and sets
#'                   (a list of subject index vectors).
#'
#' @param grpsize_setup  A set of group sizes as a data.table row (also a list).
#'
#' @param sspace  An ordered subject subspace: a list of vectors,
#' with one vector per group containing the corresponding subject indices.
#'
#' @inheritParams match_groups
#'
#' @inheritParams search_exhaustive
#'
#' @return A list of logical vectors for the best matched groups.
#'
#' @import data.table
#' @import foreach
#' @import gmp
#' @importFrom iterators iter
#' @importFrom iterpc iterpc getnext
#'
#' @keywords internal
.check_subspaces_for_group_size_setup <- function(best,
                                                  grpsize_setup,
                                                  sspace,
                                                  condition,
                                                  covariates,
                                                  halting_test,
                                                  thresh,
                                                  print_info) {
    nz <-
        sapply(names(sspace), function(cond)
            grpsize_setup[[cond]] != 0)
    sspace <- sspace[nz]
    grpsize_setup <- grpsize_setup[, names(nz), with = FALSE]
    ci <-
        .create_Cartesian_iterable(
            sapply(names(sspace), function(cond)
                function()
                    iterpc::iterpc(length(sspace[[cond]]), grpsize_setup[[cond]]), simplify = FALSE),
            iterpc::getnext,
            sspace
        )
    if (print_info) {
        Cartesian_size <- do.call(prod, lapply(
            get("iterators", environment(ci)),
            iterpc::getlength,
            bigz = TRUE
        ))
        cat("Size of Cartesian product:",
            as.character(Cartesian_size),
            "\n")
        search_start_time <- proc.time()
    }
    # calculate halting test values for batch
    best <- .foreach(
        input = iterators::iter(ci),
        preprocess_input = unlist,
        .init = best,
        .combine = .combine_sets,
        operation = function(ind) {
            ratio <-
                halting_test(condition[ind], covariates[ind, , drop = FALSE], thresh)
            list(metric = ratio, set = ind)
        }
    )
    if (print_info) {
        search_time <- (proc.time() - search_start_time)
        cases_per_cpu_second = as.double(
            Cartesian_size /
                (search_time[["user.self"]] + search_time[["sys.self"]]))
        cases_per_wall_clock_second = as.double(
            Cartesian_size / (search_time[["elapsed"]]))
        cat("Number of cases processed per second: ",
            cases_per_cpu_second, " (cpu time) or ",
            cases_per_wall_clock_second, " (wall clock time).\n", sep = '')
    }
    best
}


#' Searches the space backwards, prefering more subjects and certain group size
#' proportions.
#'
#' @details
#' While the search is done in parallel, the search space is enormous and so
#' it can be very slow in the worst case. It is perhaps most useful as a tool
#' to study other matching procedures.
#'
#' You can calculate the maximum possible number of cases to evaluate by
#' calling estimate_exhaustive().
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
#' @import foreach
search_exhaustive <- function(condition,
                              covariates,
                              halting_test,
                              thresh,
                              props,
                              max_removed_per_cond,
                              tiebreaker = NULL,
                              min_preserved = length(levels(condition)),
                              print_info = TRUE,
                              given_args = NULL,
                              ...) {
    .warn_about_extra_params(given_args, ...)
    # Finds best p-value / threshold ratio and the corresponding subsets:
    # iterates over all group sizes, by decreasing groups with original size.
    sspace <- split(seq_along(condition), condition)
    max_removed_per_cond <- max_removed_per_cond[names(sspace)]
    best <- list(metric = -Inf, sets = list())
    grpsizes <- data.table::data.table(t(vapply(sspace, length, 0)))
    minpergrp <- vapply(sspace, length, 0) - max_removed_per_cond
    grpnames <- names(sspace)
    while (best$metric <= 0) {
        grpsizes <-
            .decrease_group_sizes(grpsizes, grpnames, minpergrp)
        if (nrow(grpsizes) == 0)
            break
        total_size <- sum(grpsizes[1, grpnames, with = FALSE])
        if (total_size < min_preserved)
            stop("No subspace found with at least specified minimum total number preserved")
        if (print_info)
            cat(
                "Created",
                nrow(grpsizes),
                "group size configurations",
                "each with a total size of",
                total_size,
                "\n"
            )
        grpsizes <- .sort_group_sizes(grpsizes, grpnames, props)
        # Goes over rows of table with group sizes and finds best suitable
        # subsets, comparing those for group sizes with the same KL divergence
        # only.
        grpsizes_row <- 1
        repeat {
            divergence <- grpsizes[grpsizes_row, divergence]
            repeat {
                if (print_info)
                    cat(paste(names(grpsizes), grpsizes[grpsizes_row], sep = ": "), "\n")
                best <- .check_subspaces_for_group_size_setup(
                    best,
                    grpsizes[grpsizes_row,],
                    sspace,
                    condition,
                    covariates,
                    halting_test,
                    thresh,
                    print_info
                )
                inc(grpsizes_row) <- 1
                if (grpsizes_row > nrow(grpsizes) ||
                    grpsizes[grpsizes_row, divergence] != divergence)
                    break
            }
            if (best$metric > 0 || grpsizes_row > nrow(grpsizes))
                break
        }
    }
    if (best$metric <= 0)
        stop("No subspace found for specified threshold")
    if (is.function(tiebreaker)) {
        best_of_best <- Reduce(function(best, ind) {
            candidate <- list(
                metric = .calc_p_thresh_ratio(condition[ind], covariates[ind, , drop = FALSE], tiebreaker, 1),
                set = ind
            )
            .combine_sets(best, candidate)
        }, best$sets, init = list(metric = -Inf, sets = list()))
        if (length(best_of_best) > 0)
            best <- best_of_best
    }
    lapply(best$sets, function(ind)
        seq_along(condition) %in% ind)
}

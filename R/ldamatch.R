#' ldamatch: Selection of Statistically Similar Research Groups.
#'
#' Select statistically similar research groups by backward selection
#' using various robust algorithms,
#' including a heuristic based on linear discriminant analysis,
#' multiple heuristics based on the test statistic,
#' and parallelized exhaustive search.
#' See the help help for function \code{\link{match_groups}}.
#'
#' @docType package
#' @name ldamatch
NULL


## Actual matching procedure.

#' Normalizes the props parameter for match_groups().
#'
#' @inheritParams match_groups
#' @param keep_all_items  If TRUE and props is a character vector, last item is not dropped.
#'
#' @return A named vector: if props contains proportions, it is the same, but
#' ordered to follow the levels of condition; if props contains names of
#' conditions, the total number of subjects for the condition names in props.
#'
#' @importFrom RUnit checkTrue
#' @importFrom utils head
#'
#' @keywords internal
.normalize_props <- function(props, condition, keep_last_item = FALSE) {
    if (is.null(props)) {
        props <- prop.table(table(condition))
    } else if (is.numeric(props)) {
        RUnit::checkTrue(sum(props) == 1.0, "The sum of props must be 1.0")
        RUnit::checkTrue(
            setequal(names(props), levels(condition)),
            paste(
                "The props vector must have the same names as the conditions"
            )
        )
        props <- props[levels(condition)]
    } else if (is.character(props)) {
        RUnit::checkTrue(length(setdiff(props, levels(condition))) == 0,
                         "Only valid condition names can be used in props")
        RUnit::checkTrue(all(table(props) == 1),
                         "Each condition name can only be listed once in props")
        if (length(props) == length(levels(condition)) && !keep_last_item
        )
            props <-
            utils::head(props, -1)  # all items are not necessary
        props <- table(condition)[props]
    } else {
        stop(
            "The props parameter must be vector of group size proportions, ",
            " or names of the conditions that you wish to keep unchanged",
            " in decreasing order of preference"
        )
    }
    props
}


#' Normalizes max_removed_per_cond parameter for match_groups() and estimate_exhaustive().
#'
#' @inheritParams match_groups
#'
#' @importFrom RUnit checkTrue
#'
#' @keywords internal
.normalize_max_removed_per_cond <-
    function(max_removed_per_cond, condition) {
        RUnit::checkTrue(
            length(names(max_removed_per_cond)) == length(max_removed_per_cond),
            "Values in max_removed_per_cond must be named"
        )
        max_removed_per_cond_ <- table(condition) - 2
        max_removed_per_cond_[names(max_removed_per_cond)] <-
            unlist(max_removed_per_cond)
        max_removed_per_cond <-
            pmax(pmin(max_removed_per_cond_, table(condition) - 1), 0)
        RUnit::checkTrue(
            identical(names(max_removed_per_cond), names(table(condition))),
            "Values in max_removed_per_cond must use the same group names as condition"
        )
        max_removed_per_cond
    }


#' Determines which arguments for a function, which is its caller by default.
#'
#' @param fun    A function; default: the caller.
#' @param ncall  The parent frame index; default: 3 (the great-grandparent).
#' @return A named boolean vector that contains whether each argument is missing.
#'
#' @importFrom methods missingArg
.get_if_args_are_missing <-
    function(fun = sys.function(-1), ncall = 3) {
        f_args <- formals(fun)
        args <- setdiff(names(f_args), "...")
        vapply(args, function(arg_name)
            methods::missingArg(
                as.name(arg_name),
                envir = parent.frame(ncall),
                eval = TRUE
            ),
            FUN.VALUE = TRUE)
    }


#' The available methods for matching.
#' @export
matching_methods = c("heuristic2",
                     "heuristic3",
                     "heuristic4",
                     "random",
                     "exhaustive")

#' The available parallelized methods for matching.
#' @export
parallelized_matching_methods = c("heuristic3", "heuristic4", "exhaustive")


#' The available nondeterministic methods for matching.
#' @export
nondeterministic_matching_methods = c("random", "heuristic3", "heuristic4")


#' Creates a matched group via backward selection.
#'
#' @details
#' The exhaustive, heuristic3, and heuristic4 search methods use the foreach
#' package to parallelize computation.
#' To take advantage of this, you must register a cluster.
#' For example, to use all but one of the CPU cores, run:
#'   \code{doParallel::registerDoParallel(cores = max(1, parallel::detectCores() - 1))}
#' To use sequential processing without getting a warning, run:
#'   \code{foreach::registerDoSEQ()}
#'
#' @param condition     A factor vector containing condition labels.
#'
#' @param covariates    A columnwise matrix containing
#'                      covariates to match the conditions on.
#'
#' @param halting_test  A function to apply to `covariates` (in matrix form)
#'                      which is TRUE iff the conditions are matched.
#'                      Signature: halting_test(condition, covariates, thresh).
#'                      The following halting tests are part of this package:
#'                      \code{\link{t_halt}}, \code{\link{U_halt}},
#'                      \code{\link{l_halt}}, \code{\link{ad_halt}},
#'                      \code{\link{ks_halt}}, \code{\link{wilks_halt}},
#'                      \code{\link{f_halt}}.
#'                      You can create the intersection of two or more halting
#'                      tests using \code{\link{create_halting_test}}.
#'
#' @param thresh        The return value of halting_test has to be greater than
#'                      or equal to thresh for the matched groups.
#'
#' @param method        The choice of search method, one of "random",
#                       "heuristic2", "heuristic3", "heuristic4", and
#                       "exhaustive". The running time increases approximately
#                       in the above order.
#'                      You can get more information about each method on the
#'                      help page for "search_<method_name>"
#'                      (e.g. "\code{\link{search_exhaustive}}").
#'
#' @param props         Either the desired proportions (percentage) of the
#'                      sample for each condition as a named vector,
#'                      or the names of the conditions
#'                      for which we prefer to preserve the subjects,
#'                      in decreasing order of preference. If not specified, the
#'                      (full) sample proportions are used.
#'                      This is preferred among configurations with the same
#                       number of total subjects by the "exhaustive" method, and
#'                      taken into account by the other methods to some extent.
#'                      For example, c(A = 0.4, B = 0.4, C = 0.2) means that
#'                      we would like the number of subjects in groups A, B, and
#'                      C to be around 40\%, 40\%, and 20\% of the total number of
#'                      subjects, respectively. Whereas c("A", "B", "C") means
#'                      that if possible, we would like to keep all subjects
#'                      in group A, and prefer keeping subjects in B, even if
#'                      it results in losing more subjects from C.
#'
#' @param replicates    The maximum number of random replications to be
#'                      performed. This is only used for the "random"
#'                      method.
#'
#' @param print_info    If TRUE, prints summary information on the input and the
#'                      results, as well as progress information for the
#'                      exhaustive search and random algorithms. Default: TRUE;
#'                      can be changed using
#'                      \code{\link{set_param}("PRINT_INFO", FALSE)}.
#'
#' @param min_preserved The minimum number of preserved subjects.
#'                      It can be used to ensure that the search will not take
#'                      forever to run, but instead fail when a solution is not
#'                      found when preserving this number of subjects.
#'
#' @param max_removed_per_cond   A named integer vector, containing the maximum number
#'                      of subjects that can be removed from each group.
#'                      Specify 0 for groups if you want to preserve
#'                      all of their subjects. If you do not specify a value
#'                      for a group, it defaults to 2 less than the group size.
#'                      Values outside the valid range of 0..(N-1)
#'                      (where N is the number of subjects in the group)
#'                      are corrected without a warning.
#'
#' @param tiebreaker    NULL, or a function similar to halting_test, used to
#'                      decide between cases for which halting_test yields
#'                      equal values.
#'
#' @param lookahead     The lookahead to use: a positive integer.
#'                      It is used by the heuristic3 and heuristic4 algorithms,
#'                      with a default of 2.
#'                      The running time is O(N ^ lookahead), wheren N is the
#'                      number of subjects.
#'
#' @param all_results   If TRUE, returns all results found by method in a list.
#'                      (A list is returned even if there is only one result.)
#'                      If FALSE (the default), it returns the first result
#'                      (a logical vector).
#'
#' @param prefer_test   If TRUE, prefers higher test statistic more than
#'                      the expected group size proportion; default is TRUE.
#'                      Used by all algorithms except exhaustive, which always
#                       prefers the group size proportion (otherwise it would
#                       need to evaluate all cases for a certain number of
#                       remaining subjects before it could decide).
#'
#' @param max_removed_per_step   The number of equivalent subjects
#'                      that can be removed in each step. (The actual allowed
#'                      number may be less depending on the p-value / theshold ratio.)
#'                      This parameters is used by the heuristic3 and heuristic4
#'                      algorithms, with a default value of 1.
#'
#' @param max_removed_percent_per_step  The percentage of remaining subjects
#'                      that can be removed in each step.
#'                      Used when max_removed_per_step > 1,
#'                      with a default value of 0.5.
#'
#' @param ratio_for_slowdown  The p-value / threshold ratio at which
#'                      it starts removing subjects one by one.
#'                      Used when max_removed_per_step > 1,
#'                      with a default value of 0.5.
#'
#' @return              A logical vector that contains TRUE for the conditions
#'                      that are in the matched groups;
#'                      or if all_results = TRUE, a list of such vectors.
#'
#' @seealso \code{\link{calc_p_value}} for calculating the test statistic for
#' a group setup.
#' @seealso \code{\link{calc_metrics}} for calculating multiple metrics about
#' the goodness of the result.
#' @seealso \code{\link{compare_ldamatch_outputs}} for comparing multiple
#' different results from this function.
#' @seealso \code{\link{search_heuristic2}, \link{search_heuristic3}, \link{search_heuristic4}, \link{search_random}, \link{search_exhaustive}} for
#  more information on the search algorithms.
#' @export
match_groups <-
    function(condition,
             covariates,
             halting_test,
             thresh = .2,
             method = ldamatch::matching_methods,
             props = prop.table(table(condition)),
             replicates = get("RND_DEFAULT_REPLICATES", .ldamatch_globals),
             min_preserved = length(levels(condition)),
             print_info = get("PRINT_INFO", .ldamatch_globals),
             max_removed_per_cond = NULL,
             tiebreaker = NULL,
             lookahead = 2,
             all_results = FALSE,
             prefer_test = TRUE,
             max_removed_per_step = 1,
             max_removed_percent_per_step = 0.5,
             ratio_for_slowdown = 0.5) {
        ## Finds function for method.
        method = match.arg(method)
        search_method <-
            try(get(paste0("search_", method), mode = "function"),
                silent = TRUE)
        if (class(search_method) == "try-error")
            stop("Search method ", method, " is not available")
        if (print_info)
            cat("Search method: ", method, "\n")
        ## Checks other arguments and create set their values if missing.
        args_missing <- .get_if_args_are_missing()
        if (length(covariates) == 0) {
            if (print_info)
                cat("No covariates specified; including all subjects in output.\n")
            is.in <- rep(TRUE, length(condition))
            return(if (all_results)
                list(is.in)
                else
                    is.in)
        }
        condition <- droplevels(condition)
        covariates <- as.matrix(covariates)
        stopifnot(
            is.factor(condition),
            is.numeric(covariates),
            is.function(halting_test),
            is.numeric(thresh),
            (
                is.function(tiebreaker) || is.null(tiebreaker) ||
                    is.na(tiebreaker)
            )
        )
        RUnit::checkTrue(
            length(condition) == nrow(covariates),
            "There must be one set of covariates for each condition"
        )
        # Checks props argument.
        props <- .normalize_props(props, condition)
        # Checks length of container arguments.
        # Checks and normalizes max_removed_per_cond argument.
        max_removed_per_cond <-
            .normalize_max_removed_per_cond(max_removed_per_cond, condition)
        # Checks for a "natural match" before setting up search.
        if (halting_test(condition, covariates, thresh)) {
            if (print_info)
                cat("Groups are already matched.\n")
            is.in <- rep(TRUE, length(condition))
            return(if (all_results)
                list(is.in)
                else
                    is.in)
        }
        ## Search.
        if (print_info) {
            grpsizes <- table(condition)
            cat("Initial group sizes: ",
                paste(
                    names(grpsizes),
                    grpsizes,
                    sep = ": ",
                    collapse = "\t"
                ),
                "\n")
            cat("Starting", method, "search.\n")
            search_start_time <- proc.time()
        }
        lis.in <- search_method(
            condition = condition,
            covariates = covariates,
            halting_test = halting_test,
            thresh = thresh,
            props = props,
            max_removed_per_cond = max_removed_per_cond,
            tiebreaker = tiebreaker,
            replicates = replicates,
            min_preserved = min_preserved,
            lookahead = lookahead,
            print_info = print_info,
            prefer_test = prefer_test,
            max_removed_per_step = max_removed_per_step,
            max_removed_percent_per_step = max_removed_percent_per_step,
            ratio_for_slowdown = ratio_for_slowdown,
            given_args = names(args_missing)[!args_missing]
        )
        if (print_info) {
            search_time <- (proc.time() - search_start_time)
            cat("Finished ",
                method,
                " search in ",
                search_time[["user.self"]] + search_time[["sys.self"]],
                " seconds (wall click time passed:", search_time[["elapsed"]],
                ").\n", sep = '')
            grpsizes <- table(condition[lis.in[[1]]])
            cat(
                "Eventual group sizes:",
                paste(
                    names(grpsizes),
                    grpsizes,
                    sep = ": ",
                    collapse = "\t"
                ),
                "\n"
            )
            grpremoved <- table(condition) - grpsizes
            cat("Removed subjects:    ",
                paste(
                    names(grpremoved),
                    grpremoved,
                    sep = ": ",
                    collapse = "\t"
                ),
                "\n")
            cat(
                "The p-value before matching:",
                calc_p_value(condition, covariates, halting_test),
                "\n"
            )
            cat(
                "The p-values after matching:",
                paste(names(sort(table(sapply(lis.in, function(b) calc_p_value(condition[b], covariates[b, , drop = FALSE], halting_test))), decreasing = TRUE)), collapse = ', '),
                "\n"
            )
        }
        if (all_results)
            lis.in
        else
            lis.in[[1]]
    }

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

#' Creates ordered subspace of subject candidates for removal, using LDA.
#' Excludes constant covariates from consideration.
#'
#' @inheritParams match_groups
#'
#' @return An ordered subject subspace: a list of vectors, with one vector per
#' group containing the corresponding subject indices.
#'
#' @importFrom MASS lda
#' @importFrom stats coef var
.create_subject_subspace_using_LDA <-
    function(condition, covariates) {
        ## Computes linear projection vector, after removing constant columns.
        covariates <- covariates[, apply(covariates, 2, stats::var, na.rm = TRUE) != 0, drop = FALSE]
        if (ncol(covariates) == 0) {
            stop("No non-constant variables in covariates")
        } else if (ncol(covariates) == 1) {
            W <- 1  # just uses the identity projection.
        } else {
            W <- stats::coef(MASS::lda(condition ~ covariates))[, 1]
        }
        projection <- as.vector(covariates %*% W)
        ## Set up search space.
        # Computes means.
        mu <- mean(projection)
        mu.conditions <- tapply(projection, condition, mean)[condition]
        # Locates observations driving condition/covariate correlation(s).
        correlates <-
            ((projection < mu.conditions) == (mu.conditions < mu))
        # Computes order, filtering out non-correlates.
        ord <- order(projection)
        ord <- ord[correlates[ord]]
        # Splits on condition and reverse for those above.
        sspace <- split(ord, condition[ord])
        for (name in names(sspace))
            if (projection[sspace[[name]][1]] > mu)
                sspace[[name]] <- rev(sspace[[name]])
        sspace
    }


#' Normalizes the props parameter for match_groups().
#'
#' @inheritParams match_groups
#'
#' @return A named vector: if props contains proportions, it is the same, but
#' ordered to follow the levels of condition; if props contains names of
#' conditions, the total number of items for the condition names in props.
#'
#' @importFrom RUnit checkTrue
#' @importFrom utils head
.normalize_props <- function(props, condition) {
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
        if (length(props) == length(levels(condition)))
            props <- utils::head(props, -1)
        props <- table(condition)[props]
    } else {
        stop(
            "The props parameter must be vector of group size proportions, ",
            " or names of the conditions that you wish to keep unchanged",
            " (in decreasing order of preference)"
        )
    }
    props
}


#' Normalizes max_removed parameter for match_groups() and estimate_exhaustive().
#'
#' @inheritParams match_groups
#'
#' @importFrom RUnit checkTrue
.normalize_max_removed <- function(max_removed, condition) {
    RUnit::checkTrue(length(names(max_removed)) == length(max_removed),
                     "Values in max_removed must be named")
    max_removed_ <- table(condition) - 1
    max_removed_[names(max_removed)] <- unlist(max_removed)
    max_removed <- pmax(pmin(max_removed_, table(condition) - 1), 0)
    RUnit::checkTrue(
        identical(names(max_removed), names(table(condition))),
        "Values in max_removed must use the same group names as condition"
    )
    max_removed
}


#' Creates a matched group via backward selection.
#'
#' @details
#' The exhaustive, heuristic3, and heuristic4 search methods use the foreach
#' package to parallelize computation.
#' To take advantage of this, you must register a cluster.
#' For example, to use all but one of the CPU cores, run:
#'   \code{doMC::registerDoMC(max(1, parallel::detectCores() - 1))}
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
#' @param method        The choice of search method, one of "heuristic1"
#'                      (formerly called "heuristic"), "random", "heuristic2",
#'                      "heuristic3", "heuristic4", and "exhaustive".
#'                      The running time increases approximately in the above
#'                      order.
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
#'                      This is enforced by the "heuristic1" method,
#'                      preferred among configurations with the same number of
#'                      total subjects by the "exhaustive" method, and
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
#' @param max_removed   A named integer vector, containing the maximum number
#'                      of subjects that can be removed from each group.
#'                      Specify 0 for groups if you want to preserve
#'                      all of their subjects. If you do not specify a value
#'                      for a group, it defaults to one less than the group size.
#'                      Values outside the valid range of 0..(N-1)
#'                      (where N is the number of subjects in the group)
#'                      are corrected without a warning.
#'
#' @param tiebreaker    NULL, or a function similar to halting_test, used to
#'                      decide between cases for which halting_test yields
#'                      equal values.
#'
#' @param lookahead     The lookahead to use: a positive integer.
#'                      it is used by the heuristic3 and heuristic4 algorithms,
#'                      with a default of 2. As you increase it,
#'                      the running time increases exponentially.
#'
#' @param all_results   If TRUE, returns all results found by method in a list.
#'                      (A list is returned even if there is only one result.)
#'                      If FALSE (the default), it returns the first result
#'                      (a logical vector).
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
#'
#' @export
match_groups <-
    function(condition,
             covariates,
             halting_test,
             thresh = .2,
             method = c("heuristic1",
                        "random",
                        "heuristic2",
                        "heuristic3",
                        "heuristic4",
                        "exhaustive"),
             props = prop.table(table(condition)),
             replicates = NULL,
             min_preserved = NULL,
             print_info = get("PRINT_INFO", .ldamatch_globals),
             max_removed = NULL,
             tiebreaker = NULL,
             lookahead = NULL,
             all_results = FALSE) {
        ## Finds function for method.
        method = match.arg(method)
        search_method <-
            try(get(paste0("search_", method), mode = "function"),
                silent = TRUE)
        if (class(search_method) == "try-error")
            stop("Search method ", method, " is not available")
        ## Checks other arguments and create set their values if missing.
        if (length(covariates) == 0) {
            if (print_info)
                cat("No covariates specified, including all in output.\n")
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
        # Checks and normalizes max_removed argument.
        max_removed <- .normalize_max_removed(max_removed, condition)
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
            max_removed = max_removed,
            tiebreaker = tiebreaker,
            replicates = replicates,
            min_preserved = min_preserved,
            lookahead = lookahead,
            print_info = print_info
        )
        if (print_info) {
            total_search_time <- (proc.time() - search_start_time)[["elapsed"]]
            cat("Finished",
                method,
                "search in",
                total_search_time,
                "seconds.\n")
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
        }
        if (all_results)
            lis.in
        else
            lis.in[[1]]
    }

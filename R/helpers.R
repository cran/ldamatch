## Helpers

`inc<-` <- function(x, value)
    x + value
`dec<-` <- function(x, value)
    x - value


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
.calc_p_thresh_ratio <- function(condition,
                                 covariates,
                                 halting_test,
                                 thresh,
                                 silent = TRUE) {
    ratio <-
        try(-halting_test(condition, covariates,-thresh), silent = silent)
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
#' @return KL divergence between actual and expected group size proportions.
#'
#' @seealso \code{\link{match_groups}} for meaning of condition parameter.
#'
#' @importFrom entropy KL.plugin
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
        } else {
            # calculates divergence from group size proportions
            entropy::KL.plugin(prop.table(table_condition), props)
        }
    }


#' Combines current best and candidate sets, keeping the highest metric value.
#'
#' @param candidate  A list containing metric (a number) and set (a vector).
#'                   Candidate is only considered if metric is not zero.
#'
#' @inheritParams .check_subspaces_for_group_size_setup
#'
#' @return A list containing the highest metric and a list of set values (sets).
.combine_sets <- function(best, candidate) {
    if (is.finite(candidate$metric) && candidate$metric &&
        candidate$metric >= best$metric) {
        if (candidate$metric > best$metric) {
            best$metric <- candidate$metric
            best$sets <- list()
        }
        best$sets <- c(best$sets, list(candidate$set))
    }
    best
}


#' Warns about extra (i.e. unused) parameters.
#'
#' @param ...           Consumes extra parameters that are not used by the
#'                      search algorithm at hand; this function gives a warning
#'                      about the ones whose value is not NULL that their value
#'                      is not used.
#'
.warn_about_extra_params <- function(...) {
    extra_params <- list(...)
    for (n in names(extra_params)) {
        if (!is.null(extra_params[[n]]))
            warning(n, " parameter ignored")
    }
}


#' Uniquifies a list.
#'
#' @param l  A list.
#'
#' @return The unique list items.
.unique_list <- function(l) {
    j <- 2
    while (j <= length(l)) {
        e <- l[[j]]
        for (i in seq_len(j - 1)) {
            if (isTRUE(all.equal(l[[i]], e))) {
                l[[j]] <- NULL
                j <- j - 1
                break
            }
        }
        j <- j + 1
    }
    l
}

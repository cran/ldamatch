## Halting tests.


#' A univariate halting test using the t-test, which must be satisfied for all condition pairs.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @export
t_halt <- function(condition, covariates, thresh) {
    .apply_crit_to_condition_pairs(covariates, .t_crit, condition, thresh)
}


#' Criterion function for t_halt.
#'
#' @param covariate     A vector containing a covariate to match the conditions on.
#'
#' @inheritParams match_groups
#'
#' @return The p-value.
#'
#' @importFrom stats t.test
.t_crit <- function(covariate, condition) {
    stats::t.test(covariate ~ condition)$p.value
}


#' A univariate halting test using the Wilcoxon test, which must be satisfied for all condition pairs.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @export
U_halt <- function(condition, covariates, thresh) {
    .apply_crit_to_condition_pairs(covariates, .U_crit, condition, thresh)
}


#' Criterion function for U_halt.
#'
#' @inheritParams .t_crit
#'
#' @return The p-value.
#'
#' @importFrom stats wilcox.test
.U_crit <- function(covariate, condition) {
    stats::wilcox.test(covariate ~ condition)$p.value
}


#' A univariate halting test using Levene's test.
#'
#' Warnings such as "ANOVA F-tests on an essentially perfect fit are unreliable"
#' are suppressed.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @export
l_halt <- function(condition, covariates, thresh) {
    .apply_crit(covariates, .l_crit, condition, thresh)
}


#' Criterion function for l_halt.
#'
#' Warnings such as "ANOVA F-tests on an essentially perfect fit are unreliable"
#' are suppressed.
#'
#' @inheritParams .t_crit
#'
#' @return The p-value.
#'
#' @importFrom car leveneTest
.l_crit <- function(covariate, condition) {
    suppressWarnings(car::leveneTest(covariate ~ condition)["group", "Pr(>F)"])
}


#' A univariate halting test using the Anderson-Darling test.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @export
ad_halt <- function(condition, covariates, thresh) {
    .apply_crit(covariates, .ad_crit, condition, thresh)
}


#' Criterion function for ad_halt.
#'
#' @inheritParams .t_crit
#'
#' @return The p-value.
#'
#' @importFrom kSamples ad.test
.ad_crit <- function(covariate, condition) {
    lkS <- kSamples::ad.test(
        split(covariate, condition),
        method = get("AD_METHOD", .ldamatch_globals),
        Nsim = get("AD_NSIM", .ldamatch_globals)
    )
    lkS$ad[[get("AD_VERSION", .ldamatch_globals), 3]]
}


#' A univariate halting test using the Kolmogorov-Smirnov Test, which must be satisfied for all condition pairs.
#'
#' The condition must have two levels.
#'
#' Note that unlike many tests, the null hypothesis is that the two samples are
#' are drawn from the same distribution.
#'
#' Warnings such as "cannot compute exact p-value with ties" are suppressed.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @importFrom RUnit checkEquals
#'
#' @export
ks_halt <- function(condition, covariates, thresh) {
    .apply_crit_to_condition_pairs(covariates, .ks_crit, condition, thresh)
}

#' Criterion function for ks_halt.
#'
#' Warnings such as "cannot compute exact p-value with ties" are suppressed.
#'
#' @inheritParams .t_crit
#'
#' @return The p-value.
#'
#' @importFrom stats ks.test
.ks_crit <- function(covariate, condition) {
    cc = split(covariate, condition)
    suppressWarnings(stats::ks.test(cc[[1]], cc[[2]])$p.value)
}


#' A multivariate halting test appropriate for more than two condition levels.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @export
wilks_halt <- function(condition, covariates, thresh) {
    p <- min(summary(manova(covariates ~ droplevels(condition)),
                     test = "Wilks")$stats[1, 6])
    if (p < thresh)
        return(0.0)
    p / thresh
}


#' A univariate halting test using Fisher's exact test.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @export
f_halt <- function(condition, covariates, thresh) {
    .apply_crit(covariates, .f_crit, condition, thresh)
}


#' Criterion function for f_halt.
#'
#' @inheritParams .t_crit
#'
#' @return The p-value.
#'
#' @importFrom stats fisher.test
.f_crit <- function(covariate, condition) {
    stats::fisher.test(covariate, condition)$p.value
}


#' Returns halting tests for names, or checks if pass functions are suitable.
#'
#' @param halting_test  The name of one halting test, or a halting test function.
#'
#' @return A vector of halting test functions.
#'
#' @importFrom RUnit checkTrue
.get_halting_test <- function(halting_test) {
    halting_test <- if (is.character(halting_test))
        get(halting_test, mode = "function")
    else
        halting_test
    RUnit::checkTrue(
        length(formals(halting_test)) == 3,
        "A halting_test must be a function with three parameters (condition, covariates, thresh)"
    )
    halting_test
}


#' Recycles threshold values for halting tests.
#'
#' @param ts  Threshold value(s).
#'
#' @param hs  Halting tests.
#'
#' @return A vector with one threshold value per halting test.
.recycle <- function(ts, hs) {
    reps <- length(hs) / length(ts)
    if (reps != round(reps))
        stop("number of thresholds not a multiple of number of halting tests")
    rep(ts, reps)
}


#' Creates halting test from multiple tests.
#'
#' The created halting test function returns the smallest p-value-to-threshold
#' ratio of the values produced by the supplied tests, or zero if any of the
#' p-values does not exceed the threshold. The resulting function expects one
#' threshold per halting test in a vector or it recycles the given value(s) to
#' get a threshold for each one.
#'
#' @param halting_tests   Either a vector of halting test functions
#'                        (or function names) with the signature
#'                        halting_test(condition, covariates, thresh)
#'                        (for the meaning of the parameters see
#'                        \code{\link{match_groups}}); or it may be a list of
#'                        list(test = halting_test, cond = subset_of_conditions,
#'                             cov = variable_selector, thresh) fields.
#'                        All fields can be left out except test, and test need
#'                        not be named if it is the first item in the list.
#'                        The subset_of_conditions can be names of the
#'                        conditions to match (a character vector or a factor).
#'                        The variable_selector can be a logical vector with as
#'                        many items as there will be columns in covariates
#'                        (recommended), or a vector of integer covariate
#'                        column indices.
#'                        Each halting_test is then only applied to the
#'                        specified subset of conditions and variables of the
#'                        covariate matrix, with the specified threshold; when
#'                        a value is not specified the defaults are used.
#'                        Note that ordering the functions does not change the
#'                        behavior, but can make the execution of the combined
#'                        function faster, as the later ones are often evaluated
#'                        only if the criteria for the earlier ones is met.
#'
#' @return A function that returns the minimum of all halting test values;
#'         the threshold value supplied to it is recycled for the individual
#'         functions.
#'
#' @importFrom RUnit checkTrue
#'
#' @export
create_halting_test <- function(halting_tests) {
    # check input parameters and convert halting test names into functions
    if (all(sapply(halting_tests, is.list))) {
        # hs: list(list(halting_test, column_index), ...)
        hs <- lapply(halting_tests, function(h) {
            # make sure only the allowed field names are used
            RUnit::checkTrue(
                all(grepl(
                    "^(|test|cond|cov|thresh).*$", names(h)
                )),
                paste0(
                    "Only the parameters with the following names can be",
                    " given for halting tests:",
                    " 'test' (can be unnamed if first in list),",
                    "'cond', 'cov', and 'thresh'"
                )
            )
            # make sure test is the first one and named
            h <- if (is.null(h$test)) {
                c(list(test = h[[1]]), h[-1])
            } else {
                c(list(test = h$test), h[!grepl("^test", names(h))])
            }
            # check fields
            RUnit::checkTrue(
                class(h$test) %in% c("function", "character"),
                paste0(
                    "The halting test must be specified as",
                    " a function or a function name"
                )
            )
            RUnit::checkTrue(
                class(h$cond) %in% c("NULL", "character", "factor"),
                paste0(
                    "The conditions for the halting test (specified in",
                    " 'cond') must be either character strings or a factor ",
                    "variable"
                )
            )
            RUnit::checkTrue(
                class(h$cov) %in% c("NULL", "logical", "integer", "numeric"),
                paste0(
                    "The variable selector for the halting test",
                    "(specified in 'cov') must be either a logical vector,",
                    " or a vector of integers"
                )
            )
            RUnit::checkTrue(
                class(h$thresh) %in% c("NULL", "numeric"),
                paste0(
                    "The threshold for the halting test must be",
                    " either a number or NULL (not specified)"
                )
            )
            # make sure test is a halting test (possibly convert test name to
            # function)
            list(
                test = .get_halting_test(h$test),
                cov = h$cov,
                cond = h$cond,
                thresh = h$thresh
            )
        })
    } else {
        hs <- sapply(halting_tests, .get_halting_test, simplify = FALSE)
        if (length(hs) == 1)
            return(hs[[1]])
        hs <-
            lapply(hs, function(h)
                list(test = .get_halting_test(h)))
    }
    function(condition, covariates, threshes) {
        sign_threshes <- unique(sign(threshes))
        stopifnot(length(sign_threshes) == 1)
        if (sign_threshes > 0.0) {
            choose_worse <- min
            worst_ratio <- Inf
        } else {
            choose_worse <- max
            worst_ratio <- -Inf
        }
        ts <- .recycle(threshes, hs)
        for (i in seq_along(hs)) {
            h <- hs[[i]]
            h$covariates <- if (is.null(h$cov)) {
                covariates
            } else {
                covariates[, h$cov, drop = FALSE]
            }
            if (is.null(h$cond)) {
                h$condition <- condition
            } else {
                b <- (condition %in% h$cond)
                h$covariates <- h$covariates[b, , drop = FALSE]
                h$condition <- droplevels(condition[b])
            }
            if (is.null(h$thresh))
                h$thresh <- ts[[i]]
            ratio <- h$test(h$condition, h$covariates, h$thresh)
            if (!ratio)
                return(0.0)
            worst_ratio <- choose_worse(worst_ratio, ratio)
        }
        worst_ratio
    }
}


#' Returns smallest halting_test-threshold ratio, or 0 if less than 1.
#'
#' @param crit        The criterion function to use, such as \code{\link{t_crit}}.
#'
#' @inheritParams match_groups
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
#' @importFrom stats manova
.apply_crit <- function(covariates, crit, condition, thresh) {
    min_p <- Inf
    condition <- droplevels(condition)
    for (i in seq_len(ncol(covariates))) {
        p <- try(crit(covariates[, i], condition), silent = TRUE)
        if (class(p) == "try-error") {
            n_levels = length(levels(condition))
            if (n_levels == 1 ||
                grepl("not enough.*observations", p) ||
                #temp for a kSamples::ad.test() error
                grepl(".Fortran(C_bvalus, n = as.integer(n),", p)) {
                return(0.0)  # problem is too few observations
            } else if (grepl("'x' and 'y' must have at least 2 levels", p) ||
                       grepl("data.*constant", p)) {
                p <- 1.0  # data is constant, so matched
            } else {
                warning(p)
                return(0.0)
            }
        } else if (!is.finite(p)) {
            # e.g. l_crit for 1 subject per group)
            return(NA)
        } else if (p < thresh) {
            return(0.0)
        }
        min_p <- min(p, min_p)
    }
    min_p / thresh
}


#' Returns smallest value from .apply_crit for all condition pairs.
#'
#' @inheritParams .apply_crit
#'
#' @return The ratio of the p-value and the threshold, or 0 if the p-value is
#'         less than the threshold.
#'
.apply_crit_to_condition_pairs <-
    function(covariates, crit, condition, thresh) {
        condition <- droplevels(condition)
        if (length(levels(condition)) <= 2)
            return(.apply_crit(covariates, crit, condition, thresh))
        if (sign(thresh) > 0.0) {
            choose_worse <- min
            worst_ratio <- Inf
        } else {
            choose_worse <- max
            worst_ratio <- -Inf
        }
        for (condition_pair in combn(levels(condition), 2, simplify = FALSE)) {
            b <- (condition %in% condition_pair)
            ratio <- .apply_crit(covariates[b, , drop = FALSE], crit,
                                 condition[b], thresh)
            if (!ratio)
                return(0.0)
            worst_ratio <- choose_worse(worst_ratio, ratio)
        }
        worst_ratio
    }

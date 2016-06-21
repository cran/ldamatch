#' Compares outputs of ldamatch runs.
#'
#' It favors, in decreasing order of priority, fewer excluded subjects,
#' better balance (i.e. subsamples that diverge less from the expected
#' proportions, which are by default the proportions of the input groups), and
#' better (i.e. larger) test statistic for the matched groups.
#'
#' @param is.in1        A logical vector for output 1, TRUE iff row is in the match.
#'
#' @param is.in2        A logical vector for output 2, TRUE iff row is in the match.
#'
#' @param prefer_test   If TRUE, prefers higher test statistic more than
#'                      the group size proportion; default is FALSE if props
#'                      is specified, TRUE if it is not.
#'
#' @inheritParams match_groups
#'
#' @return A number that is > 0 if is.in1 is a better solution than is.in2,
#'         < 0 if is.in1 is a worse solution than is.in2, or
#'         0 if the two solutions are equivalent (not necessarily identical).
#'
#' @importFrom RUnit checkTrue
#'
#' @export
compare_ldamatch_outputs <- function(is.in1,
                                     is.in2,
                                     condition,
                                     covariates = matrix(),
                                     halting_test = NA,
                                     props = NULL,
                                     prefer_test = is.null(props),
                                     tiebreaker = NULL) {
    exception1 <- (class(is.in1) == "try-error")
    exception2 <- (class(is.in2) == "try-error")
    if (exception1 || exception2)
        return(exception2 - exception1)
    covariates <- as.matrix(covariates)
    RUnit::checkTrue(
        length(is.in1) == length(is.in2),
        "The is.in1 and is.in2 parameters must have the same length"
    )
    RUnit::checkTrue(
        length(is.in1) == length(condition),
        "The is.in1 and condition parameters must have the same length"
    )
    RUnit::checkTrue(
        !is.function(halting_test) || length(condition) == nrow(covariates),
        paste0(
            "The length of the condition parameter must be the same",
            " as the number of covariates"
        )
    )
    # compare number of excluded subjects
    num_excluded1 <- sum(!is.in1)
    num_excluded2 <- sum(!is.in2)
    if (num_excluded1 != num_excluded2)
        return(num_excluded2 - num_excluded1)
    for (test_stat in c(prefer_test,!prefer_test)) {
        # compare test statistic
        if (test_stat && is.function(halting_test)) {
            p_matched1 <- calc_p_value(condition[is.in1], covariates[is.in1, , drop = FALSE], halting_test)
            p_matched2 <- calc_p_value(condition[is.in2], covariates[is.in2, , drop = FALSE], halting_test)
            if (p_matched1 != p_matched2)
                return(p_matched1 - p_matched2)
            # compare tiebreaker statistic
            if (is.function(tiebreaker)) {
                p_tiebreaker1 <- calc_p_value(condition[is.in1], covariates[is.in1, , drop = FALSE], tiebreaker)
                p_tiebreaker2 <- calc_p_value(condition[is.in2], covariates[is.in2, , drop = FALSE], tiebreaker)
                if (p_tiebreaker1 != p_tiebreaker2)
                    return(p_tiebreaker1 - p_tiebreaker2)
            }
        } else {
            # compare divergence from expected balance
            props <- .normalize_props(props, condition)
            divergence1 <- .calc_subject_balance_divergence(table(condition[is.in1]), props)
            divergence2 <- .calc_subject_balance_divergence(table(condition[is.in2]), props)
            if (divergence1 != divergence2)
                return(divergence2 - divergence1)
        }
    }
    return(0)
}


#' Calculates basic metrics about ldamatch search result.
#'
#' @param is.in        The output of match_groups(): either a logical vector,
#'                     or a list of those.
#'
#' @return A list containing: \describe{
#' \item{all.is.in}{all results as a list;}
#' \item{is.in}{simply the first item in all.is.in or the error contained in is.in;}
#' \item{num_excluded}{the number of excluded subjects), p_matched
#'    (the test statistic from halting_test for the matched groups);}
#' \item{p_tiebreaker}{the test statistic from tiebreaker for the matched groups; and}
#' \item{balance_divergence}{a value characterizing the deviation
#'   from the expected group size proportions specified in props.}
#' }
#' If the value for a field cannot be calculated, it will still be present
#' with a value of NA.
#'
#' @inheritParams match_groups
#'
#' @export
calc_metrics <- function(is.in,
                         condition,
                         covariates,
                         halting_test,
                         props = prop.table(table(condition)),
                         tiebreaker = NULL) {
    if (class(is.in) == "try-error") {
        return(
            list(
                all.is.in = list(),
                is.in = is.in,
                num_excluded = NA,
                p_matched = NA,
                p_tiebreaker = NA,
                balance_divergence = NA
            )
        )
    }
    all.is.in <- if (is.list(is.in))
        is.in
    else
        list(is.in)
    is.in <- all.is.in[[1]]
    if (length(covariates) == 0) {
        p_matched <- Inf
        p_tiebreaker <- Inf
    } else {
        covariates <- as.matrix(covariates)
        p_matched <- suppressWarnings(calc_p_value(condition[is.in], covariates[is.in, , drop = FALSE], halting_test))
        p_tiebreaker <- if (is.function(tiebreaker)) {
            suppressWarnings(calc_p_value(condition[is.in], covariates[is.in, , drop = FALSE], tiebreaker))
        } else {
            NA
        }
    }
    list(
        all.is.in = all.is.in,
        is.in = is.in,
        num_excluded = sum(!is.in),
        p_matched = p_matched,
        p_tiebreaker = p_tiebreaker,
        balance_divergence = .calc_subject_balance_divergence(table(condition[is.in]), .normalize_props(props, condition))
    )
}

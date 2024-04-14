context("ldamatch")


# SETTINGS

# set default parameter values
set_param("PRINT_INFO", FALSE)
set_param("RND_DEFAULT_REPLICATES", 10)

# halting tests to be used for numeric covariates
halting_tests_for_numbers <- list(
    U_halt = U_halt, l_halt = l_halt, ks_halt = ks_halt, wilks_halt = wilks_halt,
    t_ad_halt = create_halting_test(list(t_halt, ad_halt)))

# halting tests to be used for factor covariates
halting_tests_for_factors <- list(f_halt = f_halt)  # halting tests


# CREATE TEST DATA

# univariate data
set.seed(257)
SIZE <- 2
condition <- as.factor(c(rep("control", 2 * SIZE), rep("treatment", SIZE)))
covariate <- c(rnorm(2 * SIZE), rnorm(SIZE, 1, 2))
covariate_factor <- as.numeric(as.factor(
    c(rbinom(2 * SIZE, 1, 0.5), rbinom(SIZE, 1, 0.8))))

# multivariate data

set.seed(257)
SIZE <- 2
condition <- as.factor(c(rep("control", 2 * SIZE), rep("treatment", SIZE)))

covariate1 <- c(rnorm(2 * SIZE), rnorm(SIZE, 1, 2))
covariate2 <- c(rnorm(2 * SIZE), rnorm(SIZE, 1, 2))
covariates <- cbind(covariate1, covariate2)

covariate1_factor <- as.factor(c(rbinom(2 * SIZE, 1, 0.5), rbinom(SIZE, 1, 0.8)))
covariate2_factor <- as.factor(c(rbinom(2 * SIZE, 1, 0.5), rbinom(SIZE, 1, 0.8)))
covariates_factor <- cbind(as.numeric(covariate1_factor),
                           as.numeric(covariate2_factor))


# FUNCTIONS

test_ldamatch_with_tiebreaker <- function(
        condition, covariates, halting_test = t_halt, method = "heuristic2") {
    expect_silent(match_groups(condition, covariates, halting_test,
                               method = method, tiebreaker = t_halt))
}


test_ldamatch_with_proportions <- function(
        condition, covariates, halting_test = t_halt, method = "heuristic2") {
    expect_error(match_groups(condition, covariates, halting_test,
                              method = method, props = c(a = 0.5, b = 0.5)),
                 info = "props badly named")
    expect_error(match_groups(condition, covariates, halting_test,
                              method = method,
                              props = c(treatment = 0.5, control = 0.7)),
                 info = "props does not sum to 1")
    expect_silent(match_groups(condition, covariates, halting_test,
                               method = method,
                               props = c(treatment = 0.9, control = 0.1)))
}


test_ldamatch_with_DX_preferences <- function(
        condition, covariates, halting_test = t_halt, method = "heuristic2") {
    expect_error(match_groups(condition, covariates, halting_test,
                              method = method, props = c("a")),
                 info = "DX preference bad")
    expect_silent(match_groups(condition, covariates, halting_test,
                               method = method, props = c("treatment")))
}


test_ldamatch_methods <- function(
        condition, covariates, halting_test = t_halt) {
    # test all methods
    for (method in matching_methods) {
        # we must get output when print_info is TRUE; return first result
        expect_output(
            (is.in = match_groups(condition, covariates, halting_test, method = method,
                                  print_info = TRUE, all_results = FALSE)),
                      info = paste(method, "print info, not all results"))
        expect_is(is.in, "logical")
        expect_length(is.in, length(condition))
        # we must not get output when print_info is FALSE; return all results
        expect_silent(
            (is.in = match_groups(condition, covariates, halting_test, method = method,
                                  print_info = FALSE, all_results = TRUE)))
        expect_is(is.in, "list")
        expect_length(is.in[[1]], length(condition))
    }
}


test_ldamatch_halting_tests <- function(
    condition, covariates, halting_tests = halting_tests_for_numbers,
    method = "heuristic2") {
    # test with all halting tests
    expect_false(is.null(names(halting_tests)))
    for (halting_test_name in names(halting_tests)) {
        halting_test = halting_tests[[halting_test_name]]
        # test with different thresh values
        expect_true(all(match_groups(condition, covariates, halting_test,
                                     method = method, thresh = 0)),
                    info = paste(method, halting_test_name, "thresh = 0"))
        expect_error(match_groups(condition, covariates, halting_test,
                                  method = method, thresh = 1 + 1e-10),
                     info = paste(method, halting_test_name, "thresh = 1"))
    }
}


# TESTS

foreach::registerDoSEQ()

# univariate data

test_that("ldamatch works on univariate numeric data with tiebreaker", {
    test_ldamatch_with_tiebreaker(condition, covariate)
})

test_that("ldamatch works on univariate numeric data with proportions", {
    test_ldamatch_with_proportions(condition, covariate)
})

test_that("ldamatch works on univariate numeric data with DX preferences", {
    test_ldamatch_with_DX_preferences(condition, covariate)
})

test_that("ldamatch works on univariate numeric data with all methods", {
    test_ldamatch_methods(condition, covariate)
})

test_that("ldamatch works on univariate numeric data with all halting tests", {
    test_ldamatch_halting_tests(
        condition, covariate, halting_tests_for_numbers[
            names(halting_tests_for_numbers) != "wilks_halt"])
})

test_that("ldamatch works on univariate factor data with all methods", {
    test_ldamatch_methods(condition, covariate_factor)
})

test_that("ldamatch works on univariate factor data with all halting tests", {
    test_ldamatch_halting_tests(condition, covariate_factor,
                                halting_tests_for_factors)
})


# multivariate data

test_that("ldamatch works on multivariate numeric data with tiebreaker", {
    test_ldamatch_with_tiebreaker(condition, covariates)
})

test_that("ldamatch works on multivariate numeric data with proportions", {
    test_ldamatch_with_proportions(condition, covariates)
})

test_that("ldamatch works on multivariate numeric data with DX preferences", {
    test_ldamatch_with_DX_preferences(condition, covariates)
})

test_that("ldamatch works on multivariate numeric data with all methods", {
    test_ldamatch_methods(condition, covariates)
})

test_that("ldamatch works on multivariate numeric data with all halting tests", {
    test_ldamatch_halting_tests(condition, covariates,
                                halting_tests_for_numbers)
})

test_that("ldamatch works on multivariate factor data with all methods", {
    test_ldamatch_methods(condition, covariates_factor)
})

test_that("ldamatch works on multivariate factor data with all halting tests", {
    test_ldamatch_halting_tests(condition, covariates_factor,
                                halting_tests_for_factors)
})

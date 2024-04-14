## Globals
.ldamatch_globals <- new.env(parent = emptyenv())

## Default number of replicates for random search.
assign("RND_DEFAULT_REPLICATES", 1e4, .ldamatch_globals)

## Anderson-Darling test parameters; see kSamples::ad.test.
## AD_VERSION: 1 or 2 for the two versions of the AD test statistic
assign("AD_METHOD", "asymptotic", .ldamatch_globals)
assign("AD_NSIM", 10000, .ldamatch_globals)
assign("AD_VERSION", 1, .ldamatch_globals)

## Default setting for printing additional info.
assign("PRINT_INFO", TRUE, .ldamatch_globals)

## Whether to print progress information about parallelly processing cases.
assign("PRINT_PROGRESS", TRUE, .ldamatch_globals)

## The number of cases to be retrieved at a time from iterators for parallel processing.
assign("PROCESSED_CHUNK_SIZE", 10000, .ldamatch_globals)


#' Gets value for ldamatch global parameter.
#' @seealso \code{\link{set_param}} for parameter names.
#'
#' @param name   The name of the global parameter.
#' @return The value of the global parameter.
#'
#' @importFrom RUnit checkTrue
#'
#' @export
get_param <- function(name) {
    RUnit::checkTrue(
        name %in% ls(.ldamatch_globals),
        paste0("Unknown global parameter name ", name, "; chose one of ",
               paste(ls(.ldamatch_globals), collapse = ", ")))
    get(name, .ldamatch_globals)
}


#' Sets value for ldamatch global parameter.
#'
#' @param name   The name of the global parameter.
#' @param value  The new value of the global parameter.
#' @return The previous value of the global parameter.
#'
#' @details The names of the available parameters:
#' \describe{
#'   \item{RND_DEFAULT_REPLICATES}{random search: default number of replicates}
#'   \item{Anderson-Darling test parameters; see kSamples::ad.test for explanation}{
#'     \describe{
#'       \item{AD_METHOD}{the method parameter for ad.test; default: asymptotic}
#'       \item{AD_NSIM}{the Nsim parameter for ad.test, used when AD_METHOD is 'simulated'; default: 10000}
#'       \item{AD_VERSION}{1 or 2 for the two versions of the test statistic; default: 1}
#'     }
#'   }
#'   \item{PRINT_INFO}{print summary information, and progress information for
#'     the exhaustive search algorithm}
#'   \item{PRINT_PROGRESS}{whether to print progress information about parallel processing of cases}
#'   \item{PROCESSED_CHUNK_SIZE}{the number of cases to be retrieved at a time from iterators for parallel processing}
#' }
#'
#' @seealso \code{\link{get_param}} for retrieving the current value of a
#' parameter.
#'
#' @importFrom RUnit checkTrue
#'
#' @export
set_param <- function(name, value) {
    RUnit::checkTrue(
        name %in% ls(.ldamatch_globals),
        paste0("Unknown global parameter name ", name, "; chose one of ",
               paste(ls(.ldamatch_globals), collapse = ", ")))
    prev_value <- get(name, .ldamatch_globals)
    assign(name, value, .ldamatch_globals)
    prev_value
}

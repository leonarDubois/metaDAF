#' Generate several dataset for DAF analysis
#'
#' Call the \code{\link{build_DAF_data}} function several times and aggregate
#' all results ina list. Allow repeated and compared-repeated analysis
#'
#' @inheritParams build_DAF_data
#' @param times Number of repetitions of the DAF analysis.
#' @return a list, each element being the output of a build_DAF_data function call
#'
#' @export



generate_repeated_DAF_data <- function ( count_table, metadata,
                                         times = 10, ...) {

    aggr_res <- list()

    for (i in 1:times) {
        cat("Simulation of DAF dataset", i)
        dt <- build_DAF_data(count_table = count_table,
                            metadata = metadata, ...)
        aggr_res[[i]] <- dt
        cat(" - Done \n")
    }

    return(aggr_res)
}

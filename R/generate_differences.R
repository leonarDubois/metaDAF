#' Generate artificial dataset with differences between 2 groups
#'
#' Use the \code{\link{EDDA::generateData}} function to generate artificial count values table
#'
#' @param count_table a matrix. Used as the \code{inputCount} paramter of the
#' \code{\link{EDDA::generateData}} function
#' @param ... additional parameters for the \code{\link{EDDA::generateData}} function
#'
#' @return a list of  4 elements :
#'   \item{count_table}{ a matrix of count values.}
#'   \item{metadata}{ a data frame of at least one column countaining the group variable.}
#'   \item{true_DAF}{ a character vector. The names of the true DAF.}
#'   \item{FC_true_DAF}{ a numeric vector. The fold-change of the true DAF.}
#'
#' @import EDDA
#' @export

generate_differences <- function(count_table, ... ) {

    if( !require("EDDA")) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("EDDA", version = "3.8")
    }
    library("EDDA")

    dots_arg <- list(...)
    # test on dots_arg ?

    # example to generate data with the generateData function from the EDDA pkg
    # problem can arise if FC parameter is set to default (Norm(2,1)) that leads
    # to diffAbun Features with FC of 1 (log2FC of 0)

    if ("FC" %in% names(dots_arg)) {
        sim_dt <- EDDA::generateData(inputCount = count_table, ...)
        # while (any(sim_dt$DiffAbundList$FC == 1)) {
        #     sim_dt <- EDDA::generateData(inputCount = count_table, ...)
        # }
    } else {
        sim_dt <- EDDA::generateData(inputCount = count_table,
                                    FC = "Norm(5,1)", ...)
        # while (any(sim_dt$DiffAbundList$FC == 1)) {
        #     sim_dt <- EDDA::generateData(inputCount = count_table,
        #                                 FC = "Norm(5,1)", ...)
        # }
    }

    count_tb <- sim_dt$count
    colnames(count_tb) <- as.character(seq(1, ncol(count_tb)))
    metadata <- data.frame(group = sim_dt$dataLabel)

    data_res = list(count_table = count_tb,
                    metadata = metadata)
    data_res$true_DAF <- sim_dt$DiffAbundList$geneName
    data_res$FC_true_DAF <- sim_dt$DiffAbundList$FC

    return(data_res)
}

#' Process the DAF analysis through the DESeq2 package
#'
#' @param data the ouput of the \code{\link{build_DAF_data}} function
#' @param ... additionnal parameters of the method
#' @return the output of the DESeq() function processed through the results() function
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#'

process_DESeq2 <- function( data, ... ) {

    if (!require(DESeq2)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("DESeq2", version = "3.8" )
    }
    library(DESeq2)

    analysis <- DESeq2::DESeqDataSetFromMatrix(
        countData = data$count_table,
        colData   = data$metadata,
        design    = ~ group)


    # -----------------------------------------------------
    # Call the DESeq() function which is a wrapper for :
    # estimateSizeFactors() then
    # estimateDispersions() then
    # Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest()
    res <- do.call(DESeq2::DESeq, c(object = analysis, ...))
    # -----------------------------------------------------


    raw_output <- DESeq2::results(res)

    curated_output <- raw_output[, c("pvalue", "padj")]

    OUT <- list(raw = raw_output,
               curated = as.data.frame(curated_output))

    return(OUT)
}

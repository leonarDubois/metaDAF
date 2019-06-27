#' Process the DAF analysis through the ZIBseq package
#'
#' @inheritParams process_DESeq2
#' @return the rawouput and processed information of DAF
#'

process_ZIBseq <- function( data, ... ) {

    if (!require(ZIBseq)) {
        install.packages("ZIBseq")
    }
    library(ZIBseq)

    # the function is a bit longer than usual because the package
    # is not so easy to deal with.

    dt_used <- t(data$count_table)
    features_kept <- which(colSums(dt_used) > 2 * dim(dt_used)[1])

    raw_output <- ZIBseq(data = t(data$count_table),
                         outcome = data$metadata$group, ...)

    curated_output <- data.frame( pval = raw_output$pvalues)
    rownames(curated_output) <- names(features_kept)
    curated_output$adj.pval <- p.adjust(curated_output$pval, method = "BH")

    # Add features left behind at the beginning
    not_used <- setdiff(rownames(data$count_table), names(features_kept))

    if (length(not_used) > 0) {
        dt_not_used <- data.frame(pval = rep(1, times = length(not_used)))
        rownames(dt_not_used) <- not_used
        dt_not_used$adj.pval <- 1
        curated_output <- rbind(curated_output, dt_not_used)
    }
    OUT <- list(raw = raw_output,
                curated = as.data.frame(curated_output))

    return(OUT)
}

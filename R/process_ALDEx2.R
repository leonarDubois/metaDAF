#' Process the DAF analysis through the ALDEx2 package
#'
#' @inheritParams process_DESeq2
#' @return a list countaining the raw output of ALDEx2 analysis and a curated version
#'
process_ALDEx2 <- function(data, ... ) {

    if (!require(ALDEx2)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("ALDEx2", version = "3.8")
    }
    library(ALDEx2)

    param <- list(...)

    if ("test" %in% names(param)) {
        data_clr <- aldex.clr(reads = data$count_table,
                              conditions = data$metadata$group,
                              mc.samples = 1000)

        if (param$test == "t") {
            raw_output <- aldex.ttest(data_clr,
                        conditions = as.character(data$metadata$group), ...)
        } else if (param$test == "glm") {
            raw_output <- aldex.glm(data_clr,
                        conditions = as.character(data$metadata$group), ...)
        }
    } else {
        raw_output <- aldex(reads = data$count_table,
                        conditions = as.character(data$metadata$group), ...)
    }

    keep <- grep(pattern = "[.ep|.eBH]$", x = colnames(raw_output) )
    # Keep the 2 last matchs, correspond to wilcoxon or glm test pval + BH pval
    curated_output <- raw_output[ , tail(keep, 2) ]

    OUT <- list(raw = raw_output,
               curated = curated_output)
    return(OUT)
}

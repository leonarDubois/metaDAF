#' Process the DAF analysis through the mbzinb package
#'
#' @inheritParams process_DESeq2
#' @return a list countaining the raw output of mbzinb analysis and a curated version



process_mbzinb <- function(data, ... ) {


    if(!require(mbzinb)){
        cat("Installing mbzinb from github repo.")
        install.packages(repos=NULL, method="libcurl", dependencies = T,
                         "https://raw.githubusercontent.com/jchen1981/MicrobiomeDDA/master/mbzinb_0.2.tar.gz")
    }
    library(mbzinb)


    # -----------------------------------------------------
    analysis <- mbzinb.dataset(count = data$count_table, sample = data$metadata)
    fit <- mbzinb.test(mbzinb.data = analysis, group = "group", ...)

    raw_output <- mbzinb.results(fit, nreturn = nrow(data$count_table))

    # -----------------------------------------------------
    curated_output <- raw_output[, c("PValue", "Padj")]

    OUT <- list(raw = raw_output,
               curated = curated_output)

    return(OUT)
}

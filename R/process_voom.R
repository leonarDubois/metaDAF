#' Process the DAF analysis through the vomm method of the limma package
#'
#' @inheritParams process_DESeq2
#' @return a list countaining the raw output of voom analysis and a curated version
#'


process_voom <- function(data, ... ) {

    if (!require(limma)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("limma", version = "3.8")
    }
    library(limma)

    analysis <- DGEList(counts = data$count_table,
                        group = data$metadata$group, remove.zeros = T)
    analysis <- edgeR::calcNormFactors(analysis)
    mod_mat <- model.matrix(~ group, data$metadata)
    mod_voom <- voom(analysis, mod_mat, plot = F)

    fit <- lmFit(mod_voom, mod_mat)
    tmp <- eBayes(fit)

    # Provides log2FC, lfcSE, pvalue, adj-pvalue...
    raw_output <- topTable(tmp, sort.by = "P", n = Inf)
    raw_output <- as.data.frame(raw_output)

    curated_output <- raw_output[, c("P.Value", "adj.P.Val")]

    OUT <- list(raw = raw_output,
                curated = curated_output)

    return(OUT)
}

#' Process the DAF analysis through the edgeR package
#'
#' @inheritParams process_DESeq2
#' @return a list countaining the raw output of edgeR analysis and a curated version
#'


process_edgeR <- function(data, ... ) {

    if (!require(edgeR)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("edgeR", version = "3.8")
    }
    library(edgeR)

    analysis <- DGEList(counts = data$count_table, group = data$metadata$group)

    # Compute Normalization and dipersion estimate
    analysis <- edgeR::calcNormFactors(analysis)
    mod_mat <- model.matrix(~ group, data$metadata)
    analysis <- estimateDisp(analysis, mod_mat)

    # Actual statistics, quasi-likelihood F test
    fit <- glmQLFit(analysis, mod_mat)
    qlf <- glmQLFTest(fit, coef = 2)

    # Provides log2FC, lfcSE, pvalue, adj-pvalue...
    raw_output <- topTags(qlf, n = nrow(data$count_table))
    raw_output <- as.data.frame(raw_output)

    curated_output <- raw_output[, c("PValue", "PValue")]
    curated_output$PValue.1 <- p.adjust(curated_output$PValue.1, method = "BH")

    OUT <- list(raw = raw_output,
               curated = curated_output)

    return(OUT)
}

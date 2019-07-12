#' Process the DAF analysis through the metagenomeSeq package
#'
#' @inheritParams process_DESeq2
#' @return a list countaining the raw output of metagenomeSeq analysis and a curated version
#'

process_metagenomeSeq <- function(data, ... ) {

    if (!require(metagenomeSeq)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("metagenomeSeq", version = "3.8")
    }
    library(metagenomeSeq)

    analysis <- metagenomeSeq::newMRexperiment(counts = data$count_table,
                           phenoData = AnnotatedDataFrame(data$metadata))

    # Compute scaling factors for normalization
    p <- metagenomeSeq::cumNormStatFast(analysis)
    scaling_fct <- metagenomeSeq::cumNorm(analysis, p = p)
    pd <- pData(scaling_fct)
    mod_mat <- model.matrix(~ group, data = pd)

    dots_params <- list(...)
    # fitting Zero-inflated Log-Normal model.
    # ZIG can be fitted but the recent documentation recommends ZIL over ZIG.
    if ("covar" %in% names(dots_params)) {
            covar <- dots_params$covar
            # According to the vignette, fitZig() is relevant for
            # testing confounders and covariates
            mod_mat <- model.matrix(as.formula(paste("~ group", covar,
                                                     collapse = "+")),
                                    data = pd)
            fit <- metagenomeSeq::fitZig(scaling_fct, mod_mat)

    } else {
        fit <- metagenomeSeq::fitFeatureModel(scaling_fct, mod_mat)
    }
    # -----------------------------------------------------
    raw_output <- metagenomeSeq::MRcoefs(fit, number =  nrow(data$count_table))
    curated_output <- raw_output[, c("pvalues", "adjPvalues")]
    OUT <- list(raw = raw_output,
               curated = curated_output)

    return(OUT)
}

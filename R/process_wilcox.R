#' Process the DAF analysis through a wilcoxon test
#'
#' @param data the ouput of the \code{\link{build_DAF_data}} function
#' @param ... additionnal parameters of the method
#' @return the output of the wilcox test for each feature
#'

process_wilcox <- function( data, ... ) {

    raw_output <- data.frame(feat = rownames(data$count_table),
                             pvalue = NA)


    grp <- data$metadata$group[1]

    pvals <- apply(data$count_table, 1, function(x){

        grp1 <- x[which(data$metadata$group == grp)]
        grp2 <- x[which(data$metadata$group != grp)]

        out <- wilcox.test(grp1, grp2, alternative = "two.sided")
        out$p.value
    })

    pvals_adj <- p.adjust(pvals, method = "BH")

    raw_output$pvalue <- pvals
    raw_output$padj <- pvals_adj

    curated_output <- raw_output %>% dplyr::select(-feat)
    rownames(curated_output) <- raw_output$feat

    OUT <- list(raw = raw_output,
                curated = as.data.frame(curated_output))

    return(OUT)
}

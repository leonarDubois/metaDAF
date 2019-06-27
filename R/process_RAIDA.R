#' Process the DAF analysis through the RAIDA package
#'
#' @param data the ouput of the \code{\link{build_DAF_data}} function
#' @param ... additionnal parameters of the method
#' @return the output of the DESeq() function processed through the results() function
#'

process_RAIDA <- function( data, ... ) {

    if (!require(RAIDA)) {
        install.packages(repos=NULL, method="libcurl", dependencies = T,
            "http://cals.arizona.edu/~anling/software/RAIDA_1.0.tar.gz")
    }
    library(RAIDA)

    grp_sizes <- table(data$metadata[ , "group"])
    ord <- order(data$metadata[ , "group"])
    data$metadata <- data$metadata[ord, ]
    data$count_table <- data$count_table[ , ord]

    raw_output <- raida(c.data = as.data.frame(data$count_table),
                      n.lib = grp_sizes, ...)

    curated_output <- raw_output[, c("p", "p.adj")]

    OUT <- list(raw = raw_output,
                curated = as.data.frame(curated_output))

    return(OUT)
}

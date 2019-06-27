#' Filter DAF data
#'
#' Remove some rows and columns from the count data table so that further process
#' is made easier
#'
#' @param data list. the output of the \code{\link{build_DAF_data}} function
#' @return a list with the exact same format as the input data param
#'
#' @export



preprocess_DAF_data <- function(data,
                                taxon_lvl = c("g", "s", "t", "f",
                                              "o", "p", "c", NULL) ) {

    taxon_lvl <- match.arg(taxon_lvl)

    # is nested list ? if yes, means repeated analysis
    if (any(sapply(data[[1]], is.list)) ) {
        data <- lapply(data, function(x) {
            keep <- !duplicated(rownames(x$count_table))
            x$count_table <- x$count_table[keep, ]

            if (!is.null(taxon_lvl)) {
                keep <- grepl(pattern = paste0("^", taxon_lvl,"(__)*[A-Za-z]*"),
                              rownames(x$count_table))
                x$count_table <- x$count_table[keep, ]
            }
            return(x)
        })
    } else {
        # test in order to avoid repetitions in the uncharacterized species
        keep <- !duplicated(rownames(data$count_table))
        data$count_table <- data$count_table[keep, ]

        if (!is.null(taxon_lvl)) {
            keep <- grepl(pattern = paste0("^", taxon_lvl,"(__)*[A-Za-z]*"),
                          rownames(data$count_table))
            data$count_table <- data$count_table[keep, ]
        }
    }

    return(data)
}

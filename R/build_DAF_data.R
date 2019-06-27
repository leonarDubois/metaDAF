#' Load all data for DAF analysis
#'
#' Create a list that will contain all data needed for differential
#' features analysis.
#'
#' @param count_table a matrix of count values with features as rows and samples as columns.
#' @param metadata a data frame. Rows are samples and miscellaneous variables as columns.
#' @param group_var character. The name of \code{metadata} column defining the groups.
#' @param simulation logical. Should differences be generated ? Default is TRUE.
#' @param ... additionnal parameter passed to the \code{\link{generateData}} function
#'
#' @return a list of 2 or 4 elements :
#'   \item{count_table}{ a matrix of count values.}
#'   \item{metadata}{ a data frame of at least one column countaining the group variable.}
#'
#'   If the \code{simulation} parameter was TRUE, artificial differences were
#'   generated. It is possible to know exactly which are the true
#'   differentially abundant features (DAF).
#'
#'   \item{true_DAF}{ a character vector. The names of the true DAF.}
#'   \item{FC_true_DAF}{ a numeric vector. The fold-change of the true DAF.}
#'
#' @details If the columns names of the count table do not match the rows names of the
#' metadata table, metadata columns are screened for a variable countaining these names.
#'
#' If such variable is not found, it is assumed that individuals follow the same order
#' through count table columns and metadata rows.
#'
#' @export


build_DAF_data <- function(count_table, metadata,
                          group_var = NULL, simulation = TRUE, ...){

    # -----------------------------------------------------
    # Test values of arguments
    if (ncol(count_table) != nrow(metadata)) {
        stop(paste0("Incorrect dimensions. \"metadata\" rows and ",
            "\"count_table\" columns must have the same number of individuals"))
    }

    if (!is.null(group_var)){
        if(!group_var %in% colnames(metadata)) {
            stop(paste0("Incorrect name for \"group_var\" argument. ",
                "Must be a column name of metadata"))
        } else {
            # grp_var is not null and in colnames(metadata)
            nb_grp <- length(table(metadata[ , group_var]))
            if (nb_grp < 2) {
                stop(paste0("Incorrect number of groups defined by ",
                    "\"group_var\" argument. Must be >= 2"))
            }
        }
    } else if (!simulation) {
        stop(paste0("Arguments \"group_var\" and \"simulation\" ",
            "are missing. One should be specified."))
    }
    # -----------------------------------------------------

    if (is.null(group_var)) {
        suppressMessages(invisible(capture.output(
            DAF_data <- generate_differences(count_table, ...)
        )))

    } else {
        name_match <- match( rownames(metadata), colnames(count_table))

        if (sum(name_match) == 0 | is.na(sum(name_match))) {
            cat("Count table colnames and metadata rownames do not match \n
            Looking for ID/names matching in metadata columns\n")

            col_match <- apply(metadata, 2, function(x) {
                a <- match(x, colnames(count_table))
                return( !all(is.na(a)) )
            })
            col_match <- names(which(col_match))[1]

            if(!is.na(col_match)){
                cat("Match found in column :", col_match, "\n")
                ord <- match(colnames(count_table),
                            metadata[ , col_match] )

                metadata <- metadata[ord, ]
                rownames(metadata) <- metadata[ , col_match]
            } else {
                cat(" - No match found - \n
                    It is assumed that individuals are in the same order \n")
                rownames(metadata) <- colnames(count_table)
            }
        }

        metadata = metadata %>% dplyr::select(group_var)
        colnames(metadata) <- "group"

        metadata$group <- as.numeric(factor(metadata$group)) -1

        DAF_data = list(count_table = count_table,
                        metadata = metadata)
    }


    return(DAF_data)
}

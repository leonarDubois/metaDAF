#' Compare methods features rank
#'
#' Uses processed results object to order the features based on their p-value
#'
#' @param results curated results from DAF analysis run
#' @param adjusted logical. should adjusted p-value be used ? default is TRUE
#' @return a tibble of features names, method name and rank of the feature pvalue
#'
#' @export
rank_features <- function(results, adjusted = TRUE){

    results <- na.omit(results)


    if (adjusted) {
        results = results %>% group_by(method) %>%
                mutate(rank = rank(adjusted_pvalue)) %>%
                ungroup()
    } else {
        results = results %>% group_by(method) %>%
                mutate(rank = rank(pvalue)) %>%
                ungroup()
    }
    return(results)
}

# -------------------------------------------------------------------

#' Plot the distribution of the rank accross methods
#'
#'
#' @param data the output of the \code{\link{rank_features}} function
#' @param nmax max number of features to plot
#'
#' @return boxplot. ggplot object
#'
#' @export
#'
boxplot_rank <- function(dt, nmax = 30){

    dt_filtered <- dt %>% group_by(feat_name) %>%
        mutate(mean_rank = mean(rank)) %>% ungroup() %>%
        arrange(mean_rank)

    feat_ranks <- dt_filtered %>%
                dplyr::select(feat_name, mean_rank) %>%
                unique(.)

    to_keep <- feat_ranks[1:nmax, "feat_name"]

    dt <- dt %>% filter(feat_name %in% to_keep$feat_name)



    p <- ggplot(data = dt,
                mapping = aes(x = reorder(feat_name, rank, FUN=median),
                              y = rank)) +
        geom_boxplot(outlier.shape = NA, fill = "lightgrey") +
        geom_jitter(mapping = aes(x = reorder(feat_name, rank, FUN=median),
                                  y = rank,
                                  fill = method),
                    shape = 22, size = 2.5,
                    position = position_jitter(height = 0, width = 0.3)) +
        coord_flip() +
        theme_bw() +
        scale_color_manual(values = rainbow(length(unique(dt$method))))
    p
}

#' Plot a ROC or PR curve
#'
#' @param res_data a data frame. The ouput of the \code{\link{generate_curve_data}} function
#' @return a ggplot object with the ROC curve data
#'
#' @import ggplot2 ROC
#' @export

curve_plot <- function( curve_data,
                        type = c("roc", "pr"),
                        plot_title = "[Title]", auc_limit = 1) {
    type <- match.arg(type)

    if (!require("ROC")){
        # install.packages("ROC")
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

        BiocManager::install("ROC")
    }
    library("ROC")
    if (!require("viridis")) {
        install.packages("viridis")
    }
    library("viridis")

    p <- ggplot() +  xlim(0,1) + ylim(0,1) +
        labs(title = plot_title) + theme_bw() +
        theme(plot.title = element_text(size = 18, face = "bold"),
              legend.title = element_text(size=12),
              legend.text = element_text(size=12),
              axis.text = element_text(size=12),
              axis.title = element_text(size=12) )

    if (type == "roc") {
        curve_data <- curve_data$roc
        colnames(curve_data)[which(colnames(curve_data) == "fpr")] <- "x_aes"
        colnames(curve_data)[which(colnames(curve_data) == "tpr")] <- "y_aes"

        p <- p + geom_abline(aes(intercept = 0, slope = 1),
                             linetype = "dashed", alpha = 0.5) +
            labs(x = "False positive rate",
                 y = "True positive rate")
    } else {
        curve_data <- curve_data$pr
        colnames(curve_data)[which(colnames(curve_data) == "tpr")] <- "x_aes"
        ind <- which(colnames(curve_data) == "precision")
        colnames(curve_data)[ind] <- "y_aes"
        p <- p + labs(x = "Recall (TPR)",
                      y = "Precision")
    }

    compare <- length(unique(curve_data$method)) >= 2
    repeated <- length(unique(curve_data$run)) >= 3
    if (!repeated) {
        curve_data <- curve_data %>% filter(run != "run_cat")
    }

    if (compare) {
        if (repeated) {
            # compare repeated loop
            # plot several color curves but avg for each
            p <- plot_compare_rep(curve_data, p, auc_limit)
        } else {
            # compare non repeated loop
            # plot several color curves
            p <- plot_compare(curve_data, p, auc_limit)
        }
    } else {
        if (repeated) {
            # non compare repeated loop
            # plot all repetitions as grey thin lines
            p <- plot_rep(curve_data, p, auc_limit)
        } else {
            # non compare non repeated loop
            # plot one curve
            p <- plot_simple_curve(curve_data, p, auc_limit)
        }
    }
    p
}

# -------------------------------------------------------------------


plot_simple_curve <- function (curve_data, p, auc_limit) {

    auc <- trapezint(curve_data$x_aes, curve_data$y_aes,
                     a = 0, b = auc_limit) / auc_limit
    p <- p + geom_line(data = curve_data,
                       mapping = aes(x = x_aes, y = y_aes),
                       size = 1.2, color = "black") +
        labs(subtitle = paste("AUC", auc_limit, "=",
                              as.character(round(auc, digits = 3)) ))
    return(p)
}

# -------------------------------------------------------------------

plot_rep <- function (curve_data, p, auc_limit) {

    repeat_run <- setdiff(unique(curve_data$run), "run_cat")
    for (j in repeat_run ) {
        tmp <- curve_data %>% filter(run == j)
        p <- p + geom_line(data = tmp,
                           mapping = aes(x = x_aes, y = y_aes),
                           size = 0.8, color = "grey")
    }

    curve_data <- curve_data %>% filter(run == "run_cat")

    p <- plot_simple_curve(curve_data, p, auc_limit)
    return(p)
}

# -------------------------------------------------------------------

plot_compare <- function (curve_data, p, auc_limit) {

    for (m in unique(curve_data$method)) {
        tmp <- curve_data %>% filter(method == m)
        auc <- trapezint(tmp$x_aes, tmp$y_aes,
                         a = 0, b = auc_limit) / auc_limit

        new <- paste(m, "AUC", auc_limit, "=",
                     as.character(round(auc, digits = 3)) )
        curve_data[which(curve_data$method == m), "method"] <- new
    }
    p <-  p + geom_line(data = curve_data,
                        mapping = aes(x = x_aes, y = y_aes, color = method),
                        size = 1.4) +
        scale_color_manual(values = cividis(length(unique(curve_data$method))))

    return(p)
}

# -------------------------------------------------------------------

plot_compare_rep <- function (curve_data, p, auc_limit) {

    curve_data <- curve_data %>% filter(run == "run_cat")

    for (m in unique(curve_data$method)) {
        tmp <- curve_data %>% filter(method == m)
        auc <- trapezint(tmp$x_aes, tmp$y_aes,
                         a = 0, b = auc_limit) / auc_limit
        new <- paste(m, "AUC", auc_limit, "=",
                     as.character(round(auc, digits = 3)) )
        curve_data[which(curve_data$method == m), "method"] <- new
    }

    p <- p + geom_line(data = curve_data,
                        mapping = aes(x = x_aes, y = y_aes, color = method),
                        size = 1.4) +
    scale_color_manual(values = cividis(length(unique(curve_data$method))))

    return(p)
}

# -------------------------------------------------------------------

#' Boxplot of AUC distribution across repetitions
#'
#' @param curve_data the ouput of the \code{\link{generate_curve_data}} function
#' @return a ggplot object with the boxplot data
#'
#' @import ggplot2 ROC
#'
#' @export

boxplot_AUC <- function (curve_data, auc_limit = 1) {

    curve_data <- curve_data$roc
    curve_data <- curve_data %>% filter(run != "run_cat")

    all_auc <- NULL
    for (m in unique(curve_data$method)) {
        tmp <- curve_data %>% filter(method == m)
        for (j in unique(curve_data$run) ) {
            tmp2 <- tmp %>% filter(run == j)

            auc <- trapezint(tmp2$fpr, tmp2$tpr,
                             a = 0, b = auc_limit) / auc_limit
            all_auc <- rbind(all_auc, data.frame(auc = auc, method = m))
        }
    }

    nb_color <- length(unique(curve_data$method))

    p <- ggplot(all_auc, aes(x = method, y = auc)) +
        geom_boxplot(outlier.shape = NA, fill = rainbow(nb_color), lwd = 1) +
        geom_jitter(shape = 21, bg = "grey", col = "black",
                    size = 3 , position = position_jitter(0.2)) +
        theme_bw()
    p
}

#' Compute metrics needed for ROC and PR curves
#'
#' Based on results from differentially abundant features (DAF) analysis and a
#' variable used for threshold definition, it will compute 3 differents metrics,
#' allowing to plot ROC and precision-recall curves.
#'
#' @param dt a data frame. The output of the \code{\link{run_DAF_analysis}} function.
#'  If the function return a list countaining both "raw" and "curated" results,
#'  then curated results must be used.
#' @param threshold_var a character. The name of the column that should
#'  be used as threshold variable. By default "adjusted_pvalue".
#'
#' @return a list of 2 data frame :
#' \item{roc}{Countains TPR, FPR, method name and run number for each threshold value}
#' \item{pr}{Countains TPR, precision, method name and run number for each threshold value}
#'
#' @export

generate_curve_data <- function(dt, threshold_var = "adjusted_pvalue") {

    dt <- na.omit(dt)
    dt <- unique(dt)
    all_roc <- c()
    all_pr <- c()

    for (m in unique(dt$method) ) {

        dt_tmp <- dt %>% filter(method == m)

        for (r in c(unique(dt_tmp$run), "run_cat") ) {

            if (r == "run_cat") {
                dt_tmp2 <- dt_tmp
            } else {
                dt_tmp2 <- dt_tmp %>% filter(run == r)
            }

            thres_v <- c( sort(dt_tmp2[, threshold_var]))

            TP <- unlist(lapply(thres_v, function(x) {
                sum(dt_tmp2[, threshold_var] <= x & dt_tmp2$true_DAF)
            }))
            FP <- unlist(lapply(thres_v, function(x) {
                sum(dt_tmp2[, threshold_var] <= x & !dt_tmp2$true_DAF)
            }))
            TN <- unlist(lapply(thres_v, function(x) {
                sum(dt_tmp2[, threshold_var] > x & !dt_tmp2$true_DAF)
            }))
            FN <- unlist(lapply(thres_v, function(x) {
                sum(dt_tmp2[, threshold_var] > x & dt_tmp2$true_DAF)
            }))

            tpr <- TP / (TP + FN)
            fpr <- FP / (FP + TN)
            precision <- TP / (TP + FP)
            out_roc <- data.frame(tpr = tpr, fpr = fpr)
            out_pr <- data.frame(tpr = tpr, precision = precision)

            out_roc <- out_roc %>% group_by(fpr) %>%
                mutate(tpr = min(tpr)) %>%
                ungroup() %>% unique()
            out_pr <- out_pr %>% group_by(precision) %>%
                mutate(tpr = min(tpr)) %>%
                ungroup() %>% unique()

            out_roc$method <- dt_tmp2$method[1]
            out_roc$run <- r
            out_pr$method <- dt_tmp2$method[1]
            out_pr$run <- r

            all_roc <- rbind(all_roc, out_roc)
            all_pr <- rbind(all_pr, out_pr)
        }
    }

    return(list(roc = all_roc, pr = all_pr))
}

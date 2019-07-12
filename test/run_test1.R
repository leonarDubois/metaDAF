rm(list = ls())
library(devtools)
library(tidyverse)
devtools::load_all()

# ------------------------------------------------------------------------

devtools::document()
library(pkgdown)
pkgdown::build_site()

# ------------------------------------------------------------------------
# test without simulation

dataTest <- build_DAF_data(count_table = ZellerG_2014_table_count,
                           metadata = ZellerG_2014_metadata,
                           group_var = "study_condition",
                           simulation = NULL)

# ------------------------------------------------------------------------

if ( !require("EDDA")) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("EDDA", version = "3.8")
}
library("EDDA")


dataTest <- build_DAF_data(count_table = HMP_2012_table_count,
                           metadata    = HMP_2012_metadata,
                           group_var   = NULL,
                           simulation  = TRUE, ControlRep = 5,
                           numDataPoints = 20000)
# Then analyze
res <- run_DAF_analysis(data = dataTest, method = "deseq2")
plot_data <- generate_curve_data(dt = res$curated,
                                 threshold_var = "adjusted_pvalue")
curve_plot(curve_data = plot_data, type = "roc", plot_title = "ROC curve DESeq2")


# --------------------------------
# Test repetitions
# --------------------------------
dataTest_double <- generate_repeated_DAF_data(times     = 5,
                                            count_table = HMP_2012_table_count,
                                            metadata    = HMP_2012_metadata,
                                            ControlRep = 2,
                                            numDataPoints = 20000)
res2 <- run_DAF_analysis(data = dataTest_double, method = "deseq2")
plot_data2 <- generate_curve_data(dt = res2,
                                 threshold_var = "adjusted_pvalue")
curve_plot(curve_data = plot_data2, type = "roc", plot_title = "ROC curve DESeq2 repeated")


# --------------------------------
# Test comparisons
# --------------------------------
dataTest <- build_DAF_data(count_table = QinJ_2012_table_count,
                           metadata    = QinJ_2012_metadata,
                           group_var   = NULL,
                           simulation  = TRUE,
                           ControlRep = 5,
                           numDataPoints = 20000)
res <- run_DAF_analysis(data = dataTest, method = c("deseq2", "mbzinb",
                                                     "raida",
                                                    "metagenomeseq", "voom"))
plot_data <- generate_curve_data(dt = res,
                                 threshold_var = "adjusted_pvalue")
curve_plot(curve_data = plot_data, type = "roc", plot_title = "ROC curve Comparison")
curve_plot(curve_data = plot_data, type = "roc") + coord_cartesian(xlim = c(0,0.1))

# and test rank feature
boxplot_rank(rank_features(res, adjusted = FALSE), nmax = 10)

# --------------------------------
# Test comparisons AND repetitions
# --------------------------------
dataTest <- generate_repeated_DAF_data(times       = 3,
                                       count_table = HMP_2012_table_count,
                                       metadata    = HMP_2012_metadata,
                                       ControlRep = 5,
                                       numDataPoints = 20000)
res <- run_DAF_analysis(data = dataTest, method = c("deseq2", "mbzinb",
                                                    "raida",
                                                    "metagenomeseq", "voom"))
plot_data <- generate_curve_data(dt = res,
                                 threshold_var = "adjusted_pvalue")
curve_plot(curve_data = plot_data, type = "roc", plot_title = "ROC curve repeat and compare")

boxplot_AUC(curve_data = plot_data, auc_limit = 1)


# --------------------------------
# Test boxplot focused on AUC
# --------------------------------
dataTest <- generate_repeated_DAF_data(times = 10,
                                       count_table = HMP_2012_table_count,
                                       metadata = HMP_2012_metadata,
                                       ControlRep = 10,
                                       numDataPoints = 20000)
res <- run_DAF_analysis(data = dataTest, method = c("deseq2", "mbzinb",
                                                    "raida",
                                                    "metagenomeseq", "voom",
                                                    "wilcox", "edger"))

plot_data <- generate_curve_data(dt = res, threshold_var = "adjusted_pvalue")
boxplot_AUC(curve_data = plot_data, auc_limit = 1)

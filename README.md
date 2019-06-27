README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
TODO
----

Overall improvement that are yet to be added :

-   Filter ROC curves and results based on the FC of the feature across groups (results for low, med, high FC)
-   Add messages printed for user so that he knows what's going on (use cat())
-   install packages with remotes::install\_github()
-   Add details about method (stats, normalization...) in the roxygen (in fact add what was in the first methods.Rmd in the doc of the package)

-   ROC aucs comparable if methods include a different pvalue correction ?

metaDAF
=======

The goal of metaDAF is to wrap packages and functions for the statistical analysis of diffenrentially abundant features. If data does not already countain groups to compare, functions are provided to generate artificial differences and thus benchmark methods included.

Installation
------------

see later

Example
-------

Here is the same analysis repeated several times :

``` r
dataTest_double <- generate_repeated_DAF_data(times = 5,
                                              count_table = QinJ_2012_table_count,
                                              metadata = QinJ_2012_metadata)
res2 <- run_DAF_analysis(data = dataTest_double, method = "deseq2")
plot_data2 <- generate_curve_data(dt = res2,
                                 threshold_var = "adjusted_pvalue")
curve_plot(curve_data = plot_data2, type = "roc",
           plot_title = "Example of DESeq2 run - 5 times")
```

<img src="man/figures/README-example2-1.png" width="100%" />

Or if one want to compare multiple methods together :

``` r
curve_plot(curve_data = res_test_compare_rep_plot_data, type = "roc",
           plot_title = "Example of methods comparison")
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

``` r
boxplot_AUC(curve_data = res_test_compare_rep_plot_data, auc_limit = 1)
```

<img src="man/figures/README-unnamed-chunk-1-2.png" width="100%" />

``` r


curve_plot(curve_data = res_test_compare_rep_plot_data, type = "roc",
           plot_title = "Example of methods comparison") + coord_cartesian(xlim = c(0, 0.1))
```

<img src="man/figures/README-unnamed-chunk-1-3.png" width="100%" />

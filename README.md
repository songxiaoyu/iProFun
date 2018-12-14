
<!-- README.md is generated from README.Rmd. Please edit that file -->
`iProFun`
=========

An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA copy number alterations (CNA) and DNA methylation
---------------------------------------------------------------------------------------------------------------------------------------------------

The goal of **iProFun** is to

-   characterize functional consequences of DNA copy number and methylation alterations in tumors

-   facilitate screening for cancer drivers contributing to tumor initiation and progression, since CNAs and DNA methylations that preserve functional consequences are more likely to be cancer drivers.

### Installation

You can install the latest development version from GitHub with

``` r
install.packages("devtools")
devtools::install_github("xiaoyu/iProFun")
```

<!-- * the most recent officially-released version from CRAN with -->
<!--     ```R -->
<!--     install.packages("iProFun") -->
<!--     ```` -->
<!-- * the latest development version from GitHub with -->
<!--     ```R -->
<!--     install.packages("devtools") -->
<!--     devtools::install_github("xiaoyu/iProFun") -->
<!--     ```` -->
iProFun Integrative analysis pipeline
-------------------------------------

Below is an example of how iProFun is commonly used. A full description of the tool can be found in our MCP paper.

``` r
library(iProFun)
```

### Data summary

After preprossing and data cleaning, we have 15121 genes and 569 subjects for mRNA data, 7010 genes and 174 subjects for protein, 5685 genes and 70 subjects for phospho data, 25762 genes and 552 subjects for methylation data, 11859 genes and 560 subjects for mRNA data. The following shows the data structure for each data.

``` r
rna_normal[1:5,1:5]
```

    ##   Gene_ID TCGA-04-1331 TCGA-04-1332 TCGA-04-1335 TCGA-04-1336
    ## 1    A1BG    1.1319399   -0.9080570  -0.56521708  -1.44252482
    ## 2     A2M    0.4048469    1.8307569   1.00921057   0.36595994
    ## 3   A2ML1    1.1122897    0.8220087  -0.60655790  -0.51220277
    ## 4  A4GALT    0.7822480    0.4631650   2.21204483   0.03091712
    ## 5   A4GNT    0.3174587   -1.3253512  -0.09567758   0.46415113

### Gene-level multiple linear regression to obtain summary statistics

We use sets of separate regressions in the integrative analysis pipeline to allow for different samples being measured on different sets of molecular features.

``` r
ylist_normal = list(rna_normal, protein_normal, phospho_normal)
methy_input_1_3 <-
MultiReg_together(
ylist = ylist_normal,
xlist = list(methy, cnv),
covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
cl = cl
)
```

The following shows the results for CNA.

``` r
str(cnv_input_1_3)
```

    ## List of 7
    ##  $ betas_J   : num [1:676, 1:3] 2.15 1.31 1.71 2.3 1.18 ...
    ##  $ betas_se_J: num [1:676, 1:3] 0.136 0.131 0.116 0.122 0.106 ...
    ##  $ sigma2_J  : num [1:676, 1:3] 0.687 0.824 0.607 0.596 0.814 ...
    ##  $ dfs_J     : int [1:676, 1:3] 518 519 518 519 519 518 518 519 518 519 ...
    ##  $ v_g_J     : num [1:676, 1:3] 0.0268 0.0209 0.0223 0.025 0.0139 ...
    ##  $ xName     :Classes 'data.table' and 'data.frame': 676 obs. of  1 variable:
    ##   ..$ Gene_ID: chr [1:676] "AAAS" "ABCC1" "ABI1" "ABI2" ...
    ##   ..- attr(*, ".internal.selfref")=<externalptr> 
    ##  $ yName     : chr [1:676, 1] "AAAS.NP_001166937.1:s462" "ABCC1.NP_004987.2:s930" "ABI1.NP_001012768.1:s222" "ABI2.NP_001269854.1:s183" ...

### Primo – An integrative analysis method for detecting joint associations of DNA al- terations with multi-omics traits

With the summary association statistics obtained from equations (1), we apply an integrative analysis method – Primo – to detect joint associations of DNA variation with multi-omics traits

``` r
pi1=rep(0.05, 3)
cnv_1_3 = MultiOmics_Input(cnv_input_1_3, pi1=pi1)
```

``` r
str(cnv_1_3)
```

    ## List of 7
    ##  $ NoComputation: int 676
    ##  $ Config       : num [1:8, 1:3] 0 1 0 0 1 1 0 1 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:8] "Q" "" "" "" ...
    ##   .. ..$ : NULL
    ##  $ PostProb     :'data.frame':   676 obs. of  9 variables:
    ##   ..$ Gene_ID: chr [1:676] "AAAS" "ABCC1" "ABI1" "ABI2" ...
    ##   ..$ 1      : num [1:676] 3.85e-47 2.12e-24 1.64e-41 4.48e-65 3.36e-32 ...
    ##   ..$ 2      : num [1:676] 1.07e-02 8.82e-05 4.78e-03 9.55e-07 1.93e-08 ...
    ##   ..$ 3      : num [1:676] 3.33e-56 3.70e-31 6.18e-50 1.91e-70 1.96e-38 ...
    ##   ..$ 4      : num [1:676] 1.68e-52 3.61e-30 1.80e-47 5.96e-70 1.91e-34 ...
    ##   ..$ 5      : num [1:676] 0.355673 0.590296 0.688006 0.155601 0.000432 ...
    ##   ..$ 6      : num [1:676] 1.58e-04 5.04e-07 1.76e-05 4.27e-08 3.68e-07 ...
    ##   ..$ 7      : num [1:676] 1.95e-59 8.45e-35 9.08e-54 3.40e-73 1.50e-38 ...
    ##   ..$ 8      : num [1:676] 0.633 0.41 0.307 0.844 1 ...
    ##  $ colocProb    : num [1:8] 1.70e-02 7.31e-02 7.26e-13 3.42e-08 3.77e-01 ...
    ##  $ Tstat_L      : num [1:676, 1:3] 15.8 10.1 14.6 18.7 11.2 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:3] "moderate.t" "moderate.t" "moderate.t"
    ##  $ D0           : num [1:676, 1:3] 2.34e-46 1.96e-21 2.38e-40 2.55e-60 1.19e-25 ...
    ##  $ D1           : num [1:676, 1:3] 0.0153 0.0191 0.0162 0.0127 0.016 ...

### False discovery rate assessment

To calculate the empirical FDR, we first calculated the posterior probability of a predictor being associated with an outcome, by summing over all patterns that are consistent with the association of interest.

The following shows the results when we randomly permute the sample label of the mRNA while keeping the labels of the other two traits.

``` r
 MultiReg_cnv_lr_perm_1 = MultiReg_together_perm(
    ylist = list(rna_regression, protein_regression, phospho_regression),
    xlist = list(cnv_lr_regression, cnv_baf_regression, methy_mean_regression),
    covariates = list(purity_tumor,age, gender),
    xyCommonGeneID = xy_common_geneID,
    conditional_covariate = mutation_reg_111,
    mutation_genes = mutation_gene_111,
    xyCommonSubID = list(xrnaCommonSubID, xproteinCommonSubID, xphosphoCommonSubID),
    filename = "MultiReg_cnv_lr_together_perm_1",
    permcolum = 1,
    seed=(currind-1)*10+i
  )
```

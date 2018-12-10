
<!-- README.md is generated from README.Rmd. Please edit that file -->
`iProFun`
=========

An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA copy number alterations (CNA) and DNA methylation
---------------------------------------------------------------------------------------------------------------------------------------------------

The goal of **iProFun** is to

-   characterize functional consequences of DNA copy number and methylation alterations in tumors

-   facilitate screening for cancer drivers contributing to tumor initiation and progression, since CNAs and DNA methylations that preserve functional consequences are more likely to be cancer drivers.

### Installation

You can install:

<!-- * the most recent officially-released version from CRAN with -->
<!--     ```R -->
<!--     install.packages("iProFun") -->
<!--     ```` -->
-   the latest development version from GitHub with

    ``` r
    install.packages("devtools")
    devtools::install_github("xiaoyu/iProFun")
    ```

iProFun Integrative analysis pipeline
-------------------------------------

Below is an example of how iProFun is commonly used. A full description of the tool can be found in our MCP paper.

``` r
library(iProFun)
require(metRology)
```

    ## Loading required package: metRology

    ## 
    ## Attaching package: 'metRology'

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, rbind

``` r
require(matrixStats)
```

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

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

``` r
protein_normal[1:5,1:5]
```

    ##   Gene_ID TCGA-09-1664 TCGA-13-1484 TCGA-13-1488 TCGA-13-1489
    ## 1    A1BG    1.0657701   -1.5328933    0.1679014   0.36195190
    ## 2     A2M    0.3074049   -0.4732287    0.6293645  -0.44578651
    ## 3    AAAS   -0.5805011   -3.2799496    0.2539013  -0.19987182
    ## 4    AACS   -0.1564924    0.2862237   -0.3579192  -0.94803135
    ## 5    AAK1   -1.2095157   -0.4839372   -0.2516493   0.09203467

``` r
phospho_normal[1:5, 1:5]
```

    ##   Gene_ID                phospho_ID TCGA-13-1484 TCGA-13-1489 TCGA-13-1494
    ## 1    AAAS  AAAS.NP_001166937.1:s462    0.6523611   -0.6446869    -1.258875
    ## 2    AAK1     AAK1.NP_055726.3:t640   -0.4294914           NA           NA
    ## 3   ABCC1    ABCC1.NP_004987.2:s930    1.0693706   -1.4115053    -0.906477
    ## 4   ABCF1 ABCF1.NP_001020262.1:s109    0.7664572   -2.5787145    -2.726857
    ## 5   ABCF1 ABCF1.NP_001020262.1:s228    1.3439996           NA    -2.047463

``` r
methy[1:5,1:5]
```

    ##         Gene_ID Hybridization chr TCGA-04-1331 TCGA-04-1332
    ## 1        ATP2A1    cg00000292  16    0.9221302   0.46481857
    ## 2         SLMAP    cg00002426   3    0.1486593   0.05626750
    ## 3         MEOX2    cg00003994   7    0.0293567   0.03447924
    ## 4         HOXD3    cg00005847   2    0.8248150   0.42014257
    ## 5 ZNF425;ZNF398    cg00006414   7           NA           NA

``` r
cnv[1:5, 1:5]
```

    ##   Gene_ID TCGA-04-1331 TCGA-04-1332 TCGA-04-1335 TCGA-04-1336
    ## 1     A2M       0.1434       0.2716      -0.5880      -0.6038
    ## 2   A2ML1       0.1434       0.2716      -0.5880      -0.6038
    ## 3 A3GALT2       0.1771       0.0970       0.1854       0.2455
    ## 4  A4GALT      -0.5054      -0.2474       0.1573      -0.5769
    ## 5   A4GNT       0.1777       0.1471       0.6189       0.6611

### Gene-level multiple linear regression to obtain summary statistics

``` r
ylist_normal = list(rna_normal, protein_normal, phospho_normal)
methy_input_1_3 <-
MultiReg_together(
ylist = ylist_normal,
xlist = list(methy, cnv),
covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
cl = cl
)
cnv_input_1_3 <-
MultiReg_together(
ylist = ylist_normal,
xlist = list(cnv, methy),
covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
cl = cl
)
```

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

``` r
str(methy_input_1_3)
```

    ## List of 7
    ##  $ betas_J   : num [1:1103, 1:3] 0.461 -0.425 -6.512 2.547 -1.137 ...
    ##  $ betas_se_J: num [1:1103, 1:3] 0.166 0.978 6.349 1.605 5.993 ...
    ##  $ sigma2_J  : num [1:1103, 1:3] 1.016 0.592 0.592 0.592 0.592 ...
    ##  $ dfs_J     : int [1:1103, 1:3] 519 507 507 507 507 507 507 507 507 507 ...
    ##  $ v_g_J     : num [1:1103, 1:3] 0.0273 1.6159 68.0259 4.3475 60.6253 ...
    ##  $ xName     :Classes 'data.table' and 'data.frame': 1103 obs. of  3 variables:
    ##   ..$ Gene_ID      : chr [1:1103] "EPB41L3" "RB1" "RB1" "RB1" ...
    ##   ..$ Hybridization: chr [1:1103] "cg00027083" "cg00059930" "cg01738359" "cg08383063" ...
    ##   ..$ chr          : chr [1:1103] "18" "13" "13" "13" ...
    ##   ..- attr(*, ".internal.selfref")=<externalptr> 
    ##  $ yName     : chr [1:1103, 1] "EPB41L3.NP_001268463.1:s579" "RB1.NP_000312.2:s788" "RB1.NP_000312.2:s788" "RB1.NP_000312.2:s788" ...

### Primo – An integrative analysis method for detecting joint associations of DNA al- terations with multi-omics traits

``` r
pi1=rep(0.05, 3)
cnv_1_3 = MultiOmics_Input(cnv_input_1_3, pi1=pi1)
methy_1_3 = MultiOmics_Input(methy_input_1_3, pi1=pi1)
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

``` r
str(methy_1_3)
```

    ## List of 7
    ##  $ NoComputation: int 1103
    ##  $ Config       : num [1:8, 1:3] 0 1 0 0 1 1 0 1 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:8] "Q" "" "" "" ...
    ##   .. ..$ : NULL
    ##  $ PostProb     :'data.frame':   1103 obs. of  11 variables:
    ##   ..$ Gene_ID      : chr [1:1103] "EPB41L3" "RB1" "RB1" "RB1" ...
    ##   ..$ Hybridization: chr [1:1103] "cg00027083" "cg00059930" "cg01738359" "cg08383063" ...
    ##   ..$ chr          : chr [1:1103] "18" "13" "13" "13" ...
    ##   ..$ 1            : num [1:1103] 0.918 0.965 0.721 0.598 0.829 ...
    ##   ..$ 2            : num [1:1103] 0.0669 0.0118 0.0835 0.0372 0.0564 ...
    ##   ..$ 3            : num [1:1103] 0.000538 0.003235 0.007002 0.053386 0.006758 ...
    ##   ..$ 4            : num [1:1103] 0.01026 0.01518 0.06175 0.00725 0.0441 ...
    ##   ..$ 5            : num [1:1103] 0.00309 0.00311 0.06396 0.26165 0.03625 ...
    ##   ..$ 6            : num [1:1103] 0.000805 0.000199 0.007698 0.000485 0.003229 ...
    ##   ..$ 7            : num [1:1103] 0.000194 0.001645 0.019377 0.020895 0.011612 ...
    ##   ..$ 8            : num [1:1103] 0.000227 0.000321 0.035957 0.020804 0.012656 ...
    ##  $ colocProb    : num [1:8] 0.54511 0.19223 0.00446 0.04254 0.12468 ...
    ##  $ Tstat_L      : num [1:1103, 1:3] 2.808 -0.429 -1.014 1.569 -0.188 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:3] "moderate.t" "moderate.t" "moderate.t"
    ##  $ D0           : num [1:1103, 1:3] 0.0079 0.3636 0.2384 0.1165 0.3918 ...
    ##  $ D1           : num [1:1103, 1:3] 0.00163 0.01258 0.07835 0.02055 0.07563 ...

### False discovery rate assessment


<!-- README.md is generated from README.Rmd. Please edit that file -->
`iProFun`
=========

An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA copy number alterations (CNA) and DNA methylation
---------------------------------------------------------------------------------------------------------------------------------------------------

The goal of **iProFun** is to

-   characterize functional consequences of DNA copy number and methylation alterations in tumors

-   facilitate screening for cancer drivers contributing to tumor initiation and progression, since CNAs and DNA methylations that preserve functional consequences are more likely to be cancer drivers.

This package implement iProFun using 2 main functions. The primary function is (surprise!) `iProFun`, which first fit gene-level multiple linear regression to obtain summary statistics and then run Primo for detecting joint associations of DNA alternations with multi-omics traits. False discovery rate assessment is supported by `iProFun_permutate` function. A full description of the method can be found in our [paper](https://www.biorxiv.org/content/early/2018/12/06/488833).

### Installation

You can install the latest version directly from GitHub with [devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/Rpackage/iProFun")
```

<!-- * the most recent officially-released version from CRAN with -->
<!--     ```R -->
<!--     install.packages("iProFun") -->
<!--     ```` -->
<!-- * the latest development version from GitHub with -->
<!--     ```R -->
<!--     install.packages("devtools") -->
<!--     devtools::install_github("xiaoyu/Rpackage/iProFun") -->
<!--     ```` -->
Example of use
--------------

Below is an example of iProFun Integrative analysis pipeline.

### Data

The sample data: `cna`, `methy`, `rna`, `protein`, `phospho`, `rna_pc_1_3`, `protein_pc_1_3`, `phospho_pc_1_3` are included in the package. A full description of the sample dataset can be found in the help page of the data.

``` r
library(iProFun)
data(cna)
?cna
```

### iProFun Integrative analysis pipeline

To implement iProFun, we can use `iProFun`

``` r
iprofun_result <- iProFun(ylist = list(rna, protein, phospho), xlist = list(cna, methy),
  covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
  pi = rep(0.05, 3))
```

The result of the function contains a list with the following components. The first component shows the names for each pattern. `000` means none, `010` means only the second element in the xlist present.

``` r
iprofun_result$group
```

    ## [1] "000" "100" "010" "001" "110" "101" "011" "111"

The following two components shows the marginal posterior probability for each pattern.

``` r
iprofun_result$x_1_iProFun_marginal_probability
```

    ## [1] 1.661135e-02 8.670947e-02 2.500583e-15 1.313237e-09 3.321950e-01
    ## [6] 1.617580e-04 2.325013e-20 5.643224e-01

We can also show the posterior probability for each pattern.

``` r
head(iprofun_result$x_1_iProFun_gene_posterior_probability)
```

    ##       V1                   V2                   V3                   V4
    ## 1   AAAS 3.98651972763526e-47   0.0135167064625158 1.19589573270473e-58
    ## 2  ABCC1 2.86930837780733e-24 0.000145235497054389 1.28165875484421e-33
    ## 3   ABI1 1.65018230521891e-41  0.00582999259722584 2.16520084197717e-52
    ## 4   ABI2  1.7328182967301e-64 4.49105874050004e-06 7.40324635278225e-73
    ## 5 ABLIM1 4.36745372679193e-31 3.04920204275533e-07 8.89871350325967e-40
    ## 6   ACLY 1.65097499266829e-26  0.00737357961529641 8.16791423742635e-38
    ##                     V5                  V6                   V7
    ## 1 6.46603320228941e-54   0.329178926717965 3.82042925110791e-05
    ## 2 2.09422217933432e-31   0.526659231730965 1.84720701895405e-07
    ## 3 8.05497699145462e-49   0.621006303511799 4.95904041815281e-06
    ## 4  7.6289954142909e-71    0.15576850636726 3.44556407497044e-08
    ## 5  6.9987433104478e-36 0.00504367226101182 8.51481509625661e-08
    ## 6 3.15395853809779e-33   0.296149355355857 2.45465875285334e-05
    ##                     V8                V9
    ## 1  2.1830678128178e-63 0.657266162527008
    ## 2 1.05280474672465e-38 0.473195348051279
    ## 3 1.18948909020678e-57 0.373158744850557
    ## 4 3.66831336139173e-77 0.844226968118359
    ## 5 1.60490331945608e-42 0.994955937670633
    ## 6 1.75612955219011e-42 0.696452518441318

The last several elements in the list show the beta coefficients for each gene-level multiple linear regression.

``` r
head(iprofun_result$x_1_iProFun_gene_beta)
```

    ##       V1        1         2          3
    ## 1   AAAS 2.152949 0.7862445  0.7646263
    ## 2  ABCC1 1.314813 1.0466323  0.4915840
    ## 3   ABI1 1.709761 0.8168414 -0.2947058
    ## 4   ABI2 2.303057 1.4524548  1.0769021
    ## 5 ABLIM1 1.180402 1.0120973  1.5958646
    ## 6   ACLY 1.143833 0.5185449  0.7304007

To control false discovery rate for the result, we can use `iProFun_permutate`. To save time, I will set number of permutation to 1. In our paper, we permutate 100 times to control the empirical FDR.

``` r
iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho), xlist = list(cna, methy),
  covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
  pi = rep(0.05, 3), permutate_number = 1, thresholds =c(0.05,0.99,0.01), fdr = 0.1, filter = c(1,0))
```

    ## [1] "Finish Seed 1"

The result of the function contains a list with the following components. The first component shows the cutoff values for each group based on permutation.

``` r
iprofun_permutate_result$cutoff
```

    ##     cutoff_names cutoff
    ## 1 x_1_y_1_cutoff   0.99
    ## 2 x_1_y_2_cutoff   0.05
    ## 3 x_2_y_1_cutoff   0.99
    ## 4 x_2_y_2_cutoff   0.05
    ## 5 x_3_y_1_cutoff   0.99
    ## 6 x_3_y_2_cutoff   0.05

The following components show the genes names that are identified by iProFun with the FDR crteria we set in each group.

``` r
iprofun_permutate_result$x_2_y_2
```

    ## [1] "AP3D1"  "DNAJB6"

### Contributions

If you find small bugs, larger issues, or have suggestions, please file them using the [issue tracker](https://github.com/songxiaoyu/iProFun/issues) or email the maintainer at <Jiayi.Ji@mountsinai.org>. Contributions (via pull requests or otherwise) are welcome.

<!-- ## iProFun Integrative analysis pipeline -->
<!-- Below is an example of how iProFun is commonly used.  A full description of the tool can be found in our MCP paper. -->
<!-- ```{r, include = FALSE} -->
<!-- require(metRology) -->
<!-- require(matrixStats) -->
<!-- ``` -->
<!-- ```{r, messages = FALSE, warning = FALSE,} -->
<!-- library(iProFun) -->
<!-- ``` -->
<!-- ### Data summary -->
<!-- After preprossing and data cleaning, we have 15121 genes and 569 subjects for mRNA data, 7010 genes and 174 subjects for protein, 5685 genes and 70 subjects for phospho data, 25762 genes and 552 subjects for methylation data, 11859 genes and 560 subjects for mRNA data. The following shows the data structure for each data. -->
<!-- ```{r} -->
<!-- rna_normal[1:5,1:5] -->
<!-- ``` -->
<!-- ```{r, include = FALSE} -->
<!-- # protein_normal[1:5,1:5] -->
<!-- # phospho_normal[1:5, 1:5] -->
<!-- # methy[1:5,1:5] -->
<!-- # cnv[1:5, 1:5] -->
<!-- ``` -->
<!-- ### Gene-level multiple linear regression to obtain summary statistics -->
<!-- We use sets of separate regressions in the integrative analysis pipeline to allow for different samples being measured on different sets of molecular features. -->
<!-- ```{r, eval=FALSE} -->
<!-- ylist_normal = list(rna_normal, protein_normal, phospho_normal) -->
<!-- methy_input_1_3 <- -->
<!-- MultiReg_together( -->
<!-- ylist = ylist_normal, -->
<!-- xlist = list(methy, cnv), -->
<!-- covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), -->
<!-- cl = cl -->
<!-- ) -->
<!-- ``` -->
<!-- The following shows the results for CNA. -->
<!-- ```{r} -->
<!-- str(cnv_input_1_3) -->
<!-- ``` -->
<!-- ### Primo – An integrative analysis method for detecting joint associations of DNA al- terations with multi-omics traits -->
<!-- With the summary association statistics obtained from equations (1), we apply an integrative analysis method – Primo – to detect joint associations of DNA variation with multi-omics traits -->
<!-- ```{r, results=FALSE} -->
<!-- pi1=rep(0.05, 3) -->
<!-- cnv_1_3 = MultiOmics_Input(cnv_input_1_3, pi1=pi1) -->
<!-- cnv_1_3_tidy <- MultiOmics_Input(cnv_result , pi1=pi1) -->
<!-- cnv_1_3$colocProb *100 -->
<!-- cnv_1_3_tidy$colocProb*100 -->
<!-- ``` -->
<!-- ```{r} -->
<!-- str(cnv_1_3) -->
<!-- ``` -->
<!-- ### False discovery rate assessment -->
<!-- To calculate the empirical FDR, we first calculated the posterior probability of a predictor being associated with an outcome, by summing over all patterns that are consistent with the association of interest. -->
<!-- The following shows the results when we randomly permute the sample label of the mRNA while keeping the labels of the other two traits. -->
<!-- ```{r, eval=FALSE} -->
<!--  MultiReg_cnv_lr_perm_1 = MultiReg_together_perm( -->
<!--     ylist = list(rna_regression, protein_regression, phospho_regression), -->
<!--     xlist = list(cnv_lr_regression, cnv_baf_regression, methy_mean_regression), -->
<!--     covariates = list(purity_tumor,age, gender), -->
<!--     xyCommonGeneID = xy_common_geneID, -->
<!--     conditional_covariate = mutation_reg_111, -->
<!--     mutation_genes = mutation_gene_111, -->
<!--     xyCommonSubID = list(xrnaCommonSubID, xproteinCommonSubID, xphosphoCommonSubID), -->
<!--     filename = "MultiReg_cnv_lr_together_perm_1", -->
<!--     permcolum = 1, -->
<!--     seed=(currind-1)*10+i -->
<!--   ) -->
<!-- ``` -->


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
devtools::install_github("songxiaoyu/iProFun")
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

The result of the function contains a list with three tibbles. The first tibble shows the marginal posterior probability for each pattern.

``` r
iprofun_result$`Marginal Probability`
```

    ## # A tibble: 8 x 3
    ##   group                  CNA      Methy
    ##   <chr>                <dbl>      <dbl>
    ## 1 None              3.35e- 2 0.610     
    ## 2 RNA only          2.00e- 7 0.0465    
    ## 3 Protein only      1.01e-13 0.00000214
    ## 4 Phosphosite only  2.53e- 6 0.00778   
    ## 5 RNA & Protein     5.16e- 1 0.294     
    ## 6 RNA & Phospho     5.69e-10 0.00796   
    ## 7 Protein & Phospho 5.84e-16 0.0000181 
    ## 8 All three         4.50e- 1 0.0346

The second tibble shows the posterior probability for each pattern per methylation site.

``` r
iprofun_result$`Gene Posterior Probability`
```

    ## # A tibble: 1,779 x 11
    ##    X     Gene_ID Hybridization    None `RNA only` `Protein only`
    ##    <chr> <chr>   <chr>           <dbl>      <dbl>          <dbl>
    ##  1 Meth… AAAS    cg00559473    8.36e-1   0.0254    0.00000356   
    ##  2 Meth… AAAS    cg23032316    9.69e-1   0.0157    0.000000729  
    ##  3 Meth… ABCC1   cg17199483    9.81e-1   0.0132    0.000000282  
    ##  4 Meth… ABI1    cg17918201    2.73e-4   0.0983    0.00000000187
    ##  5 Meth… ABI1    cg23087130    6.06e-1   0.00559   0.0000309    
    ##  6 Meth… ABI2    cg09845946    2.08e-1   0.119     0.000000905  
    ##  7 Meth… ABLIM1  cg05064181    9.98e-1   0.000806  0.000000292  
    ##  8 Meth… ACLY    cg24124398    7.96e-1   0.0260    0.00000386   
    ##  9 Meth… ACLY    cg25687894    9.86e-1   0.00437   0.00000154   
    ## 10 Meth… ACTL6A  cg03667091    8.51e-1   0.0481    0.00000137   
    ## # ... with 1,769 more rows, and 5 more variables: `Phosphosite
    ## #   only` <dbl>, `RNA & Protein` <dbl>, `RNA & Phospho` <dbl>, `Protein &
    ## #   Phospho` <dbl>, `All three` <dbl>

The third tibble shows the beta coefficients for each gene-level multiple linear regression.

``` r
iprofun_result$Beta
```

    ## # A tibble: 1,779 x 9
    ##    Gene_ID estimate_cnv_rna estimate_cnv_pr… estimate_cnv_ph… Hybridization
    ##    <chr>              <dbl>            <dbl>            <dbl> <chr>        
    ##  1 AAAS               2.15             2.01             0.765 <NA>         
    ##  2 ABCC1              1.31             1.21             0.492 <NA>         
    ##  3 ABI1               1.71             1.96            -0.295 <NA>         
    ##  4 ABI2               2.30             2.07             1.08  <NA>         
    ##  5 ABLIM1             1.18             1.46             1.60  <NA>         
    ##  6 ACLY               1.14             0.702            0.705 <NA>         
    ##  7 ACTL6A             1.61             1.57             0.742 <NA>         
    ##  8 ACTR2              1.32             1.36            -0.207 <NA>         
    ##  9 ADAM17             1.41             1.68            -0.784 <NA>         
    ## 10 ADD2               0.954            0.776           -0.767 <NA>         
    ## # ... with 1,769 more rows, and 4 more variables: X <chr>,
    ## #   estimate_methy_rna <dbl>, estimate_methy_protein <dbl>,
    ## #   estimate_methy_phospho <dbl>

To control false discovery rate for the result, we can use `iProFun_permutate`. To save time, I will set number of permutation to 1. In our paper, we permutate 100 times to control the empirical FDR.

``` r
iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho), xlist = list(cna, methy),
  covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
  pi = rep(0.05, 3), permutate_number = 1)
```

    ## [1] "Finish Seed 1"

The result of the function contains a list with three tibbles. The first tibble shows the cutoff values for each group based on permutation.

``` r
iprofun_permutate_result$`Posterior Probability Cutoff`
```

    ## # A tibble: 6 x 2
    ##   Group         `Cut Off probability`
    ##   <chr>                         <dbl>
    ## 1 CNV RNA                        0.75
    ## 2 CNV Protein                    0.75
    ## 3 CNV Phospho                    0.75
    ## 4 Methy RNA                      0.75
    ## 5 Methy Protein                  0.75
    ## 6 Methy Phospho                  0.75

The second tibble shows whether a gene is identified by iProFun or not with the FDR crteria we set.

``` r
iprofun_permutate_result$`iProFun Result`
```

    ## # A tibble: 676 x 7
    ##    Gene  `CNV RNA` `Methy RNA` `Methy Protein` `CNV Protein`
    ##    <chr>     <dbl>       <dbl>           <dbl>         <dbl>
    ##  1 AAAS          1           0               0             1
    ##  2 ABCC1         1           0               0             1
    ##  3 ABI2          1           0               0             1
    ##  4 ABLI…         1           0               0             1
    ##  5 ACLY          1           0               0             1
    ##  6 ACTL…         1           0               0             1
    ##  7 ADD3          1           0               0             1
    ##  8 ADNP          1           0               0             1
    ##  9 AHNAK         1           0               0             1
    ## 10 AKAP1         1           0               0             1
    ## # ... with 666 more rows, and 2 more variables: `Methy Phospho` <dbl>,
    ## #   `CNV Phospho` <dbl>

The third tibble shows add the direction of the assocation to the second tibble.

``` r
iprofun_permutate_result$`iProFun Result (Negative/Positive)`
```

    ## # A tibble: 676 x 7
    ##    Gene  `CNV RNA` `Methy RNA` `Methy Protein` `CNV Protein`
    ##    <chr>     <dbl>       <dbl>           <dbl>         <dbl>
    ##  1 AAAS          1           0               0             1
    ##  2 ABCC1         1           0               0             1
    ##  3 ABI2          1           0               0             1
    ##  4 ABLI…         1           0               0             1
    ##  5 ACLY          1           0               0             1
    ##  6 ACTL…         1           0               0             1
    ##  7 ADD3          1           0               0             1
    ##  8 ADNP          1           0               0             1
    ##  9 AHNAK         1           0               0             1
    ## 10 AKAP1         1           0               0             1
    ## # ... with 666 more rows, and 2 more variables: `Methy Phospho` <dbl>,
    ## #   `CNV Phospho` <dbl>

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

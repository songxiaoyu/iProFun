
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


<!-- README.md is generated from README.Rmd. Please edit that file -->
`iProFun`
=========

An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA copy number alterations (CNA) and DNA methylation
---------------------------------------------------------------------------------------------------------------------------------------------------

The goal of **iProFun** is to

-   characterize functional consequences of DNA copy number and methylation alterations in tumors

-   facilitate screening for cancer drivers contributing to tumor initiation and progression, since CNAs and DNA methylations that preserve functional consequences are more likely to be cancer drivers.

This package implement iProFun using 2 main functions. The primary function is (surprise!) `iProFun`, which first fit multiple linear regression to obtain summary statistics and then use the statistics to detect joint associations of DNA alternations with multi-omic traits. False discovery rate assessment is supported by `iProFun_permutate` function. A full description of the method can be found in our [paper](https://www.biorxiv.org/content/early/2018/12/06/488833).

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
iprofun_result <- iProFun(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3,
protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3))
```

    ## [1] "No IDs are specified. The first column of each data type in ylist and xlist are used."
    ## [1] "A total of 676 genes/proteins/peptides/etc. are considered."
    ## [1] "A total of 525 samples for Y1 data type."
    ## [2] "A total of 171 samples for Y2 data type."
    ## [3] "A total of 68 samples for Y3 data type." 
    ## [1] "Obtaining regression summaries for data type 1"
    ## [1] "Obtaining regression summaries for data type 2"
    ## [1] "Obtaining regression summaries for data type 3"
    ## [1] "Run iProFun"
    ## [1] "Completed"

The result of the function contains a list. The length of the output list equals to the number of element in the xlist. In each element of the list, there is a sublist containing the results associated with the data in xlist.

The following componentshows the marginal posterior probability for each pattern associated with the first element in xlist.

``` r
iprofun_result$`iProFun output for xlist 1`$colocProb
```

    ## [1] 1.620700e-02 7.908244e-02 1.803997e-10 3.075494e-07 4.505220e-01
    ## [6] 1.341961e-03 3.813558e-13 4.528463e-01

We can also show the posterior probability for each pattern associated with the first element in xlist.

``` r
head(cbind(iprofun_result$`iProFun output for xlist 1`$yName_J, iprofun_result$`iProFun output for xlist 1`$PostProb))
```

    ##             phospho_ID                                         
    ## 1  "AAAS"   "AAAS.NP_001166937.1:s462"   "3.9237564406578e-47" 
    ## 2  "ABCC1"  "ABCC1.NP_004987.2:s930"     "2.56137582110385e-24"
    ## 4  "ABI1"   "ABI1.NP_001012768.1:s222"   "1.40677992841573e-41"
    ## 5  "ABI2"   "ABI2.NP_001269854.1:s183"   "1.87012825266245e-64"
    ## 9  "ABLIM1" "ABLIM1.NP_001309811.1:s420" "4.94291355528985e-31"
    ## 16 "ACLY"   "ACLY.NP_001087.2:s481"      "1.76491191405199e-26"
    ##                                                                        
    ## 1  "0.0123749091873888"   "8.96383959025961e-54" "1.54463510952455e-51"
    ## 2  "0.000120595712748723" "8.71286778050195e-29" "4.49085550402044e-29"
    ## 4  "0.00462301449129315"  "1.40567545207946e-47" "1.64341316443298e-46"
    ## 5  "4.50847967166922e-06" "6.08461639892011e-68" "2.0230991904436e-68" 
    ## 9  "3.20999177143844e-07" "7.66963265890669e-35" "2.0376595190263e-33" 
    ## 16 "0.00733202484142715"  "6.64945949316372e-33" "7.31464611758034e-31"
    ##                                                                       
    ## 1  "0.449574541972872"   "0.000326191725840336" "3.64364112491387e-56"
    ## 2  "0.652360027007038"   "1.41577651564761e-06" "1.57737356460954e-31"
    ## 4  "0.734601494376669"   "3.61620391871919e-05" "1.69559964669517e-50"
    ## 5  "0.233270372819702"   "3.26575018606226e-07" "6.7966772468376e-70" 
    ## 9  "0.00792069194509253" "8.86052756345061e-07" "3.26468154769303e-35"
    ## 16 "0.439293617050499"   "0.000203470415596139" "2.84560285735056e-35"
    ##                       
    ## 1  "0.537724357113899"
    ## 2  "0.347517961503697"
    ## 4  "0.26073932909285" 
    ## 5  "0.766724792125608"
    ## 9  "0.992078101002974"
    ## 16 "0.553170887692478"

The following shows the beta coefficients for each gene-level multiple linear regression.

``` r
head(cbind(iprofun_result$`iProFun output for xlist 1`$yName_J, iprofun_result$`iProFun output for xlist 1`$betas_J))
```

    ##             phospho_ID                                     
    ## 1  "AAAS"   "AAAS.NP_001166937.1:s462"   "2.15294856564555"
    ## 2  "ABCC1"  "ABCC1.NP_004987.2:s930"     "1.31481344483779"
    ## 4  "ABI1"   "ABI1.NP_001012768.1:s222"   "1.70976139319942"
    ## 5  "ABI2"   "ABI2.NP_001269854.1:s183"   "2.30305705519527"
    ## 9  "ABLIM1" "ABLIM1.NP_001309811.1:s420" "1.18040186037167"
    ## 16 "ACLY"   "ACLY.NP_001087.2:s481"      "1.14383323720261"
    ##                                           
    ## 1  "0.786244511025172" "0.764626299444886"
    ## 2  "1.04663228802151"  "0.491584008481416"
    ## 4  "0.816841442720532" "-0.2947058215859" 
    ## 5  "1.45245476052363"  "1.07690206117779" 
    ## 9  "1.01209731356699"  "1.59586460265822" 
    ## 16 "0.518544940554173" "0.704512877526666"

To control false discovery rate for the result, we can use `iProFun_permutate`. To save time, I will set number of permutation to 1. In our paper, we permutate 100 times to control the empirical FDR.

``` r
iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho),
xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
pi = rep(0.05, 3), permutate_number = 1, fdr = 0.1, PostCut = 0.75, filter <- c(1,0),
grids = c(seq(0.75, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001)),
seed=123)
```

    ## [1] "Finish permutation 1"

    ## Warning in min(which(f2 < fdr)): no non-missing arguments to min; returning
    ## Inf

The result of the function contains a list with the following components. The first component shows the cutoff values for each group based on permutation.

``` r
iprofun_permutate_result$fdr_cutPob
```

    ## [[1]]
    ## [1] 0.75 0.79 0.75
    ## 
    ## [[2]]
    ## [1] 0.75 0.75   NA

The following components show the genes names that are identified by iProFun with the FDR crteria we set in each group.

``` r
head(iprofun_permutate_result$Gene_fdr[[1]])
```

    ##               Y1  Y2  Y3 
    ## [1,] "AAAS"   "1" "1" "0"
    ## [2,] "ABCC1"  "1" "1" "0"
    ## [3,] "ABI1"   "0" "0" "0"
    ## [4,] "ABI2"   "1" "1" "1"
    ## [5,] "ABLIM1" "1" "1" "1"
    ## [6,] "ACLY"   "1" "1" "0"

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

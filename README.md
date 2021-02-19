
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `iProFun`

## An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA copy number alterations (CNA) and DNA methylation

The goal of **iProFun** is to

  - characterize functional consequences of DNA copy number and
    methylation alterations in tumors

  - facilitate screening for cancer drivers contributing to tumor
    initiation and progression, since CNAs and DNA methylations that
    preserve functional consequences are more likely to be cancer
    drivers.

This package implement iProFun using 2 main functions. The primary
function is (surprise\!) `iProFun`, which first fit multiple linear
regression to obtain summary statistics and then use the statistics to
detect joint associations of DNA alternations with multi-omic traits.
False discovery rate assessment is supported by `iProFun_permutate`
function. A full description of the method can be found in our
[paper](https://www.mcponline.org/content/early/2019/06/21/mcp.RA118.001220).

### Installation

You can install the latest version directly from GitHub with
[devtools](https://github.com/hadley/devtools):

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

## Example of use

Below is an example of iProFun Integrative analysis pipeline.

### Data

The sample data: `cna`, `methy`, `rna`, `protein`, `phospho`,
`rna_pc_1_3`, `protein_pc_1_3`, `phospho_pc_1_3` are included in the
package. A full description of the sample dataset can be found in the
help page of the data.

``` r
library(iProFun)
data(cna)
?cna
```

### iProFun Integrative analysis pipeline

For analysis with overlapping genes, use:

``` r
yList = list(rna, protein, phospho); xList = list(cna, methy)
covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3)
pi1 = 0.05
```

Alternatively, we pretend the outcomes have some genes only avaiable in
some data types:

``` r
yList = list(rna[1:400,], protein, phospho[100:1590,])
```

And we pretend there are four predictor data types with a large or small
number of genes.

``` r
xList = list(cna[1:300,], cna[301:310,], methy[21:1103,],  methy[1:20,])
```

Regression on one outcome data
type:

``` r
ft1=iProFun.reg.1y(yList.1y=yList[[1]], xList=xList, covariates.1y=covariates[[1]],
                    var.ID=c("Gene_ID"),
                    var.ID.additional=c("phospho_ID", "Hybridization", "chr"))
```

Regression on all three outcome data types

``` r
reg.all=iProFun.reg(yList=yList, xList=xList, covariates=covariates,
                    var.ID=c("Gene_ID"), var.ID.additional=c("Gene_ID", "phospho_ID",
                    "Hybridization", "chr"))
```

Reformat the regression summaries for iProFun

``` r
summ=multi.omic.reg.summary(reg.out.list=reg.all, var.ID="Gene_ID")
```

FWER controlled identification for the data type with few genes

``` r
FWER=iProFun.FWER(Reg.Sum=summ, FWER.Index=c(2,4), filter=c(0,0))
```

Calculate the posterior probabilities of association patterns via
iProFun

``` r
prob=iProFun.prob(Reg.Sum=summ, NoProbXIndex=c(2,4), pi1=pi1)
```

Summarize posterior probabilities for one outcome of interest

``` r
prob1y=iProFun.sum.prob.1y(prob=prob, which.y=1, NoProbXIndex=c(2,4))
```

Fast FDR calculation (may not be
accurate)

``` r
fastFDR=estFDR(ProbPattern1y1x=prob1y[[1]], grids = seq(0.01, 0.99, by=0.01))
```

Empirical FDR controlled discoveries for one
outcome

``` r
eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=1, yList=yList, xList=xList,
                      covariates=covariates, pi1=pi1,
                      NoProbXIndex=c(2,4), filter=c(1, -1),
                      permutate_number=2, var.ID=c("Gene_ID"),
                      grids = seq(0.01, 0.99, by=0.01),fdr = 0.1, PostCut=0.75,
                      var.ID.additional=c( "phospho_ID", "Hybridization", "chr"))
```

Empirical FDR controlled discoveries for all
outcomes

``` r
eFDR=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
                 NoProbXIndex=c(2,4),filter=c(1, -1),
                 permutate_number=2, var.ID=c("Gene_ID"),
                 grids = seq(0.01, 0.99, by=0.01),fdr = 0.1, PostCut=0.75,
                  var.ID.additional=c( "phospho_ID", "Hybridization", "chr"), seed=NULL)
```

    ## [1] "Outcome" "1"      
    ## [1] "Outcome" "2"      
    ## [1] "Outcome" "3"

### Contributions

If you find small bugs, larger issues, or have suggestions, please file
them using the [issue
tracker](https://github.com/songxiaoyu/iProFun/issues) or email the
maintainer at <Jiayi.Ji@mountsinai.org>. Contributions (via pull
requests or otherwise) are
welcome.

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

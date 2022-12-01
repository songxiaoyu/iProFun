
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `iProFun`

## An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA-level alterations (e.g. somatic mutation, copy number variation (CNV) and DNA methylation)

The goal of **iProFun** is to characterize multi-omics functional
consequences of DNA-level alterations in tumor.

-   iProFun starts with linear regressions that consider multiple -omic
    outcomes (e.g. mRNA, protein and phosphoprotein) separately, so it
    allows different genes and samples for different outcomes. This
    analysis can be performed by `iProFun.reg` function.

-   iProFun uses the summary statistics from multiple regressions to
    jointly detect their associations with DNA-level alternations.

For data types with few events (e.g. somatic mutations), `iProFun`
provides estimate, standard error, Student’s t-test p-value, family-wise
error rate (FWER), multi-omic directional filtering, and whether it’s
identified by iProFun.

-   The `iProFun` identification here is determined by FWER\< a cutoff &
    pass of a directional filtering criterion.

For data types with many genes, where parallel features of the genes can
be learned from each other to boost study power, `iProFun` provides
estimate, standard error, Student’s t-test p-value, posterior
association probability, empirical false discovery rate (eFDR),
multi-omic directional filtering, and whether it’s identified by
iProFun.

-   The `iProFun` identification here is determined by

    1.  posterior association probability > a cutoff (e.g. 75%),
    2.  eFDR \< a cutoff (e.g. 0.1), and
    3.  pass of a directional filtering criterion (e.g. no filter, all
        positive/negative associations across all outcome data types,
        consistent association directions across all outcome data
        types).

A full description of the method can be found in our
[paper](https://pubmed.ncbi.nlm.nih.gov/31227599/).

### Installation

You can install the latest version directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/iProFun")
```

## Example of use

Below is an example of iProFun Integrative analysis pipeline.

### Sample data

The preprocessed data from National Cancer Institute’s Clinical
Proteomic Tumor Analysis Consortium Lung Squamous Cell Carcinoma (lscc)
study are included in the package, including `cnv`, `mut`, `rna`,
`protein`, `phospho` and `cov`.

-   Data reference: ’Satpathy, Shankha, et al. “A proteogenomic portrait
    of lung squamous cell carcinoma.” Cell 184.16 (2021): 4348-4371.

A brief description of the sample dataset can be found in the help page
of the data.

``` r
library(iProFun)
data(lscc_iProFun_Data) # load all data
objects() # list all loaded data
```

    ## [1] "cnv"     "cov"     "mut"     "phospho" "protein" "rna"

``` r
?cnv # help file is available for each individual data
```

### iProFun integrative analysis pipeline

#### - Example 1: rna/protein/phospho \~ mutation +cnv + covariates

Analysis should specify multi-omics outcome data types, predictor data
types, covariates and prior association probability. Here we consider
RNA, protein and phosphoprotein as outcomes, and mutation and cnv as
predictors, and use the same set of covariates for three outcomes. We
use a conservative prior `pi1 = 0.05`. Note, the impact of prior is
minimal on the results.

``` r
yList = list(rna, protein, phospho); xList = list(mut, cnv)
covariates = list(cov, cov, cov) # iProFun allows different covariates for different regressions, and here we repeat the same covariates for simplicity
pi1 = 0.05 # prior association probability. 
```

Try regression on one outcome data type for checking the implementation

``` r
ft1=iProFun.reg.1y(yList.1y=yList[[1]], xList=xList, covariates.1y=covariates[[1]], 
                   var.ID=c("geneSymbol"))
```

The result `ft1` is a list, which contains

-   xName (Predictor variable name corresponds to each predictor-outcome
    pair),
-   yName (Outcome variable name corresponds to each predictor-outcome
    pair),
-   betas (Coefficient estimate for predictors),
-   betas_se (Coefficent SE for predictors),
-   sigma2 (Regression error terms for predictors),
-   dfs (Regression degrees of freedom for predictors),
-   v_g ((X^T X)^-1 projection on predictors). Each of the element
    contains a sublist, which stores the corresponding results for each
    of the element in xList.

For multi-omic iProFun analysis, we need regression on all three outcome
data types:

``` r
reg.all=iProFun.reg(yList=yList, xList=xList, covariates=covariates, 
                    var.ID=c("geneSymbol"), var.ID.additional=c("id"))
```

The result `reg.all` is a list with length equals to the length of
`yList`. The first element `reg.all[[1]]` is essentially the same as
`ft1` since both of them store the results between `yList[[1]]` and
`xList`. Similarly the second element `reg.all[[2]]` stores the results
between `yList[[2]]` and `xList` and so on and so forth.

If one is interested to save regression results in a single table, one
can use this function to output the result table

``` r
reg.tab=iProFun.reg.table(reg.all=reg.all, xType = c("mutation", "cnv"), 
                          yType = c("rna", "protein", "phospho"))
```

This function calculates FWER. It’s preferred to be used for the data
type with few genes, such as somatic mutation. Mutation is the first
element in the xList, so we use `FWER.Index=c(1)` to calculate FWER for
this element.

``` r
FWER.all=iProFun.FWER(reg.all=reg.all, FWER.Index=c(1))
```

This function calculates posterior association probability and eFDR rate
for the predictor data types on one outcome. It’s preferred to be used
for data types with many genes, such as cnv. As, we don’t want to
calculate the probabilities of association patterns between the mutation
(1st element of xList) and yList and we set `NoProbXIndex = c(1)`. For a
fast demonstration, we permute the data only twice using
`permutate_number=2`.

``` r
eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=1, yList=yList, xList=xList,
                      covariates=covariates, pi1=pi1, NoProbXIndex=c(1),
                      permutate_number=2, var.ID=c("geneSymbol"), 
                      var.ID.additional=c("id"))
```

    ## [1] "perm" "1"   
    ## [1] "perm" "2"

This is an expansion of `iProFun.eFDR.1y` to calculate eFDR for all
outcomes.

``` r
eFDR.all=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
                  NoProbXIndex=c(1), permutate_number=2, var.ID=c("geneSymbol"), 
                  var.ID.additional=c("id"), seed=123)
```

    ## [1] "Outcome" "1"      
    ## [1] "perm" "1"   
    ## [1] "perm" "2"   
    ## [1] "Outcome" "2"      
    ## [1] "perm" "1"   
    ## [1] "perm" "2"   
    ## [1] "Outcome" "3"      
    ## [1] "perm" "1"   
    ## [1] "perm" "2"

The result `eFDR` is a list with length equals to the length of `yList`.
The first element `eFDR[[1]]` is essentially the same as `eFDR1` since
both of them store the results between `eFDR[[1]]` and `xList`.

This function provides `iProFun` identifications for all predictors and
outcomes, based on FWER/eFDR, association probabilities, and biological
directional filtering. The output has been reformatted to a long-format
table for usage.

``` r
# iProFun identification
# For data types with many genes, it's based on 
# (1) association probabilities > 0.75 as specified by `PostPob.cutoff=0.75`,
# (2) FDR 0.1 as specified by `fdr.cutoff = 0.1`, and 
# (3) the association direction filtering (CNV requires positive associations as specified by the second element of`filter=c(0, 1)` ).
# For data types with few genes, it's  based on 
# (1) FWER 0.1  as specified by `fwer.cutoff=0.1`, and
# (2)  the association direction filtering (mutation requires consistent association directions as specified by the first element of`filter=c(0, 1)`).
res=iProFun.detection(reg.all=reg.all, eFDR.all=eFDR.all, FWER.all=FWER.all, filter=c(0, 1),
                      NoProbButFWERIndex=1,fdr.cutoff = 0.1, fwer.cutoff=0.1, PostPob.cutoff=0.75,
                      xType=c("mutation", "cnv"), yType=c("rna", "protein", "phospho"))
```

Output some results

``` r
head(res)
```

    ##    xName yName.1                           yName.2    xType yType          est
    ## 1   TP53    TP53    NP_000537.3_S315s _1_0_314_315 mutation   rna -0.202703364
    ## 2   PTEN    PTEN NP_001291646.2_S467s _1_1_467_467 mutation   rna -0.059590144
    ## 3 CDKN2A    <NA>                              <NA> mutation   rna  1.231840027
    ## 4  KMT2D   KMT2D NP_003473.3_T1843t _1_1_1843_1843 mutation   rna -0.185339293
    ## 5 NFE2L2  NFE2L2    NP_006155.2_S215s _1_1_215_215 mutation   rna  0.211863472
    ## 6 ARID1A  ARID1A NP_006006.3_S1755s _1_1_1755_1755 mutation   rna -0.004243765
    ##          se      pvalue      FWER eFDR PostProb d.filter iProFun.identification
    ## 1 0.4010039 0.614294914 1.0000000   NA       NA        0                      0
    ## 2 0.1287821 0.644540765 1.0000000   NA       NA        0                      0
    ## 3 0.4643312 0.009243121 0.1201606   NA       NA        0                      0
    ## 4 0.1196013 0.124294431 1.0000000   NA       NA        1                      0
    ## 5 0.2408745 0.381145481 1.0000000   NA       NA        1                      0
    ## 6 0.1344626 0.974883272 1.0000000   NA       NA        0                      0

#### - Example 2: rna/protein/phospho \~ cnv

This time, we try rna/protein/phospho \~ cnv to see how iProFun works
when we need to calcualte association probabilities for all predictors.

``` r
# We still need to put cnv into a list 
yList = list(rna, protein, phospho); xList = list(cnv)
pi1 = 0.05 # prior association probability. 
```

Again, we start with regression on all three outcome data types:

``` r
reg.all=iProFun.reg(yList=yList, xList=xList, covariates=NULL, 
                    var.ID=c("geneSymbol"), var.ID.additional=c("id"))
```

To save regression results in a single table, one can use this function

``` r
reg.tab=iProFun.reg.table(reg.all=reg.all, xType = c("cnv"), 
                          yType = c("rna", "protein", "phospho"))
```

We skip the function to calculate FWER for predictors like mutation that
exists in few genes, and directly calculate posterior association
probabilities. In this case, we should specify `NoProbXIndex=NULL`.

``` r
eFDR.all=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=NULL,pi1=pi1,
                  NoProbXIndex=NULL, permutate_number=2, var.ID=c("geneSymbol"), 
                  var.ID.additional=c("id"), seed=123)
```

    ## [1] "Outcome" "1"      
    ## [1] "perm" "1"   
    ## [1] "perm" "2"   
    ## [1] "Outcome" "2"      
    ## [1] "perm" "1"   
    ## [1] "perm" "2"   
    ## [1] "Outcome" "3"      
    ## [1] "perm" "1"   
    ## [1] "perm" "2"

To summarize the results in a long-format table, we use
`iProFun.detection`.

``` r
# iProFun identification is based on 
# (1) association probabilities > 0.75 as specified by `PostPob.cutoff=0.75`,
# (2) FDR 0.1 as specified by `fdr.cutoff = 0.1`, and 
# (3) the association direction filtering (CNV requires positive associations as specified by `filter=c(1)` ).

res=iProFun.detection(reg.all=reg.all, eFDR.all=eFDR.all, FWER.all=NULL, filter=c( 1),NoProbButFWERIndex=NULL,fdr.cutoff = 0.1, fwer.cutoff=NULL, PostPob.cutoff=0.75,
                      xType=c("cnv"), yType=c("rna", "protein", "phospho"))
```

Output some results

``` r
head(res)
```

    ##   xName yName.1 yName.2 yName.3                           yName.4 xType yType
    ## 1  A1CF    A1CF    <NA>    <NA>                              <NA>   cnv   rna
    ## 2  ABI1    ABI1    ABI1    ABI1 NP_001171590.1_S183s _1_1_183_183   cnv   rna
    ## 3  ABL1    ABL1    ABL1    ABL1    NP_009297.2_S828s _1_0_823_828   cnv   rna
    ## 4  ABL2    ABL2    ABL2    ABL2    NP_009298.1_S631s _1_1_631_631   cnv   rna
    ## 5 ACKR3   ACKR3   ACKR3   ACKR3    NP_064707.1_S350s _1_1_350_350   cnv   rna
    ## 6 ACSL3   ACSL3   ACSL3   ACSL3 NP_001341087.1_S683s _1_1_683_683   cnv   rna
    ##          est        se       pvalue FWER        eFDR  PostProb d.filter
    ## 1 -1.5186404 1.1998222 2.133184e-01   NA 0.155555556 0.7863618        0
    ## 2  0.5554881 0.1080023 1.244631e-06   NA 0.006313131 1.0000000        1
    ## 3  0.9888841 0.1431169 3.756320e-10   NA 0.006313131 1.0000000        1
    ## 4  0.5337644 0.1334711 1.178884e-04   NA 0.006313131 1.0000000        1
    ## 5  1.7055492 0.6351614 8.414236e-03   NA 0.119607843 0.9066117        1
    ## 6  0.8553145 0.1804231 6.663326e-06   NA 0.006313131 1.0000000        1
    ##   iProFun.identification
    ## 1                      0
    ## 2                      1
    ## 3                      1
    ## 4                      1
    ## 5                      0
    ## 6                      1

### Contributions

If you find small bugs, larger issues, or have suggestions, please file
them using the [issue
tracker](https://github.com/songxiaoyu/iProFun/issues) or email the
maintainer at <xiaoyu.song@mountsinai.org>. Contributions (via pull
requests or otherwise) are welcome.

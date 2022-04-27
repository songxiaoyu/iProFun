\docType{package}
\name{iProFun-package}
\alias{iProFun}
\alias{iProFun-package}
\title{iProFun: Characterization of multi-omic functional consequences of DNA-level alterations in tumor}
\description{
The goal of iProFun is to characterize multi-omic functional consequences of DNA-level alterations  in tumor.

For data types with few genes, `iProFun` provides estimate, standard error, Student's t-test p-value, family-wise error rate (FWER), multi-omic directional filtering, and whether it's identified by iProFun or not.
 * `iProFun` identification is determined by FWER< a cutoff & pass of a directional filtering criterion.

For data types with many genes, where parallel features of the genes can be learned from each other to boost study power, `iProFun` provides estimate, standard error, Student's t-test p-value, posterior association probability, empirical false discovery rate (eFDR), multi-omic directional filtering, and whether it's identified by iProFun or not.
 * `iProFun`  identification is determined by  posterior association probability > a cutoff, eFDR < a cutoff & pass of a directional filtering criterion.
}
\examples{
# Load data
data(lscc_iProFun_Data)
# For analysis with overlapping genes, use:
yList = list(rna, protein, phospho); xList = list(mut, cnv)
covariates = list(cov, cov, cov)
pi1 = 0.05
# Regression on one outcome data type for checking the data
ft1=iProFun.reg.1y(yList.1y=yList[[1]], xList=xList, covariates.1y=covariates[[1]],
                   var.ID=c("geneSymbol"))
# Regression on all three outcome data types
reg.all=iProFun.reg(yList=yList, xList=xList, covariates=covariates,
                    var.ID=c("geneSymbol"), var.ID.additional=c("id"))
# Calculate FWER for data type(s) that have few number of genes
FWER.all=iProFun.FWER(reg.all=reg.all, FWER.Index=c(1))
# Calculate Empirical FDR for one outcome for checking the implementation in a data
eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=2, yList=yList, xList=xList,
                      covariates=covariates, pi1=pi1, NoProbXIndex=c(1),
                      permutate_number=2, var.ID=c("geneSymbol"),
                      var.ID.additional=c("id"))
# Calculate Empirical FDR for all outcomes
eFDR.all=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
                  NoProbXIndex=c(1), permutate_number=2, var.ID=c("geneSymbol"),
                  var.ID.additional=c("id"), seed=123)
# iProFun identification
# For data types with many genes, it's based on (1) association probabilities > 0.75,
# (2) FDR 0.1, and (3) the association direction filtering (CNV requires positive associations).
# For data types with few genes, it's  based (1) FWER 0.1, and
# (2)  the association direction filtering (mutation requires consistent association directions).
res=iProFun.detection(reg.all=reg.all, eFDR.all=eFDR.all, FWER.all=FWER.all, filter=c(0, 1),
                      NoProbButFWERIndex=1,fdr.cutoff = 0.1, fwer.cutoff=0.1, PostPob.cutoff=0.75,
                      xType=c("mutation", "cnv"), yType=c("rna", "protein", "phospho"))
# output examples
head(res)
}


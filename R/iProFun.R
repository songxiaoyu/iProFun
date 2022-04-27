#' @description
#' @useDynLib iProFun
#' @keywords internal
#' @import rlang
#' "_PACKAGE"
#' @examples
#' # Load data
#' data(lscc_iProFun_Data)
#' # For analysis with overlapping genes, use:
#' yList = list(rna, protein, phospho); xList = list(mut, cnv)
#' covariates = list(cov, cov, cov)
#' pi1 = 0.05
#' # Regression on one outcome data type
#' ft1=iProFun.reg.1y(yList.1y=yList[[1]], xList=xList, covariates.1y=covariates[[1]],
#'                    var.ID=c("geneSymbol"))
#' # Regression on all three outcome data types
#' reg.all=iProFun.reg(yList=yList, xList=xList, covariates=covariates,
#'                     var.ID=c("geneSymbol"), var.ID.additional=c("id"))
#' # Calculate FWER for data type(s) that have few number of genes
#' FWER=iProFun.FWER(reg.all=reg.all, FWER.Index=c(1), filter=c(0))
#' # Calculate Empirical FDR for one outcome
#' eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=2, yList=yList, xList=xList,
#'                       covariates=covariates, pi1=pi1, NoProbXIndex=c(1),
#'                       permutate_number=2, var.ID=c("geneSymbol"),
#'                       var.ID.additional=c("id"))
#' # Calculate Empirical FDR for all outcomes
#'eFDR=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
#'                  NoProbXIndex=c(2,4),
#'                  permutate_number=2, var.ID=c("Gene_ID"),
#'                   var.ID.additional=c( "id"), seed=123)
#' # iProFun identification
#' # For data types with abundance genes, it's  based (1) association probabilities > 0.75,
#' # (2) FDR 0.1, or (3) the association direction filtering.
#' # For data types with few genes, it's  based (1) FWER 0.1
#' # (2)  the association direction filtering.
#'
#'res=iProFun_detection(reg.all=reg.all, eFDR.all=eFDR.all, FWER.all=FWER.all, filter=c(0, 1),
#'NoProbButFWERIndex=1,fdr.cutoff = 0.1, fwer.cutoff=0.1, PostPob.cutoff=0.75,
#'xType=c("mutation", "cnv"), yType=c("rna", "protein", "phospho"))


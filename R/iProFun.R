#' @description
#' @useDynLib iProFun
#' @keywords internal
#' @import rlang
#' "_PACKAGE"
#' @examples
#' # Load data
#' data(cna, methy, rna, protein, phospho, rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3)
#' # For analysis with overlapping genes, use:
#' yList = list(rna, protein, phospho); xList = list(cna, methy)
#' covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3)
#' pi1 = 0.05
#' # Alternatively, we pretend the outcomes have some genes only avaiable in some data types:
#' yList = list(rna[1:400,], protein, phospho[100:1590,])
#' # And we pretend there are four predictor data types with a large or small number of genes.
#' xList = list(cna[1:300,], cna[301:310,], methy[21:1103,],  methy[1:20,])
#' # Regression on one outcome data type
#' ft1=iProFun.reg.1y(yList.1y=yList[[1]], xList=xList, covariates.1y=covariates[[1]],
#'                    var.ID=c("Gene_ID"),
#'                    var.ID.additional=c("phospho_ID", "Hybridization", "chr"))
#' # Regression on all three outcome data types
#' reg.all=iProFun.reg(yList=yList, xList=xList, covariates=covariates,
#'                     var.ID=c("Gene_ID"), var.ID.additional=c("Gene_ID", "phospho_ID",
#'                     "Hybridization", "chr"))
#' # Reformat the regression summaries for iProFun
#' summ=multi.omic.reg.summary(reg.out.list=reg.all, var.ID="Gene_ID")
#' # FWER controlled identification for the data type with few genes
#' FWER=iProFun.FWER(Reg.Sum=summ, FWER.Index=c(2,4), filter=c(0,0))
#' # Calculate the posterior probabilities of association patterns via iProFun
#' prob=iProFun.prob(Reg.Sum=summ, NoProbXIndex=c(2,4), pi1=pi1)
#' # Summarize posterior probabilities for one outcome of interest
#' prob1y=iProFun.sum.prob.1y(prob=prob, which.y=1, NoProbXIndex=c(2,4))
#' # Fast FDR calculation (may not be accurate)
#' fastFDR=estFDR(ProbPattern1y1x=prob1y[[1]], grids = seq(0.01, 0.99, by=0.01))
#' # Empirical FDR controlled discoveries for one outcome
#' eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=2, yList=yList, xList=xList,
#'                       covariates=covariates, pi1e=pi1,
#'                       NoProbXIndex=c(2,4), filter=c(1, -1),
#'                       permutate_number=2, var.ID=c("Gene_ID"),
#'                       grids = seq(0.01, 0.99, by=0.01),fdr = 0.1, PostCut=0.75,
#'                       var.ID.additional=c( "phospho_ID", "Hybridization", "chr"))
#' # Empirical FDR controlled discoveries for all outcomes
#'eFDR=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1e=pi1,
#'                  NoProbXIndex=c(2,4),filter=c(1, -1),
#'                  permutate_number=2, var.ID=c("Gene_ID"),
#'                  grids = seq(0.01, 0.99, by=0.01),fdr = 0.1, PostCut=0.75,
#'                   var.ID.additional=c( "phospho_ID", "Hybridization", "chr"), seed=NULL)
#'
#'

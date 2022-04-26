
#' Linear regression on one outcome data type
#'
#'  Linear regression on one outcome data type with all data types of DNA alterations.
#' @export iProFun.reg.1y
#' @param yList.1y yList is a list of data matrix for outcomes, and yList.1y is one element of the
#' list indicating the outcome on one data type.
#' @param xList xList is a list of data matrix for predictors.
#' @param covariates.1y covariates is a list of data matrix for covariates, and covariates.1y is one
#' element of the list indicating the covariates for one data type. This list should be NULL or have the
#' same No. of subjects as ylist.1y.
#' @param permutation whether to permuate the label of the outcome. permutation = F (default):
#' no permuatation and it should be used for analysis of original data. permutation = T: permutate the
#' label of outcome, which is useful in generating eFDR controlled discoveries.
#' @param var.ID var.ID gives the variable name (e.g. gene/protein name) to match different data types.
#' @param var.ID.additional var.ID.additional allows to output additional variable names from the input.
#' Often helpful if multiple rows (e.g. probes) are considered per gene to allow clear index of the rows.
#' @param seed seed allows users to externally assign seed to replicate results. Useful when permutation=T.
#' @return It contains
#' \item{xName:}{Predictor variable name corresponds to each predictor-outcome pair}
#' \item{yName:}{Outcome variable name corresponds to each predictor-outcome pair}
#' \item{betas:}{Coefficient estimate for predictors }
#' \item{betas_se:}{Coefficent SE for predictors }
#' \item{sigma2:}{Regrssion error terms for predictors }
#' \item{dfs:}{Regression degrees of freedom for predictors }
#' \item{v_g:}{ (X^T X)^{-1} projection on predictors }
#' @importFrom magrittr "%>%"
#' @import tidyr
#' @importFrom tibble as.tibble
#' @import purrr
#' @importFrom matrixStats rowMins
#' @import stats
#' @importFrom methods is
#' @import utils

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
#'                       covariates=covariates, pi1=pi1,
#'                       NoProbXIndex=c(2,4), filter=c(1, -1),
#'                       permutate_number=2, var.ID=c("Gene_ID"),
#'                       grids = seq(0.01, 0.99, by=0.01),fdr = 0.1, PostCut=0.75,
#'                       var.ID.additional=c( "phospho_ID", "Hybridization", "chr"))
#' # Empirical FDR controlled discoveries for all outcomes
#'eFDR=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
#'                  NoProbXIndex=c(2,4),filter=c(1, -1),
#'                  permutate_number=2, var.ID=c("Gene_ID"),
#'                  grids = seq(0.01, 0.99, by=0.01),fdr = 0.1, PostCut=0.75,
#'                   var.ID.additional=c( "phospho_ID", "Hybridization", "chr"), seed=NULL)
#'
#'
iProFun.reg.1y<-function(yList.1y, xList, covariates.1y, permutation=F,
                         var.ID=c("Gene_ID"),
                         var.ID.additional=NULL, seed=NULL){
  # var.ID.additional=c("phospho_ID", "Hybridization", "chr")
  # yList.1y=yList[[1]]; xList=xList; covariates.1y=covariates[[1]];
  # permutation=F; var.ID=c("Gene_ID");
  # var.ID.additional=c("phospho_ID", "Hybridization", "chr"); seed=NULL
  # ----- Data cleaning and quality check ----- #

  # identify varable ID and other non sample variables to save
  
  SubIDExclude=unique(c(var.ID, var.ID.additional))
  Gene_ID_y=yList.1y[var.ID][,1]
  Gene_ID_x=lapply(xList, function(f) f[var.ID][,1])

  # overlapping samples for regression
  ySubID= Reduce(setdiff, SubIDExclude, colnames(yList.1y))
  xSubID=lapply(xList, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
  if (is.null(covariates.1y)) {
    a=append(list(ySubID), xSubID)
    CommonSubID=Reduce(intersect, a)
  } else {
    zSubID=Reduce(setdiff, SubIDExclude, colnames(covariates.1y))
    a=append(append(list(ySubID), xSubID), list(zSubID))
    CommonSubID=Reduce(intersect, a)
  }
  # overlapping variables for regression
  xUnionGeneID=Reduce(union, Gene_ID_x)
  xyCommonGeneID=Reduce(intersect, list(xUnionGeneID, Gene_ID_y))

  # A code for permutation
  if (permutation==T) {
    if(is.null(seed)==F) set.seed(seed);
    CommonSubID_permutate <- sample(CommonSubID)}

  # ----- Linear Regression ----- #
  xlength=length(xList);
  betas=vector("list", xlength); betas_se=vector("list", xlength);
  sigma2=vector("list", xlength); dfs=vector("list", xlength); v_g=vector("list", xlength);
  xName=vector("list", xlength);  yName=vector("list", xlength)

  for (i in 1:length(xyCommonGeneID)) {
    # print(i)
    # if (i%%100==1) {print(i)}

    # print(i)
    if (permutation==T){
      y=yList.1y[which(Gene_ID_y==xyCommonGeneID[i]), CommonSubID_permutate]
    } else {
      y=yList.1y[which(Gene_ID_y==xyCommonGeneID[i]), CommonSubID]
    }
    y=t(y)

    x_i= lapply(1:xlength, function(f) xList[[f]][which(Gene_ID_x[[f]]==xyCommonGeneID[i]), CommonSubID])
    x_index=sapply(x_i, nrow)
    #print(x_index)
    x=t(do.call(rbind,x_i))

    if (is.null(covariates.1y)) {z=NULL} else {
      z= t(matrix(as.matrix(covariates.1y)[, CommonSubID], ncol=length(CommonSubID))) }
    xx=as.matrix(cbind(1, x, z))
    zz = as.matrix(cbind(rep(1, nrow(xx)), z))
    p_xx=ncol(xx)
    p_x=p_xx-ncol(zz)

    # ANOVA to select best Y if multiple Y exist for one variable ID
    if (ncol(y)>1) {
      y_complete = as.matrix(y)[complete.cases(xx), ,drop = F]
      zz_complete = zz[complete.cases(xx), ]
      xx_complete = xx[complete.cases(xx), ]
      yindex=which.min(sapply(1:ncol(y), function(f) anova(lm(y_complete[,f]~xx_complete-1), lm(y_complete[,f]~zz_complete-1))$"Pr(>F)"[2]))
    } else {yindex=1}
    yy=y[,yindex]

    # Self coding of regression
    index=complete.cases(xx) & complete.cases(yy)
    yy_complete2=yy[index];
    xx_complete2=xx[index,];
    n=length(yy_complete2)

    if(p_x==1) {C=as.matrix(c(0, 1, rep(0, p_xx-2)))}
    if (p_x>1) {
      C=matrix(0, nrow=p_xx, ncol=p_x)
      diag(C[-1,])[1:p_x]=1
    }

    ft=lm(yy_complete2~xx_complete2-1)
    naIdx=is.na(ft$coefficients)==F
    # beta
    temp=rep(0, p_xx);
    temp[naIdx]=summary(ft)$coefficients[,1]
    betas_single= t(C) %*% temp
    betas_single[which(betas_single==0)]=NA

    # beta_se
    temp=rep(0, p_xx);
    temp[naIdx]=summary(ft)$coefficients[,2]
    betas_se_single= t(C) %*% temp
    betas_se_single[which(betas_se_single==0)]=NA

    #
    sigma2_single=as.matrix(rep(mean(ft$residuals^2), p_x))
    dfs_single= as.matrix(rep(ft$df.residual, p_x))
    v=vcov(ft)
    v[which(is.na(v))]=0
    v_g_single=as.matrix(diag(t(C) %*% v %*% C))
    v_g_single[which(v_g_single==0)]=NA


    # annotation
    yName_single=yList.1y[which(Gene_ID_y==xyCommonGeneID[i])[yindex],
                            Reduce(intersect, list(SubIDExclude, colnames(yList.1y)))]

    # output as different list for different x platforms  # Here has a problem! check
    for (p in 1:xlength) {
      index=seq(sum(x_index[0:(p-1)])+1,length=x_index[p])

      betas[[p]]=c(betas[[p]], betas_single[index])
      betas_se[[p]]=c(betas_se[[p]], betas_se_single[index])
      sigma2[[p]]=c(sigma2[[p]], sigma2_single[index])
      dfs[[p]]=c(dfs[[p]], dfs_single[index])
      v_g[[p]]=c(v_g[[p]], v_g_single[index])

      xName_single=xList[[p]][which(Gene_ID_x[[p]]==xyCommonGeneID[i]), Reduce(intersect, list(SubIDExclude, colnames(xList[[p]])))]
      xName[[p]] <- rbind(xName[[p]], as.matrix(xName_single))
      if(x_index[p]!=0) {
        yName[[p]] <- rbind(yName[[p]], do.call("rbind", replicate(x_index[p], yName_single, simplify = F)))
      }
    } # end of output list

  } # end of regression

  return(list(yName=yName, xName=xName,betas=betas, betas_se=betas_se,
              sigma2=sigma2, dfs=dfs, v_g=v_g))

} # end of function



#'  Linear regression on all outcome data types (unformatted)
#'
#' Linear regression on all outcome data types with all types of DNA alterations (results not formatted for iProFun input)
#' @export iProFun.reg
#' @param yList yList is a list of data matrix for outcomes.
#' @param xList xList is a list of data matrix for predictors.
#' @param covariates covariates is a list of data matrix for covariate.
#' @param permutation.col permutation.col provides the index of the data types that should be permuated.
#' permutation.col = 0 (default): no permuatation and analysis is on original data. 0 < permutate <= length of yList:
#' permuate the label of the corresponding data type in yList. For example, permutate =2, permute the y label of second
#' data matrix.
#' @param var.ID var.ID gives the variable name (e.g. gene/protein name) to match different data types.
#' If IDs are not specified, the first columns will be considered as ID variable.
#' @param var.ID.additional var.ID.additional  allows to output additional variables from the input.
#' Often helpful if multiple rows (e.g. probes) are considered per gene to allow clear index of results.
#' @param seed seed allows users to externally assign seed to replicate results.
#' @return list with the same length as xlist. Nested within each list, it contains
#' \item{reg.out.list:}{reg.out.list returns the regression summary for each outcome data types as a list. Within each list,
#' see output of iProFun.reg.1y for details.}

iProFun.reg<-function(yList, xList, covariates, permutation.col=0,
                      var.ID=c("Gene_ID"),
                      var.ID.additional=NULL, seed=NULL) {
  ylength=length(yList)
  perm=rep(F, ylength)
  #print(c("permutation.col", permutation.col))
  if (any(permutation.col!=0)) {perm[permutation.col]=T}
  #print(perm)
  reg.out.list=vector("list", ylength)
  for (q in 1:ylength) {

    if (is.null(seed)){seed.q=NULL}
    if (is.null(seed)==F) {seed.q =as.integer((seed+q)*15)}

    reg.out.list[[q]]=iProFun.reg.1y(yList.1y=yList[[q]], xList=xList,
                           covariates.1y=covariates[[q]], permutation=perm[q],
                       var.ID=var.ID, var.ID.additional=var.ID.additional, seed=seed.q)
  }
  return(reg.out.list)
}


#' Reformat regression summaries
#'
#' Reformat regression summaries from iProFun.reg as input of iProFun analysis for association probabilities and discoveries
#' @export multi.omic.reg.summary
#' @param reg.out.list reg.out.list is the regression summaries from iProFun.reg.
#' @param var.ID var.ID gives the variable name (e.g. gene/protein name) to match different data types.
#'
#' @return list with the same length as xlist. Nested within each list, it contains
#' \item{betas_J:}{Coefficient estimate for predictors across J outcome data types}
#' \item{betas_se_J:}{Coefficent SE for predictors across J outcome data types}
#' \item{sigma2_J:}{Regrssion error term for predictors across J outcome data types}
#' \item{dfs_J:}{Regression degrees of freedom for predictors across J outcome data types}
#' \item{v_g_J:}{ (X^T X)^{-1} projection on predictors across J outcome data types}
#' \item{xName_J:}{Predictor name corresponds to each predictor-outcome pair across J outcome data types}
#' \item{yName_J:}{Outcome name corresponds to each predictor-outcome pair across J outcome data types}


multi.omic.reg.summary<-function(reg.out.list, var.ID) {

  ylength=length(reg.out.list)
  xlength=length(reg.out.list[[1]]$betas)

  betas_J <- vector("list", xlength);   betas_se_J <- vector("list", xlength);
  sigma2_J <- vector("list", xlength); dfs_J <- vector("list", xlength);
  v_g_J <- vector("list", xlength);   xName_J <- vector("list", xlength);
  yName_J <- vector("list", xlength)

  for (p in 1:xlength) {

    # creat x_variables for a data type that may or may not be analyzed with a y_variable
    x_ID= unique(do.call("rbind", lapply(1:ylength, function(f) reg.out.list[[f]]$xName[[p]])))
    colnames(x_ID)[1]=var.ID

    for (q in 1:ylength) {
      x_ID_q= reg.out.list[[q]]$xName[[p]]

      temp_betas=temp_betas_se=temp_sigma2=temp_dfs=temp_vg=rep(NA, nrow(x_ID))

      index=row.match(data.frame(x_ID_q), data.frame(x_ID))

      temp_betas[index]=reg.out.list[[q]]$betas[[p]]
      temp_betas_se[index]=reg.out.list[[q]]$betas_se[[p]]
      temp_sigma2[index]=reg.out.list[[q]]$sigma2[[p]]
      temp_dfs[index]=reg.out.list[[q]]$dfs[[p]]
      temp_vg[index]=reg.out.list[[q]]$v_g[[p]]
      temp=as.matrix(reg.out.list[[q]]$yName[[p]])
      if (ncol(temp)>1) {
        temp_yName=matrix(NA, nrow=nrow(x_ID), ncol=ncol(temp)-1)
        temp_yName[index,]=matrix(temp[, !colnames(temp)%in%var.ID], nrow=nrow(x_ID_q))
        yName_J[[p]] <- cbind(yName_J[[p]], temp_yName)
      }

      betas_J[[p]]=cbind(betas_J[[p]], temp_betas)
      betas_se_J[[p]]=cbind(betas_se_J[[p]], temp_betas_se)
      sigma2_J[[p]]=cbind(sigma2_J[[p]], temp_sigma2)
      dfs_J[[p]]=cbind(dfs_J[[p]], temp_dfs)
      v_g_J[[p]]=cbind(v_g_J[[p]], temp_vg)

    } # end q
    xName_J[[p]] <- x_ID
    yName_J[[p]] <- cbind(x_ID[,colnames(x_ID)%in%var.ID], yName_J[[p]])
  } # end p

  return(list(betas_J=betas_J, betas_se_J=betas_se_J, sigma2_J=sigma2_J, dfs_J=dfs_J,
              v_g_J=v_g_J, xName_J=xName_J, yName_J=yName_J))
}


#' Reformat the regression output into a table format. Each row is for a predictor-outcome pair. 
#' @export iProFun.reg.table
#' @param reg.sum Output from multi.omic.reg.summary.
#' @param xType A vector of string for the data types of xList, such as "mutation", "CNV" and "methylation". 
#' @param yType A vector of string for the data types of xList, such as "RNA", "protein" and "phospho". 

#' @return A table that includes gene ID, other gene info if provided, predictor data type, outcome data type, 
#' estimate, se and p-value from Student's t-test, in a long format. 

iProFun.reg.table<-function(reg.all, xType=NULL, yType=NULL, var.ID) {
  reg.sum=multi.omic.reg.summary(reg.out.list=reg.all, var.ID=var.ID)
  xlength=length(reg.sum$betas_J)
  ylength=ncol(reg.sum$betas_J[[1]])
  if(is.null(xType)) {xType=1:xlength}
  if(is.null(yType)) {xType=1:ylength}
  
  output=NULL
  for (i in 1:xlength) {
    for (j in 1:ylength) {
      beta=reg.sum$betas_J[[i]][,j]
      beta_se=reg.sum$betas_se_J[[i]][,j]
      xName=reg.sum$xName_J[[i]]
      yName=reg.sum$yName_J[[i]]
      df=reg.sum$dfs_J[[i]][,j]
      pvalue=pt(abs(beta)/beta_se , df=df, lower=F)*2
      dat=data.frame(xName=xName, yName=yName, xType=xType[i], yType=yType[j], 
                     est=beta, se=beta_se, pvalue=pvalue)
      output=output %>% bind_rows(dat)
    }
  }

  return(output)
}


#' Calculating posterior probabilities of association patterns
#'
#' This is the core function of iProFun, and it calculates the posterior probabilities of association patterns from
#' regression summaries from multiple data types. It allows some data types from regression to opt out from calculating
#' posterior probabilities, which is suitable for the data types with very few genes (e.g. somatic mutation) that
#' borrowing information across variables is not reliable.
#'@export iProFun.prob
#' @param Reg.Sum Linear regression analysis summaries formatted in multi.omic.reg.summary.
#' @param NoProbXIndex NoProbXIndex allows users to provide the index for the predictor data type(s) that are not considered
#' for calculating posterior probabilities of association patterns. NoProbXIndex = NULL (default): all data types in
#' predictors will be considered. Any 0 < NoProbXIndex <= length of xList:
#' indicates that the posterior probabilities of association patterns for the corresponding data type(s) are not calculated.
#' @param pi1 pi1 is the pre-specified prior proportion of non-null statistics. It cane be a number in (0, 1) or a vector
#' of numbers with length of ylist.
#'
#'
#' @import parallel
#' @import stats
#' @import metRology
#' @import prodlim
#' @return A list with the same length as xlist. Nested within each list, it contains
#' \item{NoComputation:}{No of variables considered for this predictor with all outcomes}
#' \item{Config:}{Corresponding association patterns. Total number 2^J for J outcome data types}
#' \item{Config.miss:}{Corresponding association patterns for variables (e.g. genes) that are missing in some outcome data type(s)}
#' \item{PostProb:}{Posterior probability for each predictor on each association pattern}
#' \item{PostProb.miss:}{Posterior probability for each predictor on each available association pattern  for variables (e.g. genes) that are missing in some outcome data type(s)}
#' \item{xName.miss:}{Predictor name for variables (e.g. genes) missing in some data type(s)}
#' \item{colocProb:}{Averaged posterior probability for predictor on each association pattern (missing data excluded)}
#' \item{Tstat_L:}{T statistics for each predictor on each outcome}
#' \item{D0:}{Estimated density under the null for predictor on each outcome}
#' \item{D1:}{Estimated density under the alternative for predictor on each outcome}

iProFun.prob = function(Reg.Sum, NoProbXIndex=NULL, pi1=0.05) {

  xlength=length(Reg.Sum$betas_J)
  ylength=ncol(Reg.Sum$betas_J[[1]])
  x.prob.index=setdiff(seq(1:xlength), NoProbXIndex)

  if (length(pi1)==1) {
    pi1=matrix(pi1, ncol=ylength, nrow=xlength)
  }

  Reg_output <- vector("list", xlength);
  x_iProFun <- vector("list", xlength);
  final_result <- vector(mode="list", xlength)


  for (p in x.prob.index) {

    # Summarize Regression
    Reg_output[[p]]= list(betas_J=Reg.Sum$betas_J[[p]], betas_se_J=Reg.Sum$betas_se_J[[p]],
                          sigma2_J=Reg.Sum$sigma2_J[[p]],dfs_J=Reg.Sum$dfs_J[[p]],
                          v_g_J = Reg.Sum$v_g_J[[p]], xName_J = Reg.Sum$xName_J[[p]],
                          yName_J = Reg.Sum$yName_J[[p]])
    # run iProFun
    x_iProFun[[p]]= iProFun.1x.prob(input=Reg_output[[p]],pi1 = pi1[p,])

    # Summarize output
    names(final_result)[p] <- paste0("iProFun output for xlist ", p)
    final_result[[p]] <- vector(mode="list", length=9)
    final_result[[p]] = append(Reg_output[[p]], x_iProFun[[p]])
  }
  return(final_result)
}


#' Calculate FWER 
#'
#' This function calculates FWER across all the variables for the same predictor data type (e.g. mutation), directly
#' from regression results without considering posterior probabilities of association patterns.
#' This strategy is preferred for data types that have few variables that cannot reliably infer association patterns.
#' @export iProFun.FWER
#' @param reg.all Linear regression analysis summaries from iProFun.reg.
#' @param FWER.Index  Index the predictor data types that calculate FWER directly
#' @return    It returns original p-values and FWER with and without filter
#' \item{pvalue:}{The p-value of linear regression with Student's t distribution}
#' \item{FWER:}{The FWER with Bonferroni correction for all genes}
#' \item{xName:}{The gene name of the predictors}

#' @import stats

iProFun.FWER= function(reg.all, FWER.Index=0, var.ID) {
  Reg.Sum=multi.omic.reg.summary(reg.out.list=reg.all, var.ID=var.ID)
  # examples: mutation that we only have few genes, and would like to calculate FWER
  # without iProFun posterior prob calculation.

  l=length(FWER.Index)
  # add filter

  final_result <- vector(mode="list", l)
  
  for (j in 1:l) {
    k=FWER.Index[j]

    Reg_output= list(betas_J=Reg.Sum$betas_J[[k]], betas_se_J=Reg.Sum$betas_se_J[[k]],
                          sigma2_J=Reg.Sum$sigma2_J[[k]],dfs_J=Reg.Sum$dfs_J[[k]],
                          v_g_J = Reg.Sum$v_g_J[[k]], xName_J = Reg.Sum$xName_J[[k]],
                          yName_J = Reg.Sum$yName_J[[k]])
    betas_J=Reg_output$betas_J

    t=betas_J/Reg_output$betas_se_J
    pvalue=sapply(1:ncol(t), function(f) sapply(1:nrow(t), function(g)
      pt(abs(t[g,f]), df=Reg_output$dfs_J[g,f], lower.tail = F)*2 ))

    FWER=apply(pvalue, 2, function(f) p.adjust(f, method="bonferroni"))
   
    result=list(pvalue=pvalue, FWER=FWER, xName=Reg_output$xName_J)
    final_result[[j]]= result
  }
  return(final_result)

}







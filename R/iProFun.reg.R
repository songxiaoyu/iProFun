
#' Linear regression on one outcome data type
#'
#'  Linear regression on one outcome data type with all data types of DNA-level alterations.
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
#' FWER.all=iProFun.FWER(reg.all=reg.all, FWER.Index=c(1))
#' # Calculate Empirical FDR for one outcome
#' eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=2, yList=yList, xList=xList,
#'                       covariates=covariates, pi1=pi1, NoProbXIndex=c(1),
#'                       permutate_number=2, var.ID=c("geneSymbol"),
#'                       var.ID.additional=c("id"))
#' # Calculate Empirical FDR for all outcomes
#'eFDR.all=iProFun.eFDR(reg.all=reg.all, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
#'                  NoProbXIndex=c(1),
#'                  permutate_number=2, var.ID=c("geneSymbol"),
#'                   var.ID.additional=c( "id"), seed=123)
#' # iProFun identification
#' # For data types with abundance genes, it's based on (1) association probabilities > 0.75,
#' # (2) FDR 0.1, and (3) the association direction filtering.
#' # For data types with few genes, it's  based (1) FWER 0.1, and
#' # (2)  the association direction filtering.
#'
#'res=iProFun.detection(reg.all=reg.all, eFDR.all=eFDR.all, FWER.all=FWER.all, filter=c(0, 1),
#'                      NoProbButFWERIndex=1,fdr.cutoff = 0.1, fwer.cutoff=0.1, PostPob.cutoff=0.75,
#'                      xType=c("mutation", "cnv"), yType=c("rna", "protein", "phospho"))



iProFun.reg.1y<-function(yList.1y, xList, covariates.1y, permutation=F,
                         var.ID=c("Gene_ID"),
                         var.ID.additional=NULL, seed=NULL){

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



#'  Linear regression on all outcome data types
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
#' @return list with the same length as xList. Nested within each list, it contains
#' \item{betas_J:}{Coefficient estimate for predictors across J outcome data types}
#' \item{betas_se_J:}{Coefficent SE for predictors across J outcome data types}
#' \item{sigma2_J:}{Regrssion error term for predictors across J outcome data types}
#' \item{dfs_J:}{Regression degrees of freedom for predictors across J outcome data types}
#' \item{v_g_J:}{ (X^T X)^{-1} projection on predictors across J outcome data types}
#' \item{xName_J:}{Predictor name corresponds to each predictor-outcome pair across J outcome data types}
#' \item{yName_J:}{Outcome name corresponds to each predictor-outcome pair across J outcome data types}


multi.omic.reg.summary<-function(reg.out.list) {

  ylength=length(reg.out.list)
  xlength=length(reg.out.list[[1]]$betas)

  betas_J <- vector("list", xlength);   betas_se_J <- vector("list", xlength);
  sigma2_J <- vector("list", xlength); dfs_J <- vector("list", xlength);
  v_g_J <- vector("list", xlength);   xName_J <- vector("list", xlength);
  yName_J <- vector("list", xlength)

  for (p in 1:xlength) {

    # creat x_variables for a data type that may or may not be analyzed with a y_variable
    x_ID= unique(do.call("rbind", lapply(1:ylength, function(f) reg.out.list[[f]]$xName[[p]])))
    # colnames(x_ID)[1]=var.ID

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
        temp_yName=matrix(NA, nrow=nrow(x_ID), ncol=ncol(temp))
        temp_yName[index,]=matrix(temp, nrow=nrow(x_ID_q))
        yName_J[[p]] <- cbind(yName_J[[p]], temp_yName)
      }

      betas_J[[p]]=cbind(betas_J[[p]], temp_betas)
      betas_se_J[[p]]=cbind(betas_se_J[[p]], temp_betas_se)
      sigma2_J[[p]]=cbind(sigma2_J[[p]], temp_sigma2)
      dfs_J[[p]]=cbind(dfs_J[[p]], temp_dfs)
      v_g_J[[p]]=cbind(v_g_J[[p]], temp_vg)

    } # end q
    xName_J[[p]] <- x_ID
    yName_J[[p]] <- t(unique(t(yName_J[[p]])))
  } # end p

  return(list(betas_J=betas_J, betas_se_J=betas_se_J, sigma2_J=sigma2_J, dfs_J=dfs_J,
              v_g_J=v_g_J, xName_J=xName_J, yName_J=yName_J))
}


#' Reformat the regression output into a table format. Each row is for a predictor-outcome pair.
#' @export iProFun.reg.table
#' @param reg.all Output from iProFun.reg.
#' @param xType A vector of string for the data types of xList, such as "mutation", "CNV" and "methylation".
#' @param yType A vector of string for the data types of xList, such as "RNA", "protein" and "phospho".

#' @return A table that includes gene ID, other gene info if provided, predictor data type, outcome data type,
#' estimate, se and p-value from Student's t-test, in a long format.

iProFun.reg.table<-function(reg.all, xType=NULL, yType=NULL) {
  reg.sum=multi.omic.reg.summary(reg.out.list=reg.all)
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
      pvalue=pt(abs(beta)/beta_se , df=df, lower.tail=F)*2
      dat=data.frame(xName=xName, yName=yName, xType=xType[i], yType=yType[j],
                     est=beta, se=beta_se, pvalue=pvalue)
      output=output %>% dplyr::bind_rows(dat)
    }
  }

  return(output)
}


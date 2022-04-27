#' Calculate FWER
#'
#' This function calculates FWER across all the outcome variables for the same predictor data type
#' (e.g. mutation), directly from regression results without considering posterior probabilities of association patterns.
#' This strategy is preferred for data types that have few genes that cannot reliably infer association patterns.
#' @export iProFun.FWER
#' @param reg.all Linear regression analysis summaries from iProFun.reg.
#' @param FWER.Index  Index the predictor data types that calculate FWER directly
#' @return    It returns original p-values and FWER with and without filter
#' \item{pvalue:}{The p-value of linear regression with Student's t distribution}
#' \item{FWER:}{The FWER with Bonferroni correction for all genes}
#' \item{xName:}{The gene name of the predictors}
#' @import stats

iProFun.FWER= function(reg.all, FWER.Index=0) {
  Reg.Sum=multi.omic.reg.summary(reg.out.list=reg.all)
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





#' iProFun eFDR assessment for one outcome data type.
#'
#'iProFun empirical false discovery rate (eFDR) assessment for one outcome data type.
#' @export iProFun.eFDR.1y
#' @param yList yList is a list of data matrix for outcomes.
#' @param xList xList is a list of data matrix for predictors.
#' @param covariates covariates is a list of data matrix for covariate.
#' @param pi1 pi1 is pre-specified prior of proportion of non-null statistics. It cane be a number in (0, 1) or a vector
#' of numbers with length of ylist.
#' @param var.ID var.ID gives the variable name (e.g. gene/protein name) to match different data types.
#' If IDs are not specified, the first columns will be considered as ID variable.
#' @param var.ID.additional var.ID.additional allows to output additional variable names from the input.
#' Often helpful if multiple rows (e.g. probes) are considered per gene to allow clear index of the rows.
#' @param permutate_number Number of permutation, default 10
#' @param reg.all The regression summary (unformatted) such as from iProFun.reg.
#' @param which.y The index for one outcome data type that the eFDR assessment is for.
#' @param NoProbXIndex NoProbXIndex allows users to provide the index for the predictor data type(s) that are not considered
#' for calculating posterior probabilities of association patterns.
#' @param seed seed allows users to externally assign seed to replicate results.
#' @return
#' \item{eFDR.grid:}{eFDR by the grid of posterior probability cutoffs.}
#' \item{fdr_cutPob:}{the cutoff values for pre-specified eFDR rate and the posterior probabilities for a pair of
#'  data types based on permutation.}
#' \item{No.Identified.filter:}{the number of identified variables for each  pair of data types.}
#' \item{No.Identified.no.filter:}{the number of identified variables for each  pair of data types.}
#' \item{Gene_fdr:}{A table summarizing the posterior probabilities (PostProb), the eFDR (eFDR.no.filter),
#' the significance under different criteria (nominal FDR, PostProb cutoffs and filter)  for each variable (e.g. gene)
#' under consideration.}
#'



iProFun.eFDR.1y= function(reg.all, which.y, yList, xList, covariates, pi1,
                          NoProbXIndex=NULL,
                          permutate_number=10,
                          var.ID=c("Gene_ID"),
                          var.ID.additional=NULL, seed=NULL) {

  xlength=length(xList)
  x.prob=setdiff(1:xlength, NoProbXIndex)
  xlength2=length(x.prob)

  # Obtain posterior probability from original data - including missing
  sum=multi.omic.reg.summary(reg.out.list=reg.all)
  prob=iProFun.prob(Reg.Sum=sum, NoProbXIndex=NoProbXIndex, pi1=pi1)
  prob1y=iProFun.sum.prob.1y(prob=prob, which.y=which.y, NoProbXIndex=NoProbXIndex)

  # Obtain posterior probability from permutation data
  perm.reg.all=reg.all
  permProb1y=vector("list", xlength2)
  for (perm in 1: permutate_number) {
    print(c("perm", perm))
    if (is.null(seed)){seed.p=NULL}
    if (is.null(seed)==F) {seed.p =as.integer((seed+perm)*978)}

     ft1=iProFun.reg.1y(yList.1y=yList[[which.y]], xList=xList,
                        covariates.1y=covariates[[which.y]], permutation=T,
                        var.ID=var.ID, var.ID.additional=var.ID.additional, seed=seed.p)

    perm.reg.all[[which.y]]=ft1
    sum1=multi.omic.reg.summary(reg.out.list=perm.reg.all)
    prob1=iProFun.prob(Reg.Sum=sum1, NoProbXIndex=NoProbXIndex, pi1=pi1)
    perm1y=iProFun.sum.prob.1y(prob=prob1, which.y=which.y, NoProbXIndex=NoProbXIndex)

    for (p in 1:xlength2) {
      permProb1y[[p]]=cbind(permProb1y[[p]], perm1y[[p]])
    }
    if( perm%%10 ==0) {print(paste0("Completed Permutation ", perm-9, " to ", perm))}
  }

  # gene eFDR
  Gene_efdr=vector(mode="list", length=2)
  for (p in 1:xlength2) {
    cut.max=1-1e-8
    cut.grid=ifelse(prob1y[[p]]<1, prob1y[[p]], cut.max )
    Gene_efdr[[p]]=    sapply(1:length(prob1y[[p]]), function(f)
      mean(permProb1y[[p]]>cut.grid[f], na.rm=T)/mean(prob1y[[p]]>cut.grid[f], na.rm=T))
  }
  # name
  xName=lapply(1:xlength2, function(p)
    prob[[x.prob[p]]]$xName_J)

  return(list(xName=xName, PostProb=prob1y, Gene_efdr=Gene_efdr ))
}


#' iProFun eFDR assessment based on permutation for multiple outcome data type.
#'
#' iProFun empirical false discovery rate (eFDR) assessment based on permutation for multiple outcome data types.

#' @export iProFun.eFDR
#' @param yList yList is a list of data matrix for outcomes.
#' @param xList xList is a list of data matrix for predictors.
#' @param covariates covariates is a list of data matrix for covariate.
#' @param pi1 pi1 is pre-specified prior of proportion of non-null statistics. It cane be a number in (0, 1) or a vector
#' of numbers with length of ylist.
#' @param var.ID var.ID gives the variable name (e.g. gene/protein name) to match different data types.
#' If IDs are not specified, the first columns will be considered as ID variable.
#' @param var.ID.additional var.ID.additional allows to output additional variable names from the input.
#' Often helpful if multiple rows (e.g. probes) are considered per gene to allow clear index of the rows.
#' @param permutate_number Number of permutation, default 10
#' @param reg.all The regression summary (unformatted) such as from iProFun.reg.
#' @param NoProbXIndex NoProbXIndex allows users to provide the index for the predictor data type(s) that are not considered
#' for calculating posterior probabilities of association patterns.
#' @param seed seed allows users to externally assign seed to replicate results.
#' @return
#' \item{xName:}{Name of the predictors.}
#' \item{PostProb:}{The association probability for each gene on each data type.}
#' \item{Gene_efdr:}{The eFDR for each gene on each data type.}


iProFun.eFDR= function(reg.all, yList, xList, covariates, pi1,
                       NoProbXIndex=NULL,
                       permutate_number=10,
                       var.ID=c("Gene_ID"),
                       var.ID.additional=NULL, seed=NULL) {

  ylength=length(yList)
  xlength=length(xList)
  x.prob=setdiff(1:xlength, NoProbXIndex)
  xlength2=length(x.prob)

  PostProb=Gene_efdr=xName=vector("list", xlength2);

  for (q in 1:ylength) {
    print(c("Outcome", q))


    if (is.null(seed)){seed.qq=NULL}
    if (is.null(seed)==F) {seed.qq =as.integer((seed+q)*4109)}

    eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=q, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
                          NoProbXIndex=NoProbXIndex, permutate_number=permutate_number,
                          var.ID=var.ID, var.ID.additional=var.ID.additional, seed=seed.qq)
    for (p in 1:xlength2 ) {
      PostProb[[p]]= cbind(PostProb[[p]], eFDR1$PostProb[[p]]);
      Gene_efdr[[p]]= cbind(Gene_efdr[[p]], eFDR1$Gene_efdr[[p]]);
      xName[[p]]=t(unique(t(cbind(xName[[p]], eFDR1$xName[[p]]))));
    }
  }
  eFDR_result=list(xName=xName, PostProb=PostProb, Gene_efdr=Gene_efdr)
  return(eFDR_result)
}


#' iProFun detection using multiple criteria

#' @export iProFun.detection
#' @param reg.all reg.all is the regression results for all outcomes from `iProFun.reg`.
#' @param eFDR.all eFDR.all is eFDR results for all outcomes from `iProFun.eFDR`.
#' @param FWER.all FWER.all is FWER results for all outceoms from `iProFun.FWER`.
#' @param NoProbButFWERIndex NoProbXIndex allows users to provide the index for the predictor data type(s) that are not considered
#' for calculating posterior probabilities of association patterns. Default=NULL. All predictors go throught posterior
#' probability calculation.
#' @param fdr.cutoff Nominal false discover rate, default at 0.1
#' @param fwer.cutoff family wise error rate cutoff to be considered as
#' significant, default at 0.1.
#' @param PostPob.cutoff Minimal posterior probability cutoff to be considered as
#' significant, default at 0.75.
#' @param filter filter is a vector with the same length of xList, taking values of 1, -1, 0 or NULL.
#' "NULL" is default imposes no filtering. 1" indicates that an association is considered for
#' significance only if its significant associations are positive across all outcome platforms. "-1" indicates
#' that an association is considered
#' for significance only if its significant associations are negative  across all outcome platforms. "0" indicates
#' that an association is considered for significance only if its significant association across all outcome
#' platforms preserve consistent directions (either positive or negative).
#' @param xType A vector of string for the data types of xList, such as "mutation", "CNV" and "methylation", for output presentation.
#' @param yType A vector of string for the data types of xList, such as "RNA", "protein" and "phospho", for output presentation.

#' @return A output table that includes gene ID, other gene info if provided, predictor data type, outcome data
#' type, estimate, se, p-value from Student's t-test, FWER, Posterior Association Probability, eFDR,
#' directional filtering, and whether it's identified in iProFun or not, in a long format.
iProFun.detection=function(reg.all, eFDR.all, FWER.all, filter=c(0, 0), NoProbButFWERIndex=NULL,
                  fdr.cutoff = 0.1, fwer.cutoff=0.1, PostPob.cutoff=0.75,
                  xType, yType) {

  reg.sum=multi.omic.reg.summary(reg.out.list=reg.all)

  xlength=length(reg.sum$betas_J)
  ylength=ncol(reg.sum$betas_J[[1]])
  x.prob.idx=setdiff(1:xlength, NoProbButFWERIndex)
  x.fwer.idx=NoProbButFWERIndex
  xlength2=length(x.prob.idx)


  output=NULL
  # for each x across all y's
  for (p in x.fwer.idx) {

    betas=betas_filter=reg.sum$betas_J[[p]]
    fwer=FWER.all[[p]]$FWER
    betas_filter[which(fwer>=fwer.cutoff)]=NA
    # all significant ones are positive association
    temp1=which(sapply(1:nrow(betas), function(f)
      all(betas_filter[f,]>0, na.rm=T) & all(is.na(betas_filter[f,]))==F ))
    # all significant ones are negative association
    temp2=which(sapply(1:nrow(betas), function(f)
      all(betas_filter[f,]<0, na.rm=T) & all(is.na(betas_filter[f,]))==F ))

    # filter
    if (is.null(filter[p]) ){ # no requirement
      x_filter_gene=seq(1, nrow(betas))
    } else {
      if (filter[p] == 1) {x_filter_gene= temp1} # all positive beta among significant results
      if (filter[p] == -1) { x_filter_gene=temp2} # all negative beta among significant results
      if (filter[p] == 0) {x_filter_gene = sort( union(temp1, temp2))} # all positive or all negative
    }
    betas_filter[-x_filter_gene,]=NA

    # save output
    for (q in 1:ylength) {
      beta=betas[,q]
      beta_se=reg.sum$betas_se_J[[p]][,q]
      xName=reg.sum$xName_J[[p]]
      yName=reg.sum$yName_J[[p]]
      df=reg.sum$dfs_J[[p]][,q]
      pvalue=pt(abs(beta)/beta_se , df=df, lower.tail=F)*2
      d.filter=ifelse(seq(1, nrow(betas)) %in% x_filter_gene, 1, 0)
      iProFun.identification=ifelse(is.na(betas_filter), 0, 1)

      dat=data.frame(xName=xName, yName=yName, xType=xType[p], yType=yType[q],
                     est=beta, se=beta_se, pvalue=pvalue, FWER=fwer[,q],
                     eFDR=NA, PostProb=NA, d.filter=d.filter, iProFun.identification=iProFun.identification[,q])
      output=output %>% dplyr::bind_rows(dat)
    }

  }

  for (p in 1:xlength2) {
    k=x.prob.idx[p]
    betas=betas_filter=reg.sum$betas_J[[k]]
    eFDR=eFDR.all$Gene_efdr[[p]]
    PostProb=eFDR.all$PostProb[[p]]
    betas_filter[which(eFDR>=fdr.cutoff)]=NA
    betas_filter[which(PostProb<=PostPob.cutoff)]=NA

    # all significant ones are positive association
    temp1=which(sapply(1:nrow(betas), function(f)
      all(betas_filter[f,]>0, na.rm=T) & all(is.na(betas_filter[f,]))==F ))
    # all significant ones are negative association
    temp2=which(sapply(1:nrow(betas), function(f)
      all(betas_filter[f,]<0, na.rm=T) & all(is.na(betas_filter[f,]))==F ))
    # filter
    if (is.null(filter[k]) ){ # no requirement
      x_filter_gene=seq(1, nrow(betas))
    } else {
      if (filter[k] == 1) {x_filter_gene= temp1} # all positive beta among significant results
      if (filter[k] == -1) { x_filter_gene=temp2} # all negative beta among significant results
      if (filter[k] == 0) {x_filter_gene = sort( union(temp1, temp2))} # all positive or all negative
    }
    betas_filter[-x_filter_gene,]=NA
    # save output
    for (q in 1:ylength) {
      beta=betas[,q]
      beta_se=reg.sum$betas_se_J[[k]][,q]
      xName=reg.sum$xName_J[[k]]
      yName=reg.sum$yName_J[[k]]
      df=reg.sum$dfs_J[[k]][,q]
      pvalue=pt(abs(beta)/beta_se , df=df, lower.tail=F)*2
      d.filter=ifelse(seq(1, nrow(betas)) %in% x_filter_gene, 1, 0)
      iProFun.identification=ifelse(is.na(betas_filter), 0, 1)

      dat=data.frame(xName=xName, yName=yName, xType=xType[k], yType=yType[q],
                     est=beta, se=beta_se, pvalue=pvalue, FWER=NA,
                     eFDR=eFDR[,q], PostProb=PostProb[,q], d.filter=d.filter, iProFun.identification=iProFun.identification[,q])
      output=output %>% dplyr::bind_rows(dat)
    }


  }

  output.narm=output[-which(is.na(output$est)),]
  return(output.narm)

}








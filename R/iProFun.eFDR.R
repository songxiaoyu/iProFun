#' Summarize the association probabilities for an outcome.
#'
#' Summarize the association probabilities for an outcome across different association patterns.
#' The summation includes variables with missingness in some data types.
#' @export iProFun.sum.prob.1y
#' @param prob The posterior probabilieis of association patterns across multiple data types,
#' including association patterns with missing data in some of the data types. The output of iProFun.prob can be used here.
#' @param which.y The index for one outcome data type that the summation is for.
#' @param NoProbXIndex  The index for the predictor data type(s) that are not considered
#' for calculating posterior probabilities of association patterns.
#' @return
#' \item{prob.all:}{The summarized posterior probabilities for all outcome variables in the data type of interest
#' (regardless of missingness in other data types).}
#'
#'

# function to get prob for 1y, with missing
iProFun.sum.prob.1y <- function (prob, which.y=2, NoProbXIndex=c(2,4)) {
  xlength=length(prob)
  x.prob=setdiff(1:xlength, NoProbXIndex)

  # calculate prob on 1y
  prob.all <- vector("list", length(x.prob))

  for (p in 1:length(x.prob)) {
    k= x.prob[p]
    Q_index=prob[[k]]$Config[,which.y]
    prob.all[[p]]=apply(prob[[k]]$PostProb[,which(Q_index==1)], 1, sum)

    config.miss=prob[[k]]$Config.miss
    for (j in 1:length(config.miss)) {
      #print(j)
      Q_index1= config.miss[[j]][, colnames(config.miss[[j]])==as.character(which.y)]
      prob.temp=prob[[k]]$PostProb.miss[[j]]

      if (length(which(Q_index1==1))==1 ) {
        prob.miss1=prob.temp[,which(Q_index1==1)]
        prob.all[[p]][row.match(data.frame(prob[[k]]$xName.miss[[j]]), data.frame(prob[[k]]$xName_J))]=prob.miss1
      }

      if (length(which(Q_index1==1))>1 &  nrow(prob.temp)>1) {
        prob.miss1=apply(as.matrix(prob.temp[,which(Q_index1==1)]), 1, sum)
        prob.all[[p]][row.match(data.frame(prob[[k]]$xName.miss[[j]]), data.frame(prob[[k]]$xName_J))]=prob.miss1
      }

    } # end (j), which fills the list of missing values
  }
  return(prob.all)
}



#' Fast estimation of FDR without permutation procedure for one data type
#'
#' @export estFDR
#' @param ProbPattern1y1x Posterior probabilities for associations between one predictor data type and one outcome data type.
#' @param grids grids specify the searching grids to find significant associations
#' @return
#' \item{estFDR:}{Estimated FDR values by different posterior probability cutoffs considered in grids}

estFDR<-function(ProbPattern1y1x, grids = seq(0.01, 0.99, by=0.01)) {
  # estFDR=  sum_i (1-Prob_i) 1(Prob_i>lambda) over No. Prob_i>lambda
  estFDR=NULL
  for (i in 1:length(grids) ) {
    denominator= sum(ProbPattern1y1x>grids[i], na.rm=T)
    numerator = sum ( (1-ProbPattern1y1x)*(ProbPattern1y1x>grids[i]), na.rm=T)
    estFDR=c(estFDR, numerator/denominator)
  }
  return(estFDR)
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
#' @param fdr norminal false discover rate, default at 0.1
#' @param PostCut PostCut specifies minimal posterior probabilty cutoff to be considered as
#' significant, default as 0.75
#' @param filter filter is a vector with the same length of xList, taking values of 1, -1, 0 or NULL.
#' "NULL" is default imposes no filtering. 1" indicates that an association is considered for
#' significance only if its significant associations are positive across all outcome platforms. "-1" indicates
#' that an association is considered
#' for significance only if its significant associations are negative  across all outcome platforms. "0" indicates
#' that an association is considered for significance only if its significant association across all outcome
#' platforms preserve consistent directions (either positive or negative).
#' @param grids Grids specify the searching grids for significant associations
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
                          NoProbXIndex=NULL,filter=NULL,
                          permutate_number=10,
                          grids = seq(0.01, 0.99, by=0.01),
                          fdr = 0.1, PostCut=0.75,
                          var.ID=c("Gene_ID"),
                          var.ID.additional=NULL, seed=NULL) {

  # NoProbXIndex=c(2,4);filter=NULL;
  # permutate_number=1;
  # grids = seq(0.01, 0.99, by=0.01);
  # fdr = 0.1; PostCut=0.75;
  # var.ID=c("Gene_ID");
  # var.ID.additional=c( "phospho_ID", "Hybridization", "chr"); seed=NULL

  xlength=length(xList)
  x.prob=setdiff(1:xlength, NoProbXIndex)
  xlength2=length(x.prob)


  # Obtain posterior probabiliy from original data - including missing
  sum=multi.omic.reg.summary(reg.out.list=reg.all, var.ID=var.ID)
  prob=iProFun.prob(Reg.Sum=sum, NoProbXIndex=NoProbXIndex, pi1=pi1)
  prob1y=iProFun.sum.prob.1y(prob=prob, which.y=which.y, NoProbXIndex=NoProbXIndex)
  count_orig_grid=lapply(1: xlength2, function(p)
    sapply(1:length(grids), function(f)
      sum(prob1y[[p]]>grids[f], na.rm=T)))


  # permutation
  perm.reg.all=reg.all
  permProb1y=vector("list", xlength2)
  for (perm in 1: permutate_number) {
    print(c("perm", perm))
    if (is.null(seed)){seed.p=NULL}
    if (is.null(seed)==F) {seed.p =as.integer((seed+perm)*978)}

     ft1=iProFun.reg.1y(yList.1y=yList[[which.y]], xList=xList,  covariates.1y=covariates[[which.y]], permutation=T,
                            var.ID=var.ID,
                            var.ID.additional=var.ID.additional, seed=seed.p)


    perm.reg.all[[which.y]]=ft1
    sum1=multi.omic.reg.summary(reg.out.list=perm.reg.all, var.ID=var.ID)
    prob1=iProFun.prob(Reg.Sum=sum1, NoProbXIndex=NoProbXIndex, pi1=pi1)
    perm1y=iProFun.sum.prob.1y(prob=prob1, which.y=which.y, NoProbXIndex=NoProbXIndex)

    for (p in 1:xlength2) {
      permProb1y[[p]]=cbind(permProb1y[[p]], perm1y[[p]])
    }
    if( perm%%10 ==0) {print(paste0("Completed Permutation ", perm-9, " to ", perm))}
  }

  # calculate the number of genes significant at each threshold level
  count_perm_grid=lapply(1:xlength2, function(p)
      sapply(1:length(grids), function(f) sum(permProb1y[[p]]>grids[f], na.rm=T)/permutate_number ))

  # calculate the eFDR on a grid of prob
  eFDR.grid=lapply(1:xlength2, function(p)
    count_perm_grid[[p]]/count_orig_grid[[p]])

  AboveCut=which(grids>=PostCut)
  grid_PostCut= AboveCut[which.min(abs(AboveCut- PostCut))]

  fdr_cut=lapply(eFDR.grid, function(f)
    max(min(which(f<fdr), na.rm=T), grid_PostCut))

  fdr_cutPob=lapply(1:xlength2, function(f) grids[fdr_cut[[f]]] )

  No.Identified.no.filter=lapply(1:xlength2, function(f)
    count_orig_grid[[f]]
    [fdr_cut[[f]] ])

  # calculate eFDR based on cutoff
  eFDR.no.filter=vector("list", xlength2);
  for (f in 1:xlength2) {
    t=length(prob1y[[f]])
    t2=rep(NA, t)
    t2[which(is.na(prob1y[[f]])==F)]=sapply(which(is.na(prob1y[[f]])==F), function(g) tail(which(prob1y[[f]][g]>grids), n=1))
    t2[sapply(t2, function(x) length(x)==0)] <- NA
    t2=unlist(t2)
    eFDR.no.filter[[f]]=eFDR.grid[[f]][t2]
    eFDR.no.filter[[f]][which(is.na(prob1y[[f]])==F & is.na(t2))] =1
  }
  # add filter
  x_filter_gene=vector("list", xlength2)
  for (j in 1: xlength2){
    k=x.prob[j]
    betas_J=prob[[k]]$betas_J

  if (is.null(filter[j]) ){ # no requirement
    x_filter_gene[[j]]=seq(1, nrow(betas_J))
  } else {

    temp1=which(sapply(1:nrow(betas_J), function(f)
      any(betas_J[f,][(eFDR.no.filter[[j]][f]<fdr)]<=0, na.rm=T))==F)

    temp2=which(sapply(1:nrow(betas_J), function(f)
      any(betas_J[f,][(eFDR.no.filter[[j]][f]<fdr)]>=0, na.rm=T))==F)

    if (filter[j] == 1) {x_filter_gene[[j]]= temp1} # all positive beta among signifiant results
    if (filter[j] == -1) { x_filter_gene[[j]]=temp2} # all negative beta among signifiant results
    if (filter[j] == 0) {x_filter_gene[[j]] = sort( union(temp1, temp2))} # all positive or all negative
    }
  }

  No.Identified.filter=lapply(1: xlength2, function(p)
      sum(prob1y[[p]][x_filter_gene[[p]]]>fdr_cutPob[[p]], na.rm=T))

  Gene_fdr=vector("list", xlength2);
  for (j in 1:xlength2) {
    k=x.prob[j]
    # no filter
    temp1=temp2=rep(0, nrow(prob[[k]]$xName_J))
    temp1[is.na(prob1y[[j]])]=temp2[is.na(prob1y[[j]])]=NA
    temp1[which(prob1y[[j]]>fdr_cutPob[[j]])]=1

    # filter
    sig=intersect(x_filter_gene[[j]], which(prob1y[[j]]>fdr_cutPob[[j]]))
    temp2[sig]=1

    Gene_fdr[[j]]=data.frame(xName=prob[[k]]$xName_J, PostProb=prob1y[[j]], eFDR.no.filter=eFDR.no.filter[[j]], sig.no.filter=temp1, sig.filter=temp2)
  }

  eFDR_result=list(eFDR.grid=eFDR.grid, fdr_cutPob=fdr_cutPob, No.Identified.filter=No.Identified.filter, No.Identified.no.filter=No.Identified.no.filter, Gene_fdr=Gene_fdr)

  return(eFDR_result)
}

#' iProFun eFDR assessment based on permutation for multiple outcome data type.
#'
#' iProFun empirical false discovery rate (eFDR) assessment based on permutation for multiple outcome data type.

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
#' @param fdr norminal false discover rate, default at 0.1
#' @param PostCut PostCut specifies minimal posterior probabilty cutoff to be considered as
#' significant, default as 0.75
#' @param filter filter is a vector with the same length of xList, taking values of 1, -1, 0 or NULL.
#' "NULL" is default imposes no filtering. 1" indicates that an association is considered for
#' significance only if its significant associations are positive across all outcome platforms. "-1" indicates
#' that an association is considered
#' for significance only if its significant associations are negative  across all outcome platforms. "0" indicates
#' that an association is considered for significance only if its significant association across all outcome
#' platforms preserve consistent directions (either positive or negative).
#' @param grids Grids specify the searching grids for significant associations
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



iProFun.eFDR= function(reg.all, yList, xList, covariates, pi1,
                          NoProbXIndex=NULL,filter=NULL,
                          permutate_number=10,
                          grids = seq(0.01, 0.99, by=0.01),
                          fdr = 0.1, PostCut=0.75,
                          var.ID=c("Gene_ID"),
                          var.ID.additional=c("phospho_ID", "Hybridization", "chr"), seed=NULL) {

  ylength=length(yList)
  xlength=length(xList)
  x.prob=setdiff(1:xlength, NoProbXIndex)
  xlength2=length(x.prob)

  eFDR.grid=fdr_cutPob= No.Identified.filter=No.Identified.no.filter=
    PostProb=eFDR.no.filter=sig.no.filter=sig.filter=xNames=vector("list", xlength2);

  for (q in 1:ylength) {
    print(c("Outcome", q))


    if (is.null(seed)){seed.qq=NULL}
    if (is.null(seed)==F) {seed.qq =as.integer((seed+q)*4109)}

    eFDR1=iProFun.eFDR.1y(reg.all=reg.all, which.y=q, yList=yList, xList=xList, covariates=covariates, pi1=pi1,
                          NoProbXIndex=NoProbXIndex,filter=filter,
                          permutate_number=permutate_number,
                          grids = grids,fdr =fdr, PostCut=PostCut,
                          var.ID=var.ID, var.ID.additional=var.ID.additional, seed=seed.qq)
   for (p in 1:xlength2 ) {
     eFDR.grid[[p]]= cbind(eFDR.grid[[p]], eFDR1$eFDR.grid[[p]]);
     fdr_cutPob[[p]]= cbind(fdr_cutPob[[p]], eFDR1$fdr_cutPob[[p]]);
     No.Identified.filter[[p]]= cbind(No.Identified.filter[[p]], eFDR1$No.Identified.filter[[p]]);
     No.Identified.no.filter[[p]]= cbind(No.Identified.no.filter[[p]], eFDR1$No.Identified.no.filter[[p]]);

     PostProb[[p]]= cbind(PostProb[[p]], eFDR1$Gene_fdr[[p]]$PostProb);
     eFDR.no.filter[[p]]= cbind(eFDR.no.filter[[p]], eFDR1$Gene_fdr[[p]]$eFDR.no.filter);
     sig.no.filter[[p]]= cbind(sig.no.filter[[p]], eFDR1$Gene_fdr[[p]]$sig.no.filter);
     sig.filter[[p]]= cbind(sig.filter[[p]], eFDR1$Gene_fdr[[p]]$sig.filter);

     rownames(eFDR.grid[[p]])=grids
     xNames[[p]]=eFDR1$Gene_fdr[[p]][,setdiff(colnames(eFDR1$Gene_fdr[[p]]), c("PostProb", "eFDR.no.filter", "sig.no.filter", "sig.filter"))]
   }

  }

  eFDR_result=list(eFDR.grid=eFDR.grid, fdr_cutPob=fdr_cutPob,
                   No.Identified.no.filter=No.Identified.no.filter,
                   No.Identified.filter=No.Identified.filter,
                   PostProb=PostProb,
                   eFDR.no.filter=eFDR.no.filter,
                   sig.no.filter=sig.no.filter,
                   sig.filter=sig.filter, xNames=xNames)
  return(eFDR_result)

}


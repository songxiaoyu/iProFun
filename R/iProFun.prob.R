
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
iProFun.sum.prob.1y <- function (prob, which.y=1, NoProbXIndex=NULL) {
  xlength=length(prob)
  x.prob=setdiff(1:xlength, NoProbXIndex)
  
  # calculate prob on 1y
  prob.all <- vector("list", length(x.prob))
  
  for (p in 1:length(x.prob)) {
    k= x.prob[p]
    Q_index=prob[[k]]$Config[,which.y]
    prob.all[[p]]=apply(prob[[k]]$PostProb[,which(Q_index==1)], 1, sum)
    
    config.miss=prob[[k]]$Config.miss
    if (is.null(config.miss)==F) {
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
        
      } 
    }# end (j), which fills the list of missing values
  }
  return(prob.all)
}


#' iProFun Integrative analysis
#'
#' @param ylist ylist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an outcome. The matrix is formated as each outcome variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable ID for integrative iProFun analysis. If more than one rows
#' are availble for one variable ID (e.g. multiple phosphosites per gene with gene ID as
#' variable ID), only the row that are most promising is selected via ANOVA for the final
#' output. Example of ylist is a list of mRNA, global protein and phosphoprotein.
#' @param xlist xlist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an predictor. All data types or platforms are simutaneouosly
#' considered in the same regression. The matrix is formated as each predictor variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable ID for integrative iProFun analysis. If more than one rows
#' are availble for one variable ID (e.g. multiple methylation sites per gene with gene ID as
#' variable ID), all rows are considered in the regression to account for the effects of each
#' other. Example of xlist is a list of CNA and DNA methylation.
#' @param covariates covariates is a list of data matrix. This list should have the same length
#' as ylist. For the regression on the ith outcome, the ith covariates matrix contains the variables
#' that to be adjusted in the regression. The matrix is formated as each covariate variable
#' as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable ID for integrative iProFun analysis. Example of covariates is
#' a list of principle components.
#' @param pi pi is pre-specified priori of proportion of non-null statistics in each set of regression.
#' iProFun is insensitive to the mis-specification of the priori within a reasonable range.
#' @param permutate whether to permuate certain data type/platform or not. permutate = 0 (default):
#' no permuatation, analysis on original data. 0 < permutate <= length of ylist: permuate the label of
#' the corresponding data type in ylist. For example, permutate =2, permute the y label of second data matrix (e.g.protein)
#' @param ID ID gives the variable name (e.g. gene/protein name) to match different data types and consider integrative
#' analysis on. If IDs are not specified, the first columns will be considered as ID variable. iProFun only
#' applies with a fair number of overlapping variables are overlappe across all data types.
#' @param x.ID x.ID specifies variable name for xlist. It provides the flexible option to use
#' different names for xlist and ylist.
#' @param y.ID y.ID specifies variable name for ylist. It provides the flexible option to use
#' different names for xlist and ylist.
#' @param sub.ID.common sub.ID.common specifies the common string for samples (e.g. TCGA). If not
#' specified, all variables except IDs and colum.to.keep will be searched for overlapping samples.
#' @param colum.to.keep colum.to.keep allows to output additional variables except IDs in the
#' final results. Often helpful if multiple rows (e.g. probes) are considered per ID variable
#' (e.g. gene) to allow clear specification of each result.
#' @param missing.rate.filter missing.rate.filter allows filtering out variables with high missing
#' rate in data.
#' @param verbose verbose=T print out the progress of iProFun analysis.
#'
#' @return list with the same length as xlist. Nested within each list, it contains
#' \item{betas_J:}{Coefficient estimate for predictor on each outcome}
#' \item{betas_se_J:}{Coefficent SE for predictor on each outcome}
#' \item{sigma2_J:}{Regrssion error term for predictor on each outcome}
#' \item{dfs_J:}{Regression degrees of freedom for predictor on each outcome}
#' \item{v_g_J:}{ (X^T X)^{-1} projection on predictor on each outcome}
#' \item{xName_J:}{Predictor name corresponds to each predictor outcome pair}
#' \item{yName_J:}{Outcome name corresponds to each predictor outcome pair}
#' \item{NoComputation:}{No of variables considered for this predictor with all outcomes}
#' \item{Config:}{Corresponding association patterns. Total number 2^J for J outcome data types}
#' \item{PostProb:}{Posterior probability for each predictor on each association pattern}
#' \item{colocProb:}{Averaged posterior probability for predictor on each association pattern}
#' \item{Tstat_L:}{T statistics for each predictor on each outcome}
#' \item{D0:}{Estimated density under the null for predictor on each outcome}
#' \item{D1:}{Estimated density under the alternative for predictor on each outcome}

#' @export iProFun
#' @importFrom magrittr "%>%"
#' @import tidyr
#' @importFrom tibble as.tibble
#' @import dplyr
#' @import purrr
#' @import metRology
#' @importFrom matrixStats rowMins
#' @examples
#' iprofun_result <- iProFun(ylist = list(rna, protein, phospho),
#' xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3))

iProFun <- function(ylist, xlist, covariates, pi, permutate = 0,
                    ID=NULL, x.ID=NULL, y.ID=NULL,  sub.ID.common="TCGA",
                    colum.to.keep=c("phospho_ID", "Hybridization", "chr"),
                    missing.rate.filter=NULL, verbose=T){
  # ----- Data cleaning and quality check ----- #
  # Check xlist, ylist and covariates data types
  ylength=length(ylist)
  xlength=length(xlist)
  zlength=length(covariates)
  if (class(ylist)!="list" | ylength<2) stop("ylist needs to be a list with at least 2 data types.")
  if (class(xlist)!="list") stop("xlist needs to be a list.")
  if ((class(covariates)!="list") | ylength!=zlength) stop("covariates need to provide a list of covariates with the same length as ylist.")

  # Create subject ID used for each regression
  if (is.null(x.ID)  & is.null(y.ID) & is.null(ID)) {SubIDExclude=unique(c(sapply(ylist, function(f) colnames(f)[1]), sapply(xlist, function(f) colnames(f)[1]), colum.to.keep))}
  if (is.null(ID)==F & is.null(y.ID)==T & is.null(x.ID)==T) {SubIDExclude=unique(c(ID, colum.to.keep))}
  if (is.null(ID)==T & is.null(x.ID)==F & is.null(y.ID)==F ) {SubIDExclude=unique(c(y.ID, x.ID, colum.to.keep))}

  if (is.null(sub.ID.common) ==F) { # e.g. TCGA
    if (class(sub.ID.common)!="character") stop("sub.ID.common needs to be a character")
    ySubID=lapply(ylist, function(f) colnames(f)[grepl(sub.ID.common, colnames(f))])
    xSubID=lapply(xlist, function(f) colnames(f)[grepl(sub.ID.common, colnames(f))])
    zSubID=lapply(covariates,  function(f) colnames(f)[grepl(sub.ID.common, colnames(f))])
  } else if (is.null(sub.ID.common)) {
    if (verbose==T) {print("All columns that are not listed in xlist and ylist as variable
                           IDs and colum.to.keep are considered as samples.")}
    ySubID=lapply(ylist, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
    xSubID=lapply(xlist, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
    zSubID=lapply(covariates, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
  }
  xCommonSubID=Reduce(intersect, xSubID)
  yzCommonSubID=lapply(1:ylength, function(f) Reduce(intersect, list(ySubID[[f]], zSubID[[f]])))
  xyzCommonSubID=lapply(yzCommonSubID, function(f) Reduce(intersect, list(f, xCommonSubID)))


  # Allow missing rate filtering
  if (is.null(missing.rate.filter)==F) {
    if (class(missing.rate.filter)!="numeric" | missing.rate.filter>=1 | missing.rate.filter<=0) stop("missing.rate.filter needs to be NULL or numeric (0,1)")
    for (i in 1:(xlength) ){
      xlist[[i]]=xlist[[i]][apply(xlist[[i]], 1, function(f) mean(is.na(f))<=missing.rate.filter),]
    }
    for (i in 1:(ylength) ){
      ylist[[i]]=ylist[[i]][apply(ylist[[i]], 1, function(f) mean(is.na(f))<=missing.rate.filter),]
    }
  }

  # Create common Gene ID
  if (is.null(x.ID)  & is.null(y.ID) & is.null(ID)) {
    if (verbose==T) {print("No IDs are specified. The first column of each data type in ylist and xlist are used.")}

    Gene_ID_y=lapply(ylist, function(f) f[,1])
    Gene_ID_x=lapply(xlist, function(f) f[,1])
  }
  if (is.null(ID)==F & is.null(y.ID)==T & is.null(x.ID)==T) {
    Gene_ID_y=lapply(ylist, function(f) f[ID][,1])
    Gene_ID_x=lapply(xlist, function(f) f[ID][,1])
  }
  if (is.null(ID)==T & is.null(x.ID)==F & is.null(y.ID)==F ) {
    Gene_ID_y=lapply(ylist, function(f) f[y.ID][,1])
    Gene_ID_x=lapply(xlist, function(f) f[x.ID][,1])
  }
  if (is.null(x.ID)==F & is.null(y.ID)==F & is.null(ID)==F) stop("Only ID or x.ID/y.ID is needed to be specified")
  if (any(is.null(x.ID) & is.null(y.ID)==F,is.null(x.ID)==F & is.null(y.ID) )) stop("Both x.ID and y.ID are needed to be specified")

  xyCommonGeneID=Reduce(intersect, append(Gene_ID_x, Gene_ID_y))

  # print parameters
  if (verbose==T) {print(paste0("A total of ", length(xyCommonGeneID), " genes/proteins/peptides/etc. are considered."))}
  if (length(xyCommonGeneID)<100) warning("Less than 100 genes/proteins/peptides/etc. are considered; the density estimation may not be robust.")
  if (verbose==T) {print(paste0("A total of ", sapply(xyzCommonSubID, length), " samples for Y", seq(1:ylength)," data type."))}


  # A code for permutation
  if (permutate > 0 & permutate<=ylength) {  xyzCommonSubID_permutate <- lapply(xyzCommonSubID, sample)}
  if (permutate<0 | permutate>ylength) stop("Permutation index can not be out of the range of ylist.")
  # colum.to.keep = ID or c(ID, "Phospho_ID") from X or Y in output



  # ----- Linear Regression ----- #

  betas_J <- vector("list", xlength);   betas_se_J <- vector("list", xlength);
  sigma2_J <- vector("list", xlength); dfs_J <- vector("list", xlength);
  v_g_J <- vector("list", xlength);   xName_J <- vector("list", xlength);
  yName_J <- vector("list", xlength)

  for ( j in 1: ylength) {
    if (verbose==T) {print(paste0("Obtaining regression summaries for data type ", j))}

    betas=vector("list", xlength); betas_se=vector("list", xlength);
    sigma2=vector("list", xlength); dfs=vector("list", xlength); v_g=vector("list", xlength);
    xName=vector("list", xlength);  yName=vector("list", xlength)

    for (i in 1:length(xyCommonGeneID)) {
    #for (i in 1:3) {
      # print(i)
      if (j == permutate){
        y=ylist[[j]][which(Gene_ID_y[[j]]==xyCommonGeneID[i]) , xyzCommonSubID_permutate[[j]]]
      } else {
        y=ylist[[j]][which(Gene_ID_y[[j]]==xyCommonGeneID[i]), xyzCommonSubID[[j]]]
      }

      y=t(y)
      x_i= lapply(1:xlength, function(f) xlist[[f]][which(Gene_ID_x[[f]]==xyCommonGeneID[i]), xyzCommonSubID[[j]] ])
      x_index=sapply(x_i, nrow)
      x=t(do.call(rbind,x_i))
      z= t(covariates[[j]][, xyzCommonSubID[[j]]])


      xx=as.matrix(cbind(1, x, z))
      zz = as.matrix(cbind(1, z))
      p_xx=ncol(xx)
      p_x=ncol(x)

      # ANOVA to select best Y if multiple Y exist for one variable ID
      if (ncol(y)>1) {
        y_complete = as.matrix(y)[complete.cases(xx), ,drop = F]
        zz_complete = zz[complete.cases(xx), ]
        xx_complete = xx[complete.cases(xx), ]
        yindex=which.min(sapply(1:ncol(y), function(f) anova(lm(y_complete[,f]~xx_complete-1), lm(y_complete[,f]~zz_complete-1))$"Pr(>F)"[2]))
      } else {yindex=1}
      yy=y[,yindex]


      ft=lm(yy~xx-1)
      alphas=as.matrix(ft$coef)
      df1=ft$df.residual
      s_g_square=sum(ft$residuals^2)/df1
      inverse_xx=vcov(ft)/s_g_square

      if(p_x==1) {C=as.matrix(c(0, 1, rep(0, p_xx-2)))}
      if (p_x>1) {
        C=matrix(0, nrow=p_xx, ncol=p_x)
        diag(C[-1,])[1:p_x]=1
      }

      betas_single=t(C) %*% alphas # variable of interest
      betas_se_single= t(C) %*% summary(ft)$coef[,2]
      sigma2_single=as.matrix(rep(s_g_square, p_x))
      dfs_single= as.matrix(rep(df1, p_x))
      v_g_single=as.matrix(diag(t(C) %*% inverse_xx %*% C))

      # annotation
         yName_single=ylist[[j]][which(Gene_ID_y[[j]]==xyCommonGeneID[i])[yindex],
                         Reduce(intersect, list(SubIDExclude, colnames(ylist[[j]])) )]

    # output as different list for different x platforms

    for (p in 1:xlength) {
      index=seq(sum(x_index[0:(p-1)])+1,length=x_index[p])

      betas[[p]]=c(betas[[p]], betas_single[index])
      betas_se[[p]]=c(betas_se[[p]], betas_se_single[index])
      sigma2[[p]]=c(sigma2[[p]], sigma2_single[index])
      dfs[[p]]=c(dfs[[p]], dfs_single[index])
      v_g[[p]]=c(v_g[[p]], v_g_single[index])

      xName_single=xlist[[p]][which(Gene_ID_x[[p]]==xyCommonGeneID[i]), Reduce(intersect, list(SubIDExclude, colnames(xlist[[p]])))]
      xName[[p]] <- rbind(xName[[p]], as.matrix(xName_single))
      yName[[p]] <- rbind(yName[[p]], do.call("rbind", replicate(x_index[p], yName_single, simplify = F)))
    }
  }

    for (p in 1:xlength) {
      betas_J[[p]]=cbind(betas_J[[p]], betas[[p]])
      betas_se_J[[p]]=cbind(betas_se_J[[p]], betas_se[[p]])
      sigma2_J[[p]]=cbind(sigma2_J[[p]], sigma2[[p]])
      dfs_J[[p]]=cbind(dfs_J[[p]], dfs[[p]])
      v_g_J[[p]]=cbind(v_g_J[[p]], v_g[[p]])
      xName_J[[p]] <- cbind(xName_J[[p]], as.matrix(xName[[p]]))
      yName_J[[p]] <- cbind(yName_J[[p]], as.matrix(yName[[p]]))
    }

}

  for (p in 1:xlength) {
    xName_J[[p]]=t(unique(t(xName_J[[p]])))
    yName_J[[p]]=t(unique(t(yName_J[[p]])))
  }


  # ------- iProFun with regression results  ------------------------------------------------------------------
  if (verbose==T) {print("Run iProFun")}
  Reg_output <- vector("list", xlength);
  x_iProFun <- vector("list", xlength);
  final_result <- vector(mode="list", xlength)
  for (p in 1:xlength) {
    # Summarize Regression
    Reg_output[[p]]= list(betas_J=betas_J[[p]], betas_se_J=betas_se_J[[p]],
                          sigma2_J=sigma2_J[[p]],dfs_J=dfs_J[[p]],
                          v_g_J = v_g_J[[p]], xName_J = xName_J[[p]],
                          yName_J = yName_J[[p]])

    # run iProFun
    x_iProFun[[p]]= MultiOmics_Input(Reg_output[[p]],pi1 = pi)

    # Summarize output
    names(final_result)[p] <- paste0("iProFun output for xlist ", p)
    final_result[[p]] <- vector(mode="list", length=6)
    final_result[[p]] = append(Reg_output[[p]], x_iProFun[[p]])
  }
  if (verbose==T) {print("Completed")}
  return(final_result)
}

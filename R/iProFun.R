#' iProFun Integrative analysis
#'
#' @param ylist ylist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an outcome. The matrix is formated as each outcome variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of ylist is
#' a list of mRNA, global protein and phosphoprotein.
#' @param xlist xlist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an predictor. All data types or platforms are simutaneouosly
#' considered in the same regression. The matrix is formated as each predictor variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of xlist is
#' a list of CNA and DNA methylation.
#' @param covariates covariates is a list of data matrix. This list should have the same length
#' as ylist. For the regression on the ith outcome, the ith covariates matrix contains the variables
#' that to be adjusted in the regression.
#' The matrix is formated as each covariate variable
#'as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of covariates is
#' a list of principle components.
#' @param pi pi is pre-specified priori of proportion of non-null statistics in each set of regression.
#' iProFun is insensitive to the mis-specification of the priori within a reasonable range.
#' @param permutate whether to permuate certain data type/platform or not. permutate = 0 (default): no permuatation, analysis
#' on original data. permutate > 0: permuate the label of outcome for the corresponding data type
#' in ylist. For example, 2 = permutate the y label of second data matrix (protein)
#'
#' @return list with 2 components
#' \item{Averaged posterior Probability:}{Averaged posterior probability for predictor on each association pattern}
#' \item{Posterior Probability:}{posterior probability for each pattern per predictor}
#' @export iProFun
#' @importFrom magrittr "%>%"
#' @import tidyr
#' @importFrom tibble as.tibble
#' @import dplyr
#' @import purrr
#' @import metRology
#' @importFrom matrixStats rowMins
#' @examples
#' iprofun_result <- iprofun(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3))

iProFun <- function(ylist, xlist, covariates, pi = rep(0.05, 3), permutate = 0,cl=NULL,
                    x.ID=NULL, y.ID=NULL, ID=NULL, sub.ID.common="TCGA", colum.to.keep=c("phospho_ID", "Hybridization", "chr"),
                    missing.rate.filter=NULL){
  
  # ----- Data cleaning and quality check ----- #
  # Check xlist, ylist and covariates data types
  ylength=length(ylist)
  xlength=length(xlist)
  zlength=length(covariates)
  if (class(ylist)!="list" | ylength<2) stop("ylist needs to be a list with at least 2 data types.")
  if (class(xlist)!="list") stop("xlist needs to be a list.")
  if ((class(covariates)!="list") | ylength!=zlength) stop("covariates need to provide a list of covariates with the same length as ylist.")
  
  # Create subject ID sub.ID.common = "TCGA" | NULL | any string
  if (is.null(x.ID)  & is.null(y.ID) & is.null(ID)) {SubIDExclude=unique(c(sapply(ylist, function(f) colnames(f)[1]), sapply(xlist, function(f) colnames(f)[1]), colum.to.keep))}
  if (is.null(ID)==F & is.null(y.ID)==T & is.null(x.ID)==T) {SubIDExclude=unique(c(ID, colum.to.keep))}
  if (is.null(ID)==T & is.null(x.ID)==F & is.null(y.ID)==F ) {SubIDExclude=unique(c(y.ID, x.ID, colum.to.keep))}

  if (is.null(sub.ID.common) ==F) { # e.g. TCGA
    if (class(sub.ID.common)!="character") stop("sub.ID.common needs to be a character")
    ySubID=lapply(ylist, function(f) colnames(f)[grepl(sub.ID.common, colnames(f))])
    xSubID=lapply(xlist, function(f) colnames(f)[grepl(sub.ID.common, colnames(f))])
    zSubID=lapply(covariates,  function(f) colnames(f)[grepl(sub.ID.common, colnames(f))])
  } else if (is.null(sub.ID.common)) {
    print("All columns that are not listed in xlist and ylist as variable IDs and colum.to.keep are considered as samples.")
    ySubID=lapply(ylist, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
    xSubID=lapply(xlist, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
    zSubID=lapply(covariates, function(f) Reduce(setdiff, SubIDExclude, colnames(f)))
  }
  xCommonSubID=Reduce(intersect, xSubID)
  yzCommonSubID=lapply(1:3, function(f) Reduce(intersect, list(ySubID[[f]], zSubID[[f]])))
  xyzCommonSubID=lapply(yzCommonSubID, function(f) Reduce(intersect, list(f, xCommonSubID)))
  

  # Allow missing rate filtering
  if (is.null(missing.rate.filter)==F) {
    if (class(missing.rate.filter)!="numeric" | missing.rate.filter>=1 | missing.rate.filter<=0) stop("missing.rate.filter needs to be NULL or numeric (0,1)")
    for (i in 1:(xlength) ){
      xlist[[i]]=xlist[[i]][parApply(cl, xlist[[i]], 1, function(f) mean(is.na(f))<=missing.rate.filter),]  
    }
    for (i in 1:(ylength) ){
      ylist[[i]]=ylist[[i]][parApply(cl, ylist[[i]], 1, function(f) mean(is.na(f))<=missing.rate.filter),]  
    }
  }

  # Create Gene ID
  if (is.null(x.ID)  & is.null(y.ID) & is.null(ID)) {
    print("No IDs are specified. The first column of each data type in ylist and xlist are used as gene/protein/peptide/etc. ID.")

    Gene_ID_y=lapply(ylist, function(f) f[,1])
    Gene_ID_x=lapply(xlist, function(f) f[,1])
  }
  if (is.null(ID)==F & is.null(y.ID)==T & is.null(x.ID)==T) {
    Gene_ID_y=lapply(ylist, function(f) f[ID])
    Gene_ID_x=lapply(xlist, function(f) f[ID])   
  }
  if (is.null(ID)==T & is.null(x.ID)==F & is.null(y.ID)==F ) {
    Gene_ID_y=lapply(ylist, function(f) f[y.ID])
    Gene_ID_x=lapply(xlist, function(f) f[x.ID])    
  }
  if (is.null(x.ID)==F & is.null(y.ID)==F & is.null(ID)==F) stop("Only ID or x.ID/y.ID is needed to be specified")
  if (any(is.null(x.ID) & is.null(y.ID)==F,is.null(x.ID)==F & is.null(y.ID) )) stop("Both x.ID and y.ID are needed to be specified")

  xyCommonGeneID=Reduce(intersect, append(Gene_ID_x, Gene_ID_y))
    
  # print parameters
  print(paste0("A total of ", length(xyCommonGeneID), " genes/proteins/peptides/etc. are considered for iProFun integrative analysis."))
  if (length(xyCommonGeneID)<100) warning("Less than 100 genes/proteins/peptides/etc. are considered; the density estimation may not be robust.")
  print(paste0("A total of ", sapply(xyCommonSubID, length), " samples are considered for Y", seq(1:ylength)," data type."))
  

  # A code for permutation
  xyCommonSubID_permutate <- lapply(xyCommonSubID, sample)
  # colum.to.keep = ID or c(ID, "Phospho_ID") from X or Y in output

  
  
  # ----- Linear Regression ----- #
  # estimate
  betas_J <- vector("list", xlength);   betas_se_J <- vector("list", xlength);   sigma2_J <- vector("list", xlength);
  dfs_J <- vector("list", xlength);   v_g_J <- vector("list", xlength);   xName_J <- vector("list", xlength); yName_J <- vector("list", xlength)

  for ( j in 1: ylength) {
    # estimate_j
    betas=NULL; betas_se=NULL; sigma2=NULL; dfs=NULL; v_g=NULL;betas_2=NULL; betas_se_2=NULL;v_g_2=NULL; sigma2_2=NULL; dfs_2=NULL
    
    for (i in 1:length(xyCommonGeneID)) {
      # print(i)
      # more than one x or y from a gene is possible: select the most significant y; run every x.

      if (j == permutate){
        # ylist[[j]]$Gene_ID does not exist, create Gene_ID
        y=ylist[[j]][which(Gene_ID_y[[j]]==xyCommonGeneID[i]) , xyCommonSubID_permutate[[j]]]
      } else {
        y=ylist[[j]][which(Gene_ID_y[[j]]==xyCommonGeneID[i]), xyCommonSubID[[j]]]
      }

      y=data.frame(t(y))
      x_i= lapply(1:xlength, function(f) xlist[[f]][which(Gene_ID_x[[f]]==xyCommonGeneID[i]), xyCommonSubID[[j]] ])
      x_index=sapply(x_i, nrow)
      x=t(do.call(rbind,x_i))
      z= t(covariates[[j]][, xyCommonSubID[[j]]])


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

      if(p_x==1) {C=as.matrix(c(0, 1, rep(0, p_xx-2)))}
      if (p_x>1) {
        C=matrix(0, nrow=p_xx, ncol=p_x)
        diag(C[-1,])[1:p_x]=1
      }

      betas=rbind(betas, t(C) %*% alphas) # variable of interest
      betas_se=rbind(betas_se, t(C) %*% summary(ft)$coef[,2])
      inverse_xx=vcov(ft)/s_g_square
      sigma2=rbind(sigma2, as.matrix(rep(s_g_square, p_x)))
      dfs=rbind(dfs, as.matrix(rep(df1, p_x)))
      v_g=rbind(v_g, as.matrix(diag(t(C) %*% inverse_xx %*% C)))

      # annotation - NEED TO FIX ANNOTATION WITH LOOP
        xName=lapply(1:xlength, function(f) xlist[[f]][which(Gene_ID_x[[f]]==xyCommonGeneID[i]), Reduce(intersect, list(SubIDExclude, colnames(xlist[[f]])))])
        yName=ylist[[j]][which(Gene_ID_y[[j]]==xyCommonGeneID[i])[yindex], 
                         Reduce(intersect, list(SubIDExclude, colnames(ylist[[j]])) )]
     
  
    }
    
    # output as different list for different x platforms
    for (i in 1:xlength) {
      betas_J[[i]]=cbind(betas_J[[i]], betas[x_index[i]])
      betas_se_J[[i]]=cbind(betas_se_J[[i]], betas_se[x_index[i]])
      sigma2_J[[i]]=cbind(sigma2_J[[i]], sigma2[x_index[i]])      
      dfs_J[[i]]=cbind(dfs_J[[i]], dfs[x_index[i]]) 
      v_g_J[[i]]=cbind(v_g_J[[i]], v_g[x_index[i]])  
      #xName_J <- cbind(xName_J[[i]], xName[[i]])  
      #yName_J <- cbind(yName_J[[i]], yName) 
    }

  }
  # STOPPED HERE TO GO HOME
  x_1_list <- list(betas_J=betas_J, betas_se_J=betas_se_J, sigma2_J=sigma2_J,dfs_J=dfs_J, v_g_J = v_g_J, xName = xName, yName = yName)
  x_2_list <- list(betas_J=betas_J_2, betas_se_J=betas_se_J_2, sigma2_J=sigma2_J_2,dfs_J=dfs_J_2, v_g_J = v_g_J_2, xName = xName_2, yName = yName_2)

  
  
  # ------- Primo with regression results  ------------------------------------------------------------------

  for (i in 1:xlength) {
    assign(paste0("x_", i, "_iProFun"), MultiOmics_Input(eval(parse(text = paste0("x_", i, "_list"))),pi1 = pi) )
  }

  Q = makeQ(1:ylength)
  group <- NULL
  for (i in 1: dim(Q)[1]){
    group = c(group, paste(Q[i,], collapse = ""))
  }
  # Output lists with the following:

  # group
  for (i in 1:xlength) {
    assign(paste0("x_", i, "_iProFun_marginal_probability"), eval(parse(text = paste0("x_",i, "_iProFun")))$colocProb)
  }
  # x_1_iProFun_marginal_probability
  # x_2_iProFun_marginal_probability

  for (i in 1:xlength) {
    assign(paste0("x_", i, "_iProFun_gene_posterior_probability"), as.data.frame(eval(parse(text = paste0("x_",i, "_iProFun")))$PostProb))
  }
  # x_1_iProFun_gene_posterior_probability %>% head

  for (i in 1:xlength) {
    assign(paste0("x_", i, "_iProFun_gene_beta"), cbind(as.data.frame(eval(parse(text = paste0("x_",i, "_list")))$xName), eval(parse(text = paste0("x_",i, "_list")))$betas_J))
  }
  # x_1_iProFun_gene_beta
  # x_2_iProFun_gene_beta

  final_result <- vector(mode="list", length=6*xlength + 1)
  names(final_result)[1] <- "group"

  for (i in 1:xlength) {
    names(final_result)[i+1] = paste0("x_",i, "_iProFun_marginal_probability")
    names(final_result)[(i+xlength+1)] <- paste0("x_",i, "_iProFun_gene_posterior_probability")
    names(final_result)[(i+2*xlength+1)] <- paste0("x_",i, "_iProFun_gene_beta")
    names(final_result)[(i+3*xlength+1)] <- paste0("x_",i, "_iProFun_gene_posterior_probability_only")
    names(final_result)[(i+4*xlength+1)] <- paste0("x_",i, "_xName")
    names(final_result)[(i+5*xlength+1)] <- paste0("x_",i, "_yName")
  }
  final_result[[1]] <- group
  for (i in 1:xlength) {
    final_result[[i+1]] = eval(parse(text = paste0("x_",i, "_iProFun_marginal_probability")))
    final_result[[(i+xlength+1)]] <- eval(parse(text = paste0("x_",i, "_iProFun_gene_posterior_probability")))
    final_result[[(i+2*xlength+1)]] <- eval(parse(text = paste0("x_",i, "_iProFun_gene_beta")))
    final_result[[(i+3*xlength+1)]] <- eval(parse(text = paste0("x_",i, "_iProFun")))$raw_postprob
    final_result[[(i+4*xlength+1)]] <- eval(parse(text = paste0("x_",i, "_iProFun")))$xname
    final_result[[(i+5*xlength+1)]] <- eval(parse(text = paste0("x_",i, "_iProFun")))$yName
  }
  return(final_result)
}

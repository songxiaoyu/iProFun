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

iProFun <- function(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3), permutate = 0){
  #ylength
  #xlength
  #zlength

  #stopif zlength!=ylength ("Add error message")
  #class(xlist)!="list" ('x ...')
  #the first column of x,y,zlist to the ID to uniquly identify samples.
  # ID no overlap; error message
  # ID few overlap; warning message
  # provide note: iProFun is running on X samples, X predictors that overlaps for ylist[[i]]
  # covariate name: rna_cov
  # multiple samples in a ID: treatment to keep other variables as addtional ID.


  ylength=length(ylist)
  xlength=length(xlist)
  # How to specify the common Gene ID in general?
  xyCommonGeneID <- ylist[[1]]$Gene_ID
  # "TCGA" need to be generalized
  xyCommonSubID <- lapply(1:ylength, function(i) names(ylist[[i]])[grepl("TCGA", names(ylist[[i]]))])
  xyCommonSubID_permutate <- lapply(xyCommonSubID, sample)

  # estimate
  betas_J=NULL; betas_se_J=NULL; sigma2_J=NULL; dfs_J=NULL; v_g_J=NULL; betas_J_2=NULL; betas_se_J_2=NULL; v_g_J_2=NULL; yName=NULL; xName=NULL; xName_2 = NULL; sigma2_J_2=NULL; dfs_J_2=NULL;yName_2=NULL;

  for ( j in 1: ylength) {
    # estimate_j
    betas=NULL; betas_se=NULL; sigma2=NULL; dfs=NULL; v_g=NULL;betas_2=NULL; betas_se_2=NULL;v_g_2=NULL; sigma2_2=NULL; dfs_2=NULL
    for (i in 1:length(xyCommonGeneID)) {
      # print(i)
      # more than one x or y from a gene is possible: select the most significant y; run every x.

      if (j == permutate){
        y=ylist[[j]][ ylist[[j]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID_permutate[[j]]]
      } else {
        y=ylist[[j]][ ylist[[j]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID[[j]]]
      }

      y=data.frame(t(y))
      x= xlist[[1]][xlist[[1]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID[[j]]]
      x=data.frame(t(x))

      z= xlist[[-1]][xlist[[-1]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID[[j]]]
      covariates_model = covariates[[j]][, xyCommonSubID[[j]]]
      z = cbind(t(z), t(covariates_model))

      xx=as.matrix(cbind(1, x, z))
      zz = as.matrix(cbind(1, t(covariates_model)))
      p_xx=ncol(xx)
      p_x=ncol(x)
      inverse_xx<-my.solve(t(xx) %*% xx)
      n=nrow(y)

      # select y
      yindex=which.min(sapply(1:ncol(y), function(f) anova(lm(y[,f]~xx-1), lm(y[,f]~zz-1))$"Pr(>F)"[2]))
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

      sigma2=rbind(sigma2, as.matrix(rep(s_g_square, p_x)))
      dfs=rbind(dfs, as.matrix(rep(df1, p_x)))
      v_g=rbind(v_g, as.matrix(diag(t(C) %*% inverse_xx %*% C)))

      ### for the other x (make the assumption that there are only 2 x)
      x_2 = xlist[[2]][xlist[[2]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID[[j]]]
      p_x_2 <- nrow(x_2)
      if(p_x_2==1) {C_2=as.matrix(c(0, rep(0,p_x), 1, rep(0, p_xx - p_x -2)))}
      if (p_x_2>1) {
        C_2=matrix(0, nrow=p_xx, ncol=p_x_2)
        diag(C_2[-(1:(p_x+1)),])[1:p_x_2]=1
      }
      betas_2=rbind(betas_2, t(C_2) %*% alphas) # variable of interest
      betas_se_2=rbind(betas_se_2, t(C_2) %*% summary(ft)$coef[,2])
      # sigma2, dfs are the same as the first x
      v_g_2=rbind(v_g_2, as.matrix(diag(t(C_2) %*% inverse_xx %*% C_2)))
      dfs_2=rbind(dfs_2, as.matrix(rep(df1, p_x_2)))
      sigma2_2=rbind(sigma2_2, as.matrix(rep(s_g_square, p_x_2)))
      # annotation

      if (j==3) {
        xName=rbind(xName, xlist[[1]][xlist[[1]]$Gene_ID==xyCommonGeneID[i],
                                      colnames(xlist[[1]])[!grepl("TCGA", colnames(xlist[[1]]))]])
        yName=rbind(yName, matrix(rep(ylist[[j]][ylist[[j]]$Gene_ID==xyCommonGeneID[i],"phospho_ID"][yindex], p_x)))
        xName_2=rbind(xName_2, xlist[[2]][xlist[[2]]$Gene_ID==xyCommonGeneID[i],
                                          colnames(xlist[[2]])[!grepl("TCGA", colnames(xlist[[2]]))]])
        yName_2=rbind(yName_2, matrix(rep(ylist[[j]][ylist[[j]]$Gene_ID==xyCommonGeneID[i],"phospho_ID"][yindex], p_x_2)))
      }
      #
    }

    betas_J=cbind(betas_J, betas)
    betas_se_J=cbind(betas_se_J, betas_se)
    betas_J_2=cbind(betas_J_2, betas_2)
    betas_se_J_2=cbind(betas_se_J_2, betas_se_2)
    sigma2_J=cbind(sigma2_J, sigma2);
    sigma2_J_2=cbind(sigma2_J_2, sigma2_2);
    dfs_J=cbind(dfs_J, dfs);
    dfs_J_2=cbind(dfs_J_2, dfs_2);
    v_g_J=cbind(v_g_J, v_g)
    v_g_J_2=cbind(v_g_J_2, v_g_2)
  }
  x_1_list <- list(betas_J=betas_J, betas_se_J=betas_se_J, sigma2_J=sigma2_J,dfs_J=dfs_J, v_g_J = v_g_J, xName = xName, yName = yName)
  x_2_list <- list(betas_J=betas_J_2, betas_se_J=betas_se_J_2, sigma2_J=sigma2_J_2,dfs_J=dfs_J_2, v_g_J = v_g_J_2, xName = xName_2, yName = yName_2)
  # Primo (The following names need to be generalized) ------------------------------------------------------------------
  cnv_iprofun <- MultiOmics_Input(x_1_list,pi1 = pi) # You don't know x is cnv and y is methy
  methy_iprofun <- MultiOmics_Input(x_2_list,pi1 = pi)
  iprofun_result <- list("Marginal Probability" = tibble(group =
                                                           c(
                                                             "None",
                                                             "RNA only",
                                                             "Protein only",
                                                             "Phosphosite only",
                                                             "RNA & Protein",
                                                             "RNA & Phospho",
                                                             "Protein & Phospho",
                                                             "All three"
                                                           ), "CNA" = cnv_iprofun$colocProb,"Methy" = methy_iprofun$colocProb),
                         "Gene Posterior Probability" = methy_iprofun$PostProb %>%
                           as.tibble() %>%
                           mutate(X = "Methylation") %>%
                           select(X, everything()) %>%
                           # magrittr::set_names(c("Gene_ID", "None","RNA only","Protein only","Phosphosite only","RNA & Protein" ,"RNA & Phospho" ,"Protein & Phospho", "All three"))
                           rename("None" = `1` ,
                                  "RNA only" = `2`,
                                  "Protein only" = `3` ,
                                  "Phosphosite only" = `4` ,
                                  "RNA & Protein" = `5` ,
                                  "RNA & Phospho" = `6`,
                                  "Protein & Phospho" = `7`,
                                  "All three" = `8`) %>%
                           bind_rows(cnv_iprofun$PostProb %>%
                                       as.tibble() %>%
                                       mutate_at(vars(V2:V9), as.numeric) %>%
                                       magrittr::set_names(c("Gene_ID", "None","RNA only","Protein only","Phosphosite only","RNA & Protein" ,"RNA & Phospho" ,"Protein & Phospho", "All three")) %>%
                                       mutate(X = "CNA") %>%
                                       # rename("None" = `V1` ,
                                       #        "RNA only" = `V2`,
                                       #        "Protein only" = `V3` ,
                                       #        "Phosphosite only" = `V4` ,
                                       #        "RNA & Protein" = `V5` ,
                                       #        "RNA & Phospho" = `V6`,
                                       #        "Protein & Phospho" = `V7`,
                                       #        "All three" = `V8`) %>%
                                       mutate(Hybridization = NA_character_, chr = NA_character_) ),
                         "Beta" = cbind(x_1_list$xName, x_1_list$betas_J) %>% as.tibble() %>%mutate_at(vars(V2:V4), as.numeric) %>% mutate(Hybridization = NA_character_, chr = NA_character_, X = "CNA") %>% rename(Gene_ID = V1, `1` =V2,  `2`  =V3, `3` = V4) %>% bind_rows(
                           cbind(x_2_list$xName, x_2_list$betas_J) %>% as.tibble() %>%
                             mutate(X = "Methylation")))
  return(iprofun_result)
}

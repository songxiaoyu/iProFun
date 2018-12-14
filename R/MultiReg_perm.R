# ------------------------------------------------------------------- #
## e.g. Multiple Regression
# ------------------------------------------------------------------- #
# ylist <- list(rna, protein, phospho)
# xlist <- list(cnv_lr, cnv_mzd, methy_mean)
# covariates <- purity_tumor


MultiReg_together_perm <- function(ylist, xlist, covariates = NULL, conditional_covariate = NULL, cl=NULL, mutation_genes = NULL, filename, permcolum){
  
  ylength=length(ylist)
  xlength=length(xlist)
  # get overlapping samples by y
  ySubID = lapply(ylist, function(x) colnames(x)[grepl("CPT", colnames(x))])
  xSubID = lapply(xlist, function(x) colnames(x)[grepl("CPT", colnames(x))])
  covSubID <- colnames(covariates)[grepl("CPT", colnames(covariates))]
  xCommonSubID = intersect(Reduce(intersect, xSubID), covSubID)
  xyCommonSubID = lapply(ySubID, function(x) Reduce(intersect, list(x, xCommonSubID)))
  allxCommomSubID = unique(unlist(xyCommonSubID))
  xlist=lapply(1:xlength, function(i) xlist[[i]][,c(colnames(xlist[[i]])[!grepl("CPT", colnames(xlist[[i]]))], 
                                                 allxCommomSubID), with=F])
  ylist=lapply(1:ylength, function(i) ylist[[i]][,c(colnames(ylist[[i]])[!grepl("CPT", colnames(ylist[[i]]))],xyCommonSubID[[i]]), with=F]) 

  
  # control genes with missing values 
  xlist=lapply(xlist, function(i) i[complete.cases(i),])
  # ylist=lapply(ylist, function(i) i[complete.cases(i),])
  # ylist=lapply(ylist, function(i) i[parApply(cl, i, 1, function(f) mean(is.na(f))<=0.5),])  
  ylist[[2]]=ylist[[2]][parApply(cl, ylist[[2]], 1, function(f) mean(is.na(f))<=0.5),]  # mRNA first that we know don't have much missing
  ylist[[3]]=ylist[[3]][parApply(cl, ylist[[3]], 1, function(f) mean(is.na(f))<=0.5),]  # mRNA first that we know don't have much missing
  
  # get overlapping genes of all platforms
  yGeneID=lapply(ylist, function(x) data.frame(x)[,"Gene_ID"])
  xGeneID=lapply(xlist, function(x) data.frame(x)[,"Gene_ID"])
  xyCommonGeneID=Reduce(intersect, c(xGeneID, yGeneID))    # try with 1st 500 genes

  xlist=lapply(xlist, function(i) i[which(i$Gene_ID %in% xyCommonGeneID), ])
  ylist=lapply(ylist, function(i) i[which(i$Gene_ID %in% xyCommonGeneID), ])
  xyCommonSubID_permutate <- lapply(xyCommonSubID, sample)
  # not ordered 
 
  # estimate 
  betas_J = NULL
  betas_se_J = NULL
  sigma2_J = NULL
  dfs_J = NULL
  v_g_J = NULL
  r_square_J = NULL
  t_J = NULL
  yName = NULL
  xName = NULL
  
  for ( j in 1: ylength) {
    # estimate_j
    betas = NULL
    betas_se = NULL
    sigma2 = NULL
    dfs = NULL
    v_g = NULL
    r_square = NULL
    ts <- NULL
    for (i in 1:length(xyCommonGeneID)) {
      # print(i)
      # more than one x or y from a gene is possible: select the most significant y; run every x.
      
      if (j == permcolum){
        y=ylist[[j]][ ylist[[j]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID_permutate[[j]], with=F]
      } else {
        y=ylist[[j]][ ylist[[j]]$Gene_ID==xyCommonGeneID[i] , xyCommonSubID[[j]], with=F]
      }
      y = data.frame(t(y))
      x = xlist[[1]][xlist[[1]]$Gene_ID == xyCommonGeneID[i] , xyCommonSubID[[j]], with = F]
      x = data.frame(t(x))
      
      z1 = xlist[[2]][xlist[[2]]$Gene_ID == xyCommonGeneID[i] , xyCommonSubID[[j]], with = F]
      z2 = xlist[[3]][xlist[[3]]$Gene_ID == xyCommonGeneID[i] , xyCommonSubID[[j]], with = F]
      # z3 = xlist[[4]][xlist[[4]]$Gene_ID == xyCommonGeneID[i] , xyCommonSubID[[j]], with = F]
      covariates_model = covariates[, xyCommonSubID[[j]], with = F]
      if (xyCommonGeneID[i] %in% mutation_genes){
        conditional_covariate_model = conditional_covariate[conditional_covariate$Gene_ID == xyCommonGeneID[i], xyCommonSubID[[j]], with = F]
        z = cbind(t(z1), t(z2), t(covariates_model), t(conditional_covariate_model))
        xx = as.matrix(cbind(1, x, z))
        zz = as.matrix(cbind(1, t(covariates_model)))
        p_xx = ncol(xx)
        p_x = ncol(x)
        inverse_xx <- my.solve(t(xx) %*% xx)
        n = nrow(y)
        
      } else {
        z = cbind(t(z1), t(z2), t(covariates_model))
        xx = as.matrix(cbind(1, x, z))
        zz = as.matrix(cbind(1, t(covariates_model)))
        p_xx = ncol(xx)
        p_x = ncol(x)
        inverse_xx <- my.solve(t(xx) %*% xx)
        n = nrow(y)
      }
      
      # select y
      yindex=which.min(sapply(1:ncol(y), function(f) anova(lm(y[,f]~xx-1), lm(y[,f]~zz-1))$"Pr(>F)"[2]))
      yy = y[, yindex]
      
      ft = lm(yy ~ xx - 1)
      alphas = as.matrix(ft$coef)
      df1 = ft$df.residual
      s_g_square = sum(ft$residuals ^ 2) / df1
      if (p_x == 1) {
        C = as.matrix(c(0, 1, rep(0, p_xx - 2)))
      }
      if (p_x > 1) {
        C = matrix(0, nrow = p_xx, ncol = p_x)
        diag(C[-1, ])[1:p_x] = 1
      }
      
      betas = rbind(betas, t(C) %*% alphas) # variable of interest
      betas_se = rbind(betas_se, t(C) %*% summary(ft)$coef[, 2])
      ts <- rbind(ts, summary(ft)$coeff[2,3])
      sigma2 = rbind(sigma2, as.matrix(rep(s_g_square, p_x)))
      dfs = rbind(dfs, as.matrix(rep(df1, p_x)))
      v_g = rbind(v_g, as.matrix(diag(t(C) %*% inverse_xx %*% C)))
      r_square = rbind(r_square, summary(ft)$r.squared)
        # annotation
         
     if (j==3) {
       xName=rbind(xName, xlist[[1]][xlist[[1]]$Gene_ID==xyCommonGeneID[i],
                                     colnames(xlist[[1]])[!grepl("CPT", colnames(xlist[[1]]))], with=F]) 
       yName=rbind(yName, matrix(rep(ylist[[j]][ylist[[j]]$Gene_ID==xyCommonGeneID[i],Gene_ID][yindex], p_x)))
       name = rep("cnv", dim(xName)[1])
      }
#       
    }
    
    betas_J = cbind(betas_J, betas)
    betas_se_J = cbind(betas_se_J, betas_se)
    sigma2_J = cbind(sigma2_J, sigma2)
    dfs_J = cbind(dfs_J, dfs)
    v_g_J = cbind(v_g_J, v_g)
    r_square_J = cbind(r_square_J, r_square)
    t_J = cbind(t_J, ts)
  }

  return(list(filename=filename, betas_J=betas_J, betas_se_J=betas_se_J, sigma2_J=sigma2_J,dfs_J=dfs_J, v_g_J = v_g_J, r_square_J = r_square_J,t_J = t_J, xName = xName, yName = yName, name = name))
}







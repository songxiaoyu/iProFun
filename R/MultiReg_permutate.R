MultiReg_permutate <- function(){ set.seed(i)
  methy_perm_1_input=MultiReg_together_perm(ylist=ylist_normal, xlist=list(methy, cnv), filename="methy_perm_1",  covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),cl=cl, permcolum = 1)

  set.seed(i)
  methy_perm_2_input=MultiReg_together_perm(ylist=ylist_normal, xlist=list(methy, cnv), filename="methy_perm_2",  covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),cl=cl, permcolum = 2)


  set.seed(i)
  methy_perm_3_input=MultiReg_together_perm(ylist=ylist_normal, xlist=list(methy, cnv), filename="methy_perm_3",  covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),cl=cl, permcolum = 3)


  set.seed(i)
  cnv_perm_1_input=MultiReg_together_perm(ylist=ylist_normal, xlist=list(cnv, methy), filename="cnv_perm_1",   covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),cl=cl, permcolum = 1)


  set.seed(i)
  cnv_perm_2_input=MultiReg_together_perm(ylist=ylist_normal, xlist=list(cnv, methy), filename="cnv_perm_2",   covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),cl=cl, permcolum = 2)


  set.seed(i)
  cnv_perm_3_input=MultiReg_together_perm(ylist=ylist_normal, xlist=list(cnv, methy), filename="cnv_perm_3",   covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),cl=cl, permcolum = 3)


  # The following is for the overlapping data

  pi1=rep(0.05, 3)

  cnv_1_3 = MultiOmics_Input(cnv_input_1_3, pi1=pi1)
  methy_1_3 = MultiOmics_Input(methy_input_1_3, pi1=pi1)

  methy_perm_1 = MultiOmics_Input(methy_perm_1_input, pi1=pi1)
  cnv_perm_1 = MultiOmics_Input(cnv_perm_1_input, pi1=pi1)

  methy_perm_2 = MultiOmics_Input(methy_perm_2_input, pi1=pi1)
  cnv_perm_2 = MultiOmics_Input(cnv_perm_2_input, pi1=pi1)

  methy_perm_3 = MultiOmics_Input(methy_perm_3_input, pi1=pi1)
  cnv_perm_3 = MultiOmics_Input(cnv_perm_3_input, pi1=pi1)

  # FDR--------------------------------------------------------------------------------------------------------------------

  thresholds <- c(seq(0.05, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001))
  methy_perm_1_original <- vector("numeric", length(thresholds))
  cnv_perm_1_original <- vector("numeric", length(thresholds))
  methy_perm_2_original <- vector("numeric", length(thresholds))
  cnv_perm_2_original <- vector("numeric", length(thresholds))
  methy_perm_3_original <- vector("numeric", length(thresholds))
  cnv_perm_3_original <- vector("numeric", length(thresholds))
  methy_perm_2_count <- vector("numeric", length(thresholds))
  cnv_perm_2_count <- vector("numeric", length(thresholds))
  methy_perm_3_count <- vector("numeric", length(thresholds))
  cnv_perm_3_count <- vector("numeric", length(thresholds))
  methy_perm_1_count <- vector("numeric", length(thresholds))
  cnv_perm_1_count <- vector("numeric", length(thresholds))

  for (j in 1:length(thresholds)) {

    methy_perm_1_original[j] <- sum(methy_1_3$PostProb[methy_1_3$PostProb$`V2`+methy_1_3$PostProb$`V5`+methy_1_3$PostProb$`V6`+methy_1_3$PostProb$`V8` > thresholds[j],] %>% pull(2) %in% methy_filter)

    cnv_perm_1_original[j] <- sum(cnv_1_3$PostProb[cnv_1_3$PostProb$`V2`+cnv_1_3$PostProb$`V5`+cnv_1_3$PostProb$`V6`+cnv_1_3$PostProb$`V8` > thresholds[j],] %>% pull(1) %in% cnv_filter)

    methy_perm_2_original[j] <- sum(methy_1_3$PostProb[methy_1_3$PostProb$`V3`+methy_1_3$PostProb$`V5`+methy_1_3$PostProb$`V7`+methy_1_3$PostProb$`V8` > thresholds[j],] %>% pull(2) %in% methy_filter)

    cnv_perm_2_original[j] <- sum(cnv_1_3$PostProb[cnv_1_3$PostProb$`V3`+cnv_1_3$PostProb$`V5`+cnv_1_3$PostProb$`V7`+cnv_1_3$PostProb$`V8` > thresholds[j],] %>% pull(1) %in% cnv_filter)

    methy_perm_3_original[j] <- sum(methy_1_3$PostProb[methy_1_3$PostProb$`V4`+methy_1_3$PostProb$`V6`+methy_1_3$PostProb$`V7`+methy_1_3$PostProb$`V8` > thresholds[j],] %>% pull(2) %in% methy_filter)

    cnv_perm_3_original[j] <- sum(cnv_1_3$PostProb[cnv_1_3$PostProb$`V4`+cnv_1_3$PostProb$`V6`+cnv_1_3$PostProb$`V7`+cnv_1_3$PostProb$`V8` > thresholds[j],] %>% pull(1) %in% cnv_filter)

    methy_perm_2_count[j] <- sum(methy_perm_2$PostProb[methy_perm_2$PostProb$`V3`+methy_perm_2$PostProb$`V5`+methy_perm_2$PostProb$`V7`+methy_perm_2$PostProb$`V8`> thresholds[j],] %>% pull(2) %in% methy_filter)

    cnv_perm_2_count[j] <- sum(cnv_perm_2$PostProb[cnv_perm_2$PostProb$`V3`+cnv_perm_2$PostProb$`V5`+cnv_perm_2$PostProb$`V7`+cnv_perm_2$PostProb$`V8` > thresholds[j],] %>% pull(1) %in% cnv_filter)

    methy_perm_3_count[j] <- sum(methy_perm_3$PostProb[methy_perm_3$PostProb$`V4`+methy_perm_3$PostProb$`V6`+methy_perm_3$PostProb$`V7`+methy_perm_3$PostProb$`V8` > thresholds[j],] %>% pull(2) %in% methy_filter)

    cnv_perm_3_count[j] <- sum(cnv_perm_3$PostProb[cnv_perm_3$PostProb$`V4`+cnv_perm_3$PostProb$`V6`+cnv_perm_3$PostProb$`V7`+cnv_perm_3$PostProb$`V8` > thresholds[j],] %>% pull(1) %in% cnv_filter)

    methy_perm_1_count[j] <- sum(methy_perm_1$PostProb[methy_perm_1$PostProb$`V2`+methy_perm_1$PostProb$`V5`+methy_perm_1$PostProb$`V6`+methy_perm_1$PostProb$`V8` > thresholds[j],] %>% pull(2) %in% methy_filter)

    cnv_perm_1_count[j] <- sum(cnv_perm_1$PostProb[cnv_perm_1$PostProb$`V2`+cnv_perm_1$PostProb$`V5`+cnv_perm_1$PostProb$`V6`+cnv_perm_1$PostProb$`V8` > thresholds[j],] %>% pull(1) %in% cnv_filter)
  }

  MR1 <- tibble("Post Greater than" = thresholds,
                "FDR Permute RNA (CNV)" = cnv_perm_1_count/cnv_perm_1_original,
                "No. of Discovery (RNA only, Protein&RNA, RNA&Phospho, ALL three) (CNV)" = cnv_perm_1_original,
                "FDR Permute Protein (CNV)" = cnv_perm_2_count/cnv_perm_2_original,
                "No. of Discovery (Protein only, Protein&RNA, Protein&Phospho, ALL three) (CNV)" = cnv_perm_2_original,
                "FDR Permute Phospho (CNV)" = cnv_perm_3_count/cnv_perm_3_original,
                "No. of Discovery (Phospho only, Phospho&RNA, Protein&Phospho, ALL three) (CNV)" = cnv_perm_3_original,
                "FDR Permute RNA (Methy)" = methy_perm_1_count/methy_perm_1_original,
                "No. of Discovery (RNA only, Protein&RNA, RNA&Phospho, ALL three) (Methy)" = methy_perm_1_original,
                "FDR Permute Protein (Methy)" = methy_perm_2_count/methy_perm_2_original,
                "No. of Discovery (Protein only, Protein&RNA, Protein&Phospho, ALL three) (Methy)" = methy_perm_2_original,
                "FDR Permute Phospho (Methy)" = methy_perm_3_count/methy_perm_3_original,
                "No. of Discovery (Phospho only, Phospho&RNA, Protein&Phospho, ALL three) (Methy)" = methy_perm_3_original)


  result[[i]] <- MR1

  print(paste("Finish Seed", i))


}

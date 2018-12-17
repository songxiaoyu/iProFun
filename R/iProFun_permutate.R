#' Calculate the empirical FDR using permutation
#'
#' @param ylist mRNA, global protein and phosphoprotein
#' @param xlist CNA and methylation
#' @param covariates sets of covariates adjusted in the regression analyses
#' @param pi prior probability
#' @param permutate_number Number of permutation, default 10
#'
#' @return
#' @export
#'
#' @examples iprofun_permutate
iprofun_permutate = function(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3), permutate_number = 10){
  for (i in 1 : permutate_number) {
    set.seed(i)
    iprofun_perm_1 <- iProFun(permutate = 1)
    cnv_perm_1 <- iprofun_perm_1$`Gene Posterior Probability` %>%
      filter(X == "CNA")
    methy_perm_1 <- iprofun_perm_1$`Gene Posterior Probability` %>%
      filter(X == "Methylation")

    iprofun_perm_2 <- iProFun(permutate = 2)
    cnv_perm_2 <- iprofun_perm_2$`Gene Posterior Probability` %>%
      filter(X == "CNA")
    methy_perm_2 <- iprofun_perm_2$`Gene Posterior Probability` %>%
      filter(X == "Methylation")

    iprofun_perm_3 <- iProFun(permutate = 3)
    cnv_perm_3 <- iprofun_perm_3$`Gene Posterior Probability` %>%
      filter(X == "CNA")
    methy_perm_3 <- iprofun_perm_3$`Gene Posterior Probability` %>%
      filter(X == "Methylation")
    iprofun = iProFun()
    cnv_1_3 = iprofun$`Gene Posterior Probability` %>%
      filter(X == "CNA")
    methy_1_3 = iprofun$`Gene Posterior Probability` %>%
      filter(X == "Methylation")

    cnv_filter <- iprofun$Beta %>%
      dplyr::filter(
        X == "CNA",
        "estimate_cnv_rna" > 0,
        "estimate_cnv_protein" > 0,
        "estimate_cnv_phospho" > 0
      ) %>%
      pull(Gene_ID)

    methy_filter <- c(
      iprofun$Beta %>%
        dplyr::filter(
          X == "Methylation",
          "estimate_methy_rna" > 0,
          "estimate_methy_protein" > 0,
          "estimate_methy_phospho" > 0
        ) %>%
        pull(Hybridization),
      iprofun$Beta %>%
        dplyr::filter(
          X == "Methylation",
          "estimate_methy_rna" < 0,
          "estimate_methy_protein" < 0,
          "estimate_methy_phospho" < 0
        ) %>%
        pull(Hybridization)
    )

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

      methy_perm_1_original[j] <- sum(methy_1_3[methy_1_3$`RNA only`+methy_1_3$"RNA & Protein"+methy_1_3$"RNA & Phospho"+methy_1_3$"All three" > thresholds[j],] %>%pull(Hybridization) %in% methy_filter)

      cnv_perm_1_original[j] <- sum(cnv_1_3[cnv_1_3$`RNA only`+cnv_1_3$"RNA & Protein"+cnv_1_3$"RNA & Phospho"+cnv_1_3$"All three" > thresholds[j],] %>%pull(Gene_ID) %in% cnv_filter)

      methy_perm_2_original[j] <- sum(methy_1_3[methy_1_3$`Protein only`+methy_1_3$"RNA & Protein"+methy_1_3$"Protein & Phospho"+methy_1_3$"All three" > thresholds[j],] %>%pull(Hybridization) %in% methy_filter)

      cnv_perm_2_original[j] <- sum(cnv_1_3[cnv_1_3$`Protein only`+cnv_1_3$"RNA & Protein"+cnv_1_3$"Protein & Phospho"+cnv_1_3$"All three" > thresholds[j],] %>%pull(Gene_ID) %in% cnv_filter)

      methy_perm_3_original[j] <- sum(methy_1_3[methy_1_3$"Phosphosite only"+methy_1_3$"RNA & Phospho"+methy_1_3$"Protein & Phospho"+methy_1_3$"All three" > thresholds[j],] %>%pull(Hybridization) %in% methy_filter)

      cnv_perm_3_original[j] <- sum(cnv_1_3[cnv_1_3$"Phosphosite only"+cnv_1_3$"RNA & Phospho"+cnv_1_3$"Protein & Phospho"+cnv_1_3$"All three" > thresholds[j],] %>%pull(Gene_ID) %in% cnv_filter)

      methy_perm_2_count[j] <- sum(methy_perm_2[methy_perm_2$`Protein only`+methy_perm_2$"RNA & Protein"+methy_perm_2$"Protein & Phospho"+methy_perm_2$"All three"> thresholds[j],] %>%pull(Hybridization) %in% methy_filter)

      cnv_perm_2_count[j] <- sum(cnv_perm_2[cnv_perm_2$`Protein only`+cnv_perm_2$"RNA & Protein"+cnv_perm_2$"Protein & Phospho"+cnv_perm_2$"All three" > thresholds[j],] %>%pull(Gene_ID) %in% cnv_filter)

      methy_perm_3_count[j] <- sum(methy_perm_3[methy_perm_3$"Phosphosite only"+methy_perm_3$"RNA & Phospho"+methy_perm_3$"Protein & Phospho"+methy_perm_3$"All three" > thresholds[j],] %>%pull(Hybridization) %in% methy_filter)

      cnv_perm_3_count[j] <- sum(cnv_perm_3[cnv_perm_3$"Phosphosite only"+cnv_perm_3$"RNA & Phospho"+cnv_perm_3$"Protein & Phospho"+cnv_perm_3$"All three" > thresholds[j],] %>%pull(Gene_ID) %in% cnv_filter)

      methy_perm_1_count[j] <- sum(methy_perm_1[methy_perm_1$`RNA only`+methy_perm_1$"RNA & Protein"+methy_perm_1$"RNA & Phospho"+methy_perm_1$"All three" > thresholds[j],] %>%pull(Hybridization) %in% methy_filter)

      cnv_perm_1_count[j] <- sum(cnv_perm_1[cnv_perm_1$`RNA only`+cnv_perm_1$"RNA & Protein"+cnv_perm_1$"RNA & Phospho"+cnv_perm_1$"All three" > thresholds[j],] %>%pull(Gene_ID) %in% cnv_filter)
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
    # write.csv(MR1, paste0("result/Multiple Regression PC3_seed_", i, ".csv"))
    print(paste("Finish Seed", i))
  }
}

#' Calculate the empirical FDR using permutation
#'
#' @param ylist mRNA, global protein and phosphoprotein
#' @param xlist CNA and methylation
#' @param covariates sets of covariates adjusted in the regression analyses
#' @param pi prior probability
#' @param permutate_number Number of permutation, default 10
#' @param fdr False Discover Rate, default as 0.1
#' @param posterior Minimal Posterior Probabilty Cutoff, default as 0.75
#' @return A list with 3 sublists, the first one is the cutoff values for each greoup based on permutation, the second is the results for iProFun
#' @export
#'
#' @examples iprofun_permutate
iprofun_permutate = function(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3), permutate_number = 10, fdr = 0.1, posterior = 0.75){
  iprofun = iProFun()
  cnv_1_3 = iprofun$`Gene Posterior Probability` %>%
    filter(X == "CNA")
  methy_1_3 = iprofun$`Gene Posterior Probability` %>%
    filter(X == "Methylation")

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

    cnv_filter <- iprofun$Beta %>%
      dplyr::filter(
        X == "CNA",
        estimate_cnv_rna > 0,
        estimate_cnv_protein > 0,
        estimate_cnv_phospho > 0
      ) %>%
      pull(Gene_ID)

    methy_filter <- c(
      iprofun$Beta %>%
        dplyr::filter(
          X == "Methylation",
          estimate_methy_rna > 0,
          estimate_methy_protein > 0,
          estimate_methy_phospho > 0
        ) %>%
        pull(Hybridization),
      iprofun$Beta %>%
        dplyr::filter(
          X == "Methylation",
          estimate_methy_rna < 0,
          estimate_methy_protein < 0,
          estimate_methy_phospho < 0
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
    result <- vector("list", permutate_number)
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
  CCRCC_permutation_mean <- bind_rows(result, .id = "Seed") %>%
    group_by(`Post Greater than`) %>%
    dplyr::select(-c(Seed)) %>%
    summarise_all(funs(mean = mean))
  rna_cnv_cutoff <- CCRCC_permutation_mean %>%
    filter(`FDR Permute RNA (CNV)_mean` < fdr, `Post Greater than` > posterior) %>% pull(`Post Greater than`) %>% .[1]
  protein_cnv_cutoff <- CCRCC_permutation_mean %>%
    filter(`FDR Permute Protein (CNV)_mean` < fdr, `Post Greater than` > posterior) %>% pull(`Post Greater than`) %>% .[1]
  phospho_cnv_cutoff <- CCRCC_permutation_mean %>%
    filter(`FDR Permute Phospho (CNV)_mean` < fdr, `Post Greater than` > posterior) %>% pull(`Post Greater than`) %>% .[1]
  rna_methy_cutoff <- CCRCC_permutation_mean %>%
    filter(`FDR Permute RNA (Methy)_mean` < fdr, `Post Greater than` > posterior) %>% pull(`Post Greater than`) %>% .[1]
  protein_methy_cutoff <- CCRCC_permutation_mean %>%
    filter(`FDR Permute Protein (Methy)_mean` < fdr, `Post Greater than` > posterior) %>% pull(`Post Greater than`) %>% .[1]
  phospho_methy_cutoff <- CCRCC_permutation_mean %>%
    filter(`FDR Permute Phospho (Methy)_mean` < fdr, `Post Greater than` > posterior) %>% pull(`Post Greater than`) %>% .[1]

  gene_cnv_RNA <- cnv_filter[cnv_filter %in% (cnv_1_3[cnv_1_3$`RNA only`+cnv_1_3$`RNA & Protein`+cnv_1_3$`RNA & Phospho`+cnv_1_3$`All three`> rna_cnv_cutoff,] %>% pull(Gene_ID) )] %>% unique()


  gene_methy_RNA <-  as.tibble(methy_filter[methy_filter %in% (methy_1_3[methy_1_3$`RNA only`+methy_1_3$`RNA & Protein`+methy_1_3$`RNA & Phospho`+methy_1_3$`All three`> rna_methy_cutoff,] %>% pull(Hybridization))]) %>% inner_join(methy_1_3, by = c("value" = "Hybridization")) %>% pull(Gene_ID) %>% unique()

  hybrid_methy_RNA <- methy_filter[methy_filter %in% (methy_1_3[methy_1_3$`RNA only`+methy_1_3$`RNA & Protein`+methy_1_3$`RNA & Phospho`+methy_1_3$`All three`> rna_methy_cutoff,] %>% pull(Hybridization))]

  # protein

  gene_cnv_protein <- cnv_filter[cnv_filter %in% (cnv_1_3[cnv_1_3$`Protein only`+cnv_1_3$`RNA & Protein`+cnv_1_3$`Protein & Phospho`+cnv_1_3$`All three`> protein_cnv_cutoff,] %>% pull(Gene_ID))] %>% unique()

  gene_methy_protein <- as.tibble(methy_filter[methy_filter %in% (methy_1_3[methy_1_3$`Protein only`+methy_1_3$`RNA & Protein`+methy_1_3$`Protein & Phospho`+methy_1_3$`All three`> protein_methy_cutoff,] %>% pull(Hybridization))]) %>% inner_join(methy_1_3, by = c("value" = "Hybridization")) %>% pull(Gene_ID) %>% unique()

  hybrid_methy_protein <- as.tibble(methy_filter[methy_filter %in% (methy_1_3[methy_1_3$`Protein only`+methy_1_3$`RNA & Protein`+methy_1_3$`Protein & Phospho`+methy_1_3$`All three`> protein_methy_cutoff,] %>% pull(Hybridization))]) %>% inner_join(methy_1_3, by = c("value" = "Hybridization")) %>% pull(1) %>% unique()

  # phospho

  gene_cnv_phospho <- cnv_filter[cnv_filter %in% (cnv_1_3[cnv_1_3$`Phosphosite only`+cnv_1_3$`RNA & Phospho`+cnv_1_3$`Protein & Phospho`+cnv_1_3$`All three`> phospho_cnv_cutoff,] %>% pull(Gene_ID) )] %>% unique()


  gene_methy_phospho <- methy_1_3[methy_1_3$`Phosphosite only`+methy_1_3$`RNA & Phospho`+methy_1_3$`Protein & Phospho`+methy_1_3$`All three`> phospho_methy_cutoff,] %>% pull(Gene_ID) %>% unique()

  venn_table <- tibble("Gene" = gene_cnv_RNA, "CNV RNA" = 1) %>%
    full_join(tibble("Gene" = gene_methy_RNA, "Methy RNA" = 1), by = "Gene") %>%
    full_join(tibble("Gene" = gene_methy_protein, "Methy Protein" = 1), by = "Gene") %>%
    full_join(tibble("Gene" = gene_cnv_protein, "CNV Protein" = 1), by = "Gene") %>%
    full_join(tibble("Gene" = gene_methy_phospho, "Methy Phospho" = 1), by = "Gene") %>%
    full_join(tibble("Gene" = gene_cnv_phospho, "CNV Phospho" = 1), by = "Gene") %>%
    mutate_all(funs(replace(., is.na(.), 0)))

  iProFun_GeneList <- venn_table %>%
    bind_rows(tibble(Gene = setdiff(cnv_1_3 %>% pull(Gene_ID), venn_table %>% pull(Gene)),
                     "CNV RNA" = 0,
                     "Methy RNA" = 0,
                     "Methy Protein" = 0,
                     "CNV Protein" = 0,
                     "Methy Phospho" = 0,
                     "CNV Phospho" = 0))

  methy_all_negative <- iprofun$Beta %>%
    dplyr::filter(
      X == "Methylation",
      estimate_methy_rna < 0,
      estimate_methy_protein < 0,
      estimate_methy_phospho < 0
    ) %>% pull(Gene_ID) # same genes but different methy sites, so...

  iProFun_GeneList_negative <- iProFun_GeneList %>%
    mutate(`Methy RNA` = case_when(Gene %in% methy_all_negative ~ -`Methy RNA`,
           TRUE ~ `Methy RNA`),
           `Methy Protein` = case_when(Gene %in% methy_all_negative ~ -`Methy Protein`,
                                       TRUE ~ `Methy Protein`),
           `Methy Phospho` = case_when(Gene %in% methy_all_negative ~ -`Methy Phospho`,
                                       TRUE ~ `Methy Phospho`))
  cutoff <- tibble(
    Group = c(
      "CNV RNA",
      "CNV Protein",
      "CNV Phospho",
      "Methy RNA",
      "Methy Protein",
      "Methy Phospho"
    ),
    "Cut Off probability" = c(
      rna_cnv_cutoff,
      protein_cnv_cutoff,
      phospho_cnv_cutoff,
      rna_methy_cutoff,
      protein_methy_cutoff,
      phospho_methy_cutoff
    )
  )

  result_permutate <- list("Posterior Probability Cutoff" = cutoff,
                 "iProFun Result" = iProFun_GeneList,
                 "iProFun Result (Negative/Positive)" = iProFun_GeneList_negative)
  return(result_permutate)
}

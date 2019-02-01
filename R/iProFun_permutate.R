#' iProFun false discovery rate assessment based on permutations
#'
#' @param ylist ylist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an outcome. The matrix is formated as each outcome variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of ylist is
#' a list of mRNA, global protein and phosphoprotein.
#' @param xlist xlist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an predictor. The matrix is formated as each outcome variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of xlist is
#' a list of cna and methylation.
#' @param covariates covariates is a list of data matrix. Each matrix is a data type or platform that
#' one would like to adjust as a covariate. The matrix is formated as each outcome variable
#'as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of covariates is
#' a list of principle components.
#' @param pi prior probability
#' @param permutate_number Number of permutation, default 10
#' @param fdr False Discover Rate, default as 0.1
#' @param posterior Minimal Posterior Probabilty Cutoff, default as 0.75
#' @param filter is a vector with the same length of ylist, taking values of 1, -1, 0 or NULL. "1"
#' indicates filter in all positive results across multi-omic platforms. XXXX .... Jiayi!!!!
#' @param thresholds xxx Jiayi
#' @param seed Jiayi
#' @return list with 3 components
#' \item{Posterior Probability Cutoff:}{the cutoff values for each group based on permutation}
#' \item{iProFun Result:}{A table indicating whether a gene is identified by iProFun or not}
#' \item{iProFun Result (Negative/Positive):}{A table indicating whether a gene is identified by iProFun or not and whether the association is positive or not}
#' @export iProFun_permutate
#'
#' @examples
#' iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3), permutate_number = 1, fdr = 0.1, posterior = 0.75)
iProFun_permutate = function(ylist, xlist, covariates,
                             pi = rep(0.05, 3), permutate_number = 10, fdr = 0.1, posterior = 0.75, thresholds = c(seq(0.75, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001))){
                             # ,filter = NULL, seed=?.Random.seed

  iprofun_result = iProFun(ylist = ylist, xlist = xlist, covariates = covariates, pi=pi, permutate = 0 )

  # need to add filter here
  result <- vector("list", permutate_number)

  for (i in 1 : permutate_number) {
    # set.seed(seed)
    set.seed(i)
    for (j in 1:length(ylist)){
      assign(paste0("iprofun_perm_", j), iProFun(ylist = ylist, xlist = xlist, covariates = covariates, permutate = j))
    }



    for (p in 1:length(thresholds)) {

      for (j in 1:length(ylist)){
        for (k in 1:length(xlist)){
          assign(paste0("x_",k, "_iprofun_perm_", j, "_count"), vector("numeric", length(thresholds)))
          assign(paste0("x_",k, "_iprofun_perm_", j, "_original"), vector("numeric", length(thresholds)))
        }
      }

      for(j in 1:length(ylist)){
        for ( k in 1:length(xlist)) {
          assign(paste0("x_",k, "_iprofun_perm_", j, "_count"), dim(iprofun_result[[paste0("x_",k, "_xName")]][apply(eval(parse(text = paste0("iprofun_perm_", j)))[[paste0("x_",k, "_iProFun_gene_posterior_probability_only")]][, Q[,j]==1],1,sum) > thresholds[p], ,drop = F])[1])
          assign(paste0("x_",k, "_iprofun_perm_", j, "_original"), dim(iprofun_result[[paste0("x_",k, "_xName")]][apply(iprofun_result[[paste0("x_",k, "_iProFun_gene_posterior_probability_only")]][, Q[,j]==1],1,sum) > thresholds[p], , drop = F])[1])
        }
      }
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
    write.csv(MR1, paste0("Breast/Rdata/breast_seed_", i, ".csv"))
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

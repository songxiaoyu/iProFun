#' iProFun Integrative analysis pipeline
#'
#' @param ylist mRNA, global protein and phosphoprotein
#' @param xlist CNA and methylation
#' @param covariates sets of covariates adjusted in the regression analyses
#' @param pi prior probability
#' @param permutate whether to permuate or not. 1 = permuatte the first column(mRNA), 2 = permutate the second column, 3 = permutate the third column
#'
#' @return A list with 2 lists: The first sublist is the result for CNA, the second sublist is the result for Methylation
#' @export iProFun
#' @importFrom magrittr "%>%"
#' @import tidyr
#' @importFrom tibble as.tibble
#' @import dplyr
#' @import purrr
#' @import metRology
#' @import matrixStats
#' @examples
iProFun <- function(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi, permutate = 0){
  rna_regression <- ylist[[1]]
  protein_regression <- ylist[[1]]
  phospho_regression <- ylist[[3]]
  cnv_regression <- xlist[[1]]
  methy_regression <- xlist[[2]]
  rna_pc_1_3 <- covariates[[1]]
  protein_pc_1_3 <- covariates[[2]]
  phospho_pc_1_3 <- covariates[[3]]
  rna_regression_long <- rna_regression %>%
    gather(subject_id, rna, - Gene_ID) %>% as.tibble()

  protein_regression_long <- protein_regression %>%
    gather(subject_id, protein, - Gene_ID) %>% as.tibble()

  phospho_regression_long <- phospho_regression %>%
    gather(subject_id, phospho, - Gene_ID, - phospho_ID) %>% as.tibble()

  cnv_regression_long <- cnv_regression %>%
    gather(subject_id, cnv, - Gene_ID) %>% as.tibble()

  methy_regression_long <- methy_regression %>%
    as.tibble() %>%
    select(-"Hybridization", -"chr") %>%
    group_by(Gene_ID) %>%
    mutate(methy = row_number()) %>%
    gather(subject_id, value, - "Gene_ID",-methy) %>%
    ungroup() %>%
    mutate(methy = paste0("methy", methy)) %>%
    spread(methy, value, fill = 0)

  rna_name_permutate <- sample(names(rna_regression)[-1])
  rna_regression_permutated <- rna_regression
  names(rna_regression_permutated)[-1] <- rna_name_permutate
  rna_regression_long_permutated <- rna_regression_permutated %>%
    gather(subject_id, rna, - Gene_ID) %>% as.tibble()

  protein_name_permutate <- sample(names(protein_regression)[-1])
  protein_regression_permutated <- protein_regression
  names(protein_regression_permutated)[-1] <- protein_name_permutate
  protein_regression_long_permutated <- protein_regression_permutated %>%
    gather(subject_id, protein, - Gene_ID) %>% as.tibble()

  phospho_name_permutate <- sample(names(phospho_regression)[-c(1,2,3)])
  phospho_regression_permutated <- phospho_regression
  names(phospho_regression_permutated)[-c(1,2,3)] <- phospho_name_permutate
  phospho_regression_long_permutated <- phospho_regression_permutated %>%
    gather(subject_id, phospho, - Gene_ID, - phospho_ID)

  # rna ---------------------------------------------------------------------
if (permutate == 1){
  rna_regression_combine <- rna_regression_long_permutated %>%
    inner_join(cnv_regression_long) %>%
    inner_join(methy_regression_long) %>%
    inner_join(rna_pc_1_3)
} else {
  rna_regression_combine <- rna_regression_long %>%
    inner_join(cnv_regression_long) %>%
    inner_join(methy_regression_long) %>%
    inner_join(rna_pc_1_3)
}
  rna_regression_model <- rna_regression_combine %>%
    group_by(Gene_ID) %>%
    nest() %>%
    mutate(model = map(data, ~lm(rna ~ cnv + methy1 +methy2+methy3+methy4+methy5+methy6+methy7+methy8+methy9+methy10+methy11+methy12+methy13+methy14 + rna_pc1 + rna_pc2 + rna_pc3, data = .)))

  rna_regression_tidy <- rna_regression_model %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest

  rna_cnv <- rna_regression_tidy %>%
    filter(term == "cnv") %>%
    rename_if(is.numeric, ~paste0(., "_cnv_rna")) %>%
    dplyr::select(-term)

  rna_methy <- rna_regression_tidy %>%
    filter(term %in% paste0(rep("methy", 14), 1:14)) %>%
    rename_if(is.numeric, ~paste0(., "_methy_rna")) %>%
    dplyr::select(-term) %>%
    arrange(Gene_ID) %>%
    bind_cols(methy_regression %>% arrange(Gene_ID) %>% select(Hybridization))

  rna_info <- rna_cnv %>%
    inner_join(rna_methy)

  rna_regression_glance <- rna_regression_model %>%
    mutate(model_info = map(model, broom::glance)) %>%
    dplyr::select(-data, -model) %>%
    unnest %>%
    rename_if(is.numeric, ~paste0(., "_rna"))

  # protein ---------------------------------------------------------------------
  if (permutate == 2){
    protein_regression_combine <- protein_regression_long_permutated %>%
      inner_join(cnv_regression_long) %>%
      inner_join(methy_regression_long) %>%
      inner_join(protein_pc_1_3)
  } else {
    protein_regression_combine <- protein_regression_long %>%
      inner_join(cnv_regression_long) %>%
      inner_join(methy_regression_long) %>%
      inner_join(protein_pc_1_3)
  }

  protein_regression_model <- protein_regression_combine %>%
    group_by(Gene_ID) %>%
    nest() %>%
    mutate(model = map(data, ~lm(protein ~ cnv + methy1 +methy2+methy3+methy4+methy5+methy6+methy7+methy8+methy9+methy10+methy11+methy12+methy13+methy14 + protein_pc1 + protein_pc2 + protein_pc3, data = .)))

  protein_regression_tidy <- protein_regression_model %>%
    mutate(model_info = map(model, broom::tidy)) %>%
    dplyr::select(-data, -model) %>%
    unnest

  protein_cnv <- protein_regression_tidy %>%
    filter(term == "cnv") %>%
    rename_if(is.numeric, ~paste0(., "_cnv_protein")) %>%
    dplyr::select(-term)

  protein_methy <- protein_regression_tidy %>%
    filter(term %in% paste0(rep("methy", 14), 1:14)) %>%
    rename_if(is.numeric, ~paste0(., "_methy_protein")) %>%
    dplyr::select(-term) %>%
    arrange(Gene_ID) %>%
    bind_cols(methy_regression %>% arrange(Gene_ID) %>% select(Hybridization))

  protein_info <- protein_cnv %>%
    inner_join(protein_methy)

  protein_regression_glance <- protein_regression_model %>%
    mutate(model_info = map(model, broom::glance)) %>%
    dplyr::select(-data, -model) %>%
    unnest %>%
    rename_if(is.numeric, ~paste0(., "_protein"))

  # phospho ---------------------------------------------------------------------
  if (permutate == 3){
    phospho_regression_combine <- phospho_regression_long_permutated %>%
      inner_join(cnv_regression_long) %>%
      inner_join(methy_regression_long) %>%
      inner_join(phospho_pc_1_3)
  } else {
    phospho_regression_combine <- phospho_regression_long %>%
      inner_join(cnv_regression_long) %>%
      inner_join(methy_regression_long) %>%
      inner_join(phospho_pc_1_3)
  }

  phospho_regression_model <- phospho_regression_combine %>%
    group_by(Gene_ID, phospho_ID) %>%
    nest() %>%
    mutate(model_full = map(data, ~lm(phospho ~ cnv + methy1 +methy2+methy3+methy4+methy5+methy6+methy7+methy8+methy9+methy10+methy11+methy12+methy13+methy14 + phospho_pc1 + phospho_pc2 + phospho_pc3, data = .)),
           model_small = map(data, ~lm(phospho ~  phospho_pc1 + phospho_pc2 + phospho_pc3, data = .)))

  phospho_model_smallest_p <- phospho_regression_model %>%
    mutate(anova = map2(model_small, model_full, anova)) %>%
    dplyr::select(-data, -model_full, -model_small) %>%
    unnest() %>%
    drop_na() %>%
    group_by(Gene_ID) %>%
    top_n(-1, `Pr(>F)`) %>%
    ungroup

  phospho_model_final <- phospho_regression_model %>%
    dplyr::select(-model_small, -data) %>%
    inner_join(phospho_model_smallest_p) %>%
    filter(phospho_ID != "CALD1.NP_149129.2:s783")

  phospho_regression_tidy <- phospho_model_final %>%
    mutate(model_info = map(model_full, broom::tidy)) %>%
    dplyr::select(Gene_ID, phospho_ID, model_info ) %>%
    unnest

  phospho_cnv <- phospho_regression_tidy %>%
    filter(term == "cnv") %>%
    rename_if(is.numeric, ~paste0(., "_cnv_phospho")) %>%
    dplyr::select(-term)

  phospho_methy <- phospho_regression_tidy %>%
    filter(term %in% paste0(rep("methy", 14), 1:14)) %>%
    rename_if(is.numeric, ~paste0(., "_methy_phospho")) %>%
    dplyr::select(-term) %>%
    arrange(Gene_ID) %>%
    bind_cols(methy_regression %>% arrange(Gene_ID) %>% select(Hybridization))

  phospho_info <- phospho_cnv %>%
    inner_join(phospho_methy)

  phospho_regression_glance <- phospho_model_final %>%
    mutate(model_info = map(model_full, broom::glance)) %>%
    dplyr::select(Gene_ID, phospho_ID, model_info) %>%
    unnest %>%
    rename_if(is.numeric, ~paste0(., "_phospho"))

  # combine results ---------------------------------------------------------

  cnv_estimate <- rna_cnv %>%
    dplyr::select(Gene_ID, estimate_cnv_rna) %>%
    inner_join(protein_cnv %>%
                 dplyr::select(Gene_ID, estimate_cnv_protein)) %>%
    inner_join(phospho_cnv %>%
                 dplyr::select(phospho_ID, Gene_ID, estimate_cnv_phospho))

  cnv_std.error <- rna_cnv %>%
    dplyr::select(Gene_ID, std.error_cnv_rna) %>%
    inner_join(protein_cnv %>%
                 dplyr::select(Gene_ID, std.error_cnv_protein)) %>%
    inner_join(phospho_cnv %>%
                 dplyr::select(phospho_ID, Gene_ID, std.error_cnv_phospho))

  cnv_sigma2 <- rna_regression_glance %>%
    dplyr::select(Gene_ID, sigma2_rna = sigma_rna) %>%
    inner_join(protein_regression_glance %>%
                 dplyr::select(Gene_ID, sigma2_protein = sigma_protein)) %>%
    inner_join(phospho_regression_glance %>%
                 dplyr::select(phospho_ID, Gene_ID, sigma2_phospho = sigma_phospho)) %>%
    mutate_if(is.numeric, ~(.^2))


  cnv_df.residual <- rna_regression_glance %>%
    dplyr::select(Gene_ID, df.residual_rna) %>%
    inner_join(protein_regression_glance %>%
                 dplyr::select(Gene_ID, df.residual_protein)) %>%
    inner_join(phospho_regression_glance %>%
                 dplyr::select(phospho_ID, Gene_ID, df.residual_phospho))

  cnv_r.squared <- rna_regression_glance %>%
    dplyr::select(Gene_ID, r.squared_rna) %>%
    inner_join(protein_regression_glance %>%
                 dplyr::select(Gene_ID, r.squared_protein)) %>%
    inner_join(phospho_regression_glance %>%
                 dplyr::select(phospho_ID, Gene_ID, r.squared_phospho))

  cnv_v <- cnv_std.error %>%
    inner_join(rna_regression_glance %>%
                 dplyr::select(Gene_ID, sigma_rna) %>%
                 inner_join(protein_regression_glance %>%
                              dplyr::select(Gene_ID, sigma_protein)) %>%
                 inner_join(phospho_regression_glance %>%
                              dplyr::select(phospho_ID, Gene_ID, sigma_phospho))) %>%
    mutate(v_rna = (std.error_cnv_rna/sigma_rna) ^2,
           v_protein = (std.error_cnv_protein/sigma_protein) ^2,
           v_phospho = (std.error_cnv_phospho/sigma_phospho) ^2) %>%
    dplyr::select(phospho_ID, Gene_ID, starts_with("v"))

  # cbind(cnv_input_1_3$xName, cnv_input_1_3$v_g_J) %>%
  #   inner_join(cnv_v) %>%
  #   as.tibble %>%
  #   arrange(V1) %>%
  #   mutate(diff = V1-v_rna) %>%
  #   pull(diff) %>% sum


  cnv_statistic <- rna_cnv %>%
    dplyr::select(Gene_ID, statistic_cnv_rna) %>%
    inner_join(protein_cnv %>%
                 dplyr::select(Gene_ID, statistic_cnv_protein)) %>%
    inner_join(phospho_cnv %>%
                 dplyr::select(phospho_ID, Gene_ID, statistic_cnv_phospho))


  cnv_full <- cnv_std.error %>%
    inner_join(cnv_estimate) %>%
    inner_join(cnv_sigma2) %>%
    inner_join(cnv_df.residual) %>%
    inner_join(cnv_v) %>%
    inner_join(cnv_r.squared) %>%
    inner_join(cnv_statistic)

  cnv_result <- list(
    betas_J = cnv_full %>% dplyr::select(starts_with("estimate")) %>% as.matrix(),
    betas_se_J = cnv_full %>% dplyr::select(starts_with("std.error")) %>% as.matrix(),
    sigma2_J = cnv_full %>% dplyr::select(starts_with("sigma2")) %>% as.matrix(),
    dfs_J = cnv_full %>% dplyr::select(starts_with("df")) %>% as.matrix(),
    v_g_J = cnv_full %>% dplyr::select(starts_with("v")) %>% as.matrix(),
    r_square_J = cnv_full %>% dplyr::select(starts_with("r.squared")) %>% as.matrix(),
    t_J = cnv_full %>% dplyr::select(starts_with("statistic")) %>% as.matrix(),
    xName = cnv_full %>% dplyr::select(phospho_ID, Gene_ID),
    yName = cnv_full %>% dplyr::select(Gene_ID),
    name = rep("cnv", 676)
  )

  methy_estimate <- rna_methy %>%
    dplyr::select(Gene_ID, Hybridization, estimate_methy_rna) %>%
    inner_join(protein_methy %>%
                 dplyr::select(Gene_ID, Hybridization, estimate_methy_protein)) %>%
    inner_join(phospho_methy %>%
                 dplyr::select(Gene_ID, Hybridization,Gene_ID, estimate_methy_phospho))

  methy_std.error <- rna_methy %>%
    dplyr::select(Gene_ID, Hybridization, std.error_methy_rna) %>%
    inner_join(protein_methy %>%
                 dplyr::select(Gene_ID, Hybridization, std.error_methy_protein)) %>%
    inner_join(phospho_methy %>%
                 dplyr::select(Gene_ID, Hybridization,Gene_ID, std.error_methy_phospho))

  methy_sigma2 <- rna_regression_glance %>%
    dplyr::select(Gene_ID,sigma2_rna = sigma_rna) %>%
    inner_join(protein_regression_glance %>%
                 dplyr::select(Gene_ID, sigma2_protein = sigma_protein)) %>%
    inner_join(phospho_regression_glance %>%
                 dplyr::select(Gene_ID, phospho_ID, sigma2_phospho = sigma_phospho)) %>%
    mutate_if(is.numeric, ~(.^2))

  methy_df.residual <- rna_regression_glance %>%
    dplyr::select(Gene_ID, df.residual_rna) %>%
    inner_join(protein_regression_glance %>%
                 dplyr::select(Gene_ID, df.residual_protein)) %>%
    inner_join(phospho_regression_glance %>%
                 dplyr::select(Gene_ID, phospho_ID, df.residual_phospho))


  methy_r.squared <- rna_regression_glance %>%
    dplyr::select(Gene_ID, r.squared_rna) %>%
    inner_join(protein_regression_glance %>%
                 dplyr::select(Gene_ID, r.squared_protein)) %>%
    inner_join(phospho_regression_glance %>%
                 dplyr::select(Gene_ID, phospho_ID, r.squared_phospho))

  methy_v <- methy_std.error %>%
    inner_join(rna_regression_glance %>%
                 dplyr::select(Gene_ID, sigma_rna) %>%
                 inner_join(protein_regression_glance %>%
                              dplyr::select(Gene_ID, sigma_protein)) %>%
                 inner_join(phospho_regression_glance %>%
                              dplyr::select(Gene_ID,phospho_ID, sigma_phospho))) %>%
    mutate(v_rna = (std.error_methy_rna/sigma_rna) ^2,
           v_protein = (std.error_methy_protein/sigma_protein) ^2,
           v_phospho = (std.error_methy_phospho/sigma_phospho) ^2) %>%
    dplyr::select(Gene_ID, Hybridization,Gene_ID, starts_with("v"))

  methy_statistic <- rna_methy %>%
    dplyr::select(Gene_ID, Hybridization, statistic_methy_rna) %>%
    inner_join(protein_methy %>%
                 dplyr::select(Gene_ID, Hybridization, statistic_methy_protein)) %>%
    inner_join(phospho_methy %>%
                 dplyr::select(Gene_ID, Hybridization,Gene_ID, statistic_methy_phospho))


  methy_full <- methy_std.error %>%
    inner_join(methy_estimate) %>%
    inner_join(methy_sigma2) %>%
    inner_join(methy_df.residual) %>%
    inner_join(methy_v) %>%
    inner_join(methy_r.squared) %>%
    inner_join(methy_statistic)

  methy_result <- list(
    betas_J = methy_full %>% dplyr::select(starts_with("estimate")) %>% as.matrix(),
    betas_se_J = methy_full %>% dplyr::select(starts_with("std.error")) %>% as.matrix(),
    sigma2_J = methy_full %>% dplyr::select(starts_with("sigma2")) %>% as.matrix(),
    dfs_J = methy_full %>% dplyr::select(starts_with("df")) %>% as.matrix(),
    v_g_J = methy_full %>% dplyr::select(starts_with("v")) %>% as.matrix(),
    r_square_J = methy_full %>% dplyr::select(starts_with("r.squared")) %>% as.matrix(),
    t_J = methy_full %>% dplyr::select(starts_with("statistic")) %>% as.matrix(),
    xName = methy_full %>% dplyr::select(Gene_ID, Hybridization,Gene_ID),
    yName = methy_full %>% dplyr::select(Gene_ID, Hybridization),
    name = rep("methy", 1103)
)

# Primo ------------------------------------------------------------------
  cnv_iprofun <- MultiOmics_Input(cnv_result,pi1 = pi)
  methy_iprofun <- MultiOmics_Input(methy_result,pi1 = pi)
  iprofun_result <- list(cnv_iprofun, methy_iprofun)
  return(iprofun_result)
}

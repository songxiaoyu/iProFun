# Find overlapping sample for x, y and covariate---------------------------

rna_subid <- colnames(rna_normal)[grepl("TCGA", colnames(rna_normal))]
protein_subid <- colnames(protein_normal)[grepl("TCGA", colnames(protein_normal))]
phospho_subid <- colnames(phospho_normal)[grepl("TCGA", colnames(phospho_normal))]
cnv_subid <- colnames(cnv)[grepl("TCGA", colnames(cnv))]
methy_subid <- colnames(methy)[grepl("TCGA", colnames(methy))]


xrnaCommonSubID <- intersect(intersect(rna_subid, cnv_subid), methy_subid)
xproteinCommonSubID <- intersect(intersect(protein_subid, cnv_subid), methy_subid)
xphosphoCommonSubID <- intersect(intersect(phospho_subid, cnv_subid), methy_subid)
x_common_subid <- unique(c(xrnaCommonSubID, xproteinCommonSubID, xphosphoCommonSubID))
xyCommonSubID <- list(xrnaCommonSubID, xproteinCommonSubID, xphosphoCommonSubID)
x_common_subid <- unique(c(xrnaCommonSubID, xproteinCommonSubID))

cnv <- cnv %>%
  select("Gene_ID", x_common_subid)
methy <- methy %>%
  select("Gene_ID", "Hybridization", "chr", x_common_subid)

rna_normal <- rna_normal %>%
  select("Gene_ID", x_common_subid)

protein_normal <- protein_normal %>%
  select("Gene_ID", intersect(x_common_subid, names(protein_normal)))
phospho_normal <- phospho_normal %>%
  select("Gene_ID", "phospho_ID", intersect(x_common_subid, names(phospho_normal)))
# rna_normal <- rna_normal[apply(rna_normal, 1, function(f) mean(is.na(f))<=0.5),]
protein_normal <- protein_normal[apply(protein_normal, 1, function(f) mean(is.na(f))<=0.5),]
phospho_normal <- phospho_normal[apply(phospho_normal, 1, function(f) mean(is.na(f))<=0.5),]
# cnv <- cnv[apply(cnv, 1, function(f) mean(is.na(f))<=0.5),]
cnv <- cnv[complete.cases(cnv),]
# methy <- methy[apply(methy, 1, function(f) mean(is.na(f))<=0.5),]
methy <- methy[complete.cases(methy),]
methy_hyb <- methy_regression %>% pull(Hybridization)
# Find overlapping genes for x and y---------------------------
rnaGeneID <- data.frame(rna_normal)[,"Gene_ID"]
proteinGeneID <- data.frame(protein_normal)[,"Gene_ID"]
phosphoGeneID <- data.frame(phospho_normal)[,"Gene_ID"]
cnv_GeneID <- data.frame(cnv)[,"Gene_ID"]
methy_GeneID <- data.frame(methy)[,"Gene_ID"]
xy_common_geneID <- intersect(intersect(intersect(intersect(methy_GeneID, rnaGeneID), proteinGeneID), phosphoGeneID),cnv_GeneID)

rna_regression <- rna_normal %>%
  filter(Gene_ID %in% xy_common_geneID) %>%
  dplyr::select(Gene_ID, xrnaCommonSubID)

protein_regression <- protein_normal %>%
  filter(Gene_ID %in% xy_common_geneID) %>%
  dplyr::select(Gene_ID, xproteinCommonSubID)

phospho_regression <- phospho_normal %>%
  filter(Gene_ID %in% xy_common_geneID) %>%
  dplyr::select(Gene_ID, phospho_ID, xphosphoCommonSubID)

cnv_regression <- cnv %>%
  filter(Gene_ID %in% xy_common_geneID) %>%
  dplyr::select(Gene_ID, x_common_subid)

methy_regression <- methy %>%
  filter(Gene_ID %in% xy_common_geneID) %>%
  select("Gene_ID", "Hybridization", "chr", x_common_subid)

# change to long format ---------------------------------------------------


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


protein_pc_1_3 <- protein_pc_1_3 %>%
  mutate(PC = c("protein_pc1", "protein_pc2", "protein_pc3")) %>%
  gather(subject_id, pc,-PC) %>% as.tibble %>%
  spread(PC, pc)

rna_pc_1_3 <- rna_pc_1_3 %>%
  mutate(PC = c("rna_pc1", "rna_pc2", "rna_pc3")) %>%
  gather(subject_id, pc,-PC) %>% as.tibble %>%
  spread(PC, pc)

phospho_pc_1_3 <- phospho_pc_1_3 %>%
  mutate(PC = c("phospho_pc1", "phospho_pc2", "phospho_pc3")) %>%
  gather(subject_id, pc,-PC) %>% as.tibble %>%
  spread(PC, pc)

# rna ---------------------------------------------------------------------

rna_regression_combine <- rna_regression_long %>%
  inner_join(cnv_regression_long) %>%
  inner_join(methy_regression_long) %>%
  inner_join(rna_pc_1_3)

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

protein_regression_combine <- protein_regression_long %>%
  inner_join(cnv_regression_long) %>%
  inner_join(methy_regression_long) %>%
  inner_join(protein_pc_1_3)

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

phospho_regression_combine <- phospho_regression_long %>%
  inner_join(cnv_regression_long) %>%
  inner_join(methy_regression_long) %>%
  inner_join(phospho_pc_1_3)

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
# rna_cnv %>%
#   dplyr::select(Gene_ID, estimate_cnv_rna) %>%
#   inner_join(protein_cnv %>%
#                dplyr::select(Gene_ID, estimate_cnv_protein)) %>%
#   inner_join(phospho_cnv %>%
#                dplyr::select(Gene_ID, estimate_cnv_phospho)) %>%
#   inner_join(cbind(cnv_input_1_3$xName, cnv_input_1_3$betas_J)) %>%
#   dplyr::select(estimate_cnv_rna, V1) %>%
#   mutate(diff = estimate_cnv_rna - V1) %>% pull(diff) %>% sum
#
# sum(rna_methy %>% pull(2)) - sum(methy_input_1_3$betas_J[,1])

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

# cbind(cnv_input_1_3$xName,cnv_input_1_3$sigma2_J) %>%
#   inner_join(cnv_sigma2) %>% as.tibble()
#   mutate(diff = V3-sigma2_phospho) %>% pull(diff) %>% sum

cnv_df.residual <- rna_regression_glance %>%
  dplyr::select(Gene_ID, df.residual_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, df.residual_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(phospho_ID, Gene_ID, df.residual_phospho))

# cbind(cnv_input_1_3$xName, cnv_input_1_3$dfs_J) %>%
#   inner_join(cnv_df.residual) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V3-df.residual_phospho) %>% pull(diff) %>% sum

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

# save(cnv_result, file = "CCRCC/RData/cnv_regression_unimputed_x_filterNA_50%_new_phospho.RData")




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
  dplyr::select(Gene_ID, Hybridization, sigma2_rna = sigma_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, Hybridization, sigma2_protein = sigma_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(Gene_ID, Hybridization,Gene_ID, sigma2_phospho = sigma_phospho)) %>%
  mutate_if(is.numeric, ~(.^2))

# cbind(methy_input_1_3$xName,methy_input_1_3$sigma2_J) %>%
#   inner_join(methy_sigma2) %>% as.tibble()
#   mutate(diff = V3-sigma2_phospho) %>% pull(diff) %>% sum

methy_df.residual <- rna_regression_glance %>%
  dplyr::select(Gene_ID, Hybridization, df.residual_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, Hybridization, df.residual_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(Gene_ID, Hybridization,Gene_ID, df.residual_phospho))

# cbind(methy_input_1_3$xName, methy_input_1_3$dfs_J) %>%
#   inner_join(methy_df.residual) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V3-df.residual_phospho) %>% pull(diff) %>% sum

methy_r.squared <- rna_regression_glance %>%
  dplyr::select(Gene_ID, Hybridization, r.squared_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, Hybridization, r.squared_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(Gene_ID, Hybridization,Gene_ID, r.squared_phospho))

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

cbind(methy_input_1_3$xName, methy_input_1_3$v_g_J) %>%
  inner_join(methy_v) %>%
  as.tibble %>%
  arrange(`1`) %>%
  mutate(diff = `1`-v_rna) %>%
  pull(diff) %>% sum


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

cbind(methy_input_1_3$xName,methy_input_1_3$betas_se_J) %>%
  inner_join(methy_std.error) %>% as.tibble() %>%
  mutate(diff = `3` - std.error_methy_phospho) %>% pull(diff) %>% sum
cbind(methy_input_1_3$xName,methy_input_1_3$sigma2_J) %>%
  inner_join(methy_sigma2) %>% as.tibble()
cbind(methy_input_1_3$xName, methy_input_1_3$v_g_J) %>%
  inner_join(methy_v) %>%
  as.tibble %>% mutate(diff = v_phospho - `3`) %>% pull(diff) %>% sum ###?????

### And the methy hyb

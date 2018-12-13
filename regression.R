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
xyCommonSubID <- list(xrnaCommonSubID, xproteinCommonSubID)

# Find overlapping genes for x and y---------------------------
rnaGeneID <- data.frame(rna_normal)[,"Gene_ID"]
proteinGeneID <- data.frame(protein_normal)[,"Gene_ID"]
phosphoGeneID <- data.frame(phospho_normal)[,"Gene_ID"]

methy_GeneID <- data.frame(methy)[,"Gene_ID"]
xy_common_geneID <- intersect(intersect(intersect(methy_GeneID, rnaGeneID), proteinGeneID), phosphoGeneID)

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
  dplyr::select(Gene_ID, x_common_subid)

# change to long format ---------------------------------------------------


rna_regression_long <- rna_regression %>%
  gather(subject_id, rna, - Gene_ID) %>% as.tibble()

protein_regression_long <- protein_regression %>%
  gather(subject_id, protein, - Gene_ID) %>% as.tibble()

phospho_regression_long <- phospho_regression %>%
  gather(subject_id, phospho, - Gene_ID, - index, - peptide) %>% as.tibble()

cnv_regression_long <- cnv_regression %>%
  gather(subject_id, cnv, - Gene_ID) %>% as.tibble()

methy_regression_long <- methy_regression %>%
  gather(subject_id, methy, - Gene_ID) %>% as.tibble()

protein_pc_1_3 %>% gather(subject_id, pc_1_3) %>% as.tibble

regression_long <- rna_regression_long %>%
  inner_join(protein_regression_long) %>%
  inner_join(cnv_baf_regression_long) %>%
  inner_join(cnv_regression_long) %>%
  inner_join(methy_regression_long) %>%
  inner_join(purity_tumor_long) %>%
  inner_join(age_long) %>%
  inner_join(gender_long)

regression_long_with_mutation <- regression_long %>%
  inner_join(mutation_reg_111_long)


# rna ---------------------------------------------------------------------
rna_mutation_model <- regression_long_with_mutation %>%
  group_by(Gene_ID) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rna ~ cnv + cnv_baf + methy + purity + age + gender + age + gender + mutation, data = .)))

rna_mutation_model_gene <- rna_mutation_model %>% pull(1)

rna_regression_model <- regression_long %>%
  group_by(Gene_ID) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rna ~ cnv + cnv_baf + methy + purity + age + gender + age + gender, data = .)))

rna_regression_model <- rna_regression_model %>%
  filter(!(Gene_ID %in% rna_mutation_model_gene)) %>%
  bind_rows(rna_mutation_model)


rna_regression_tidy <- rna_regression_model %>%
  mutate(model_info = map(model, broom::tidy)) %>%
  dplyr::select(-data, -model) %>%
  unnest

rna_cnv <- rna_regression_tidy %>%
  filter(term == "cnv") %>%
  rename_if(is.numeric, ~paste0(., "_cnv_rna")) %>%
  dplyr::select(-term)

rna_cnv_baf <- rna_regression_tidy %>%
  filter(term == "cnv_baf") %>%
  rename_if(is.numeric, ~paste0(., "_cnv_baf_rna")) %>%
  dplyr::select(-term)

rna_methy <- rna_regression_tidy %>%
  filter(term == "methy") %>%
  rename_if(is.numeric, ~paste0(., "_methy_rna")) %>%
  dplyr::select(-term)

rna_mutation_regression_tidy <- rna_mutation_model %>%
  mutate(model_info = map(model, broom::tidy)) %>%
  dplyr::select(Gene_ID, model_info) %>%
  unnest

rna_mutation <- rna_mutation_regression_tidy %>%
  filter(term == "mutation") %>%
  rename_if(is.numeric, ~paste0(., "_mutation_rna")) %>%
  dplyr::select(-term)

rna_info <- rna_cnv %>%
  inner_join(rna_cnv_baf) %>%
  inner_join(rna_methy)

rna_mutation_regression_glance <- rna_mutation_model %>%
  mutate(model_info = map(model, broom::glance)) %>%
  dplyr::select(Gene_ID, model_info) %>%
  unnest %>%
  rename_if(is.numeric, ~paste0(., "_rna"))

rna_regression_glance <- rna_regression_model %>%
  mutate(model_info = map(model, broom::glance)) %>%
  dplyr::select(-data, -model) %>%
  unnest %>%
  rename_if(is.numeric, ~paste0(., "_rna"))

# protein ---------------------------------------------------------------------
protein_mutation_model <- regression_long_with_mutation %>%
  group_by(Gene_ID) %>%
  nest() %>%
  mutate(model = map(data, ~lm(protein ~ cnv + cnv_baf + methy + purity + age + gender + age + gender + mutation, data = .)))

protein_mutation_model_gene <- protein_mutation_model %>% pull(1)

protein_regression_model <- regression_long %>%
  group_by(Gene_ID) %>%
  nest() %>%
  mutate(model = map(data, ~lm(protein ~ cnv + cnv_baf + methy + purity + age + gender + age + gender, data = .)))

protein_regression_model <- protein_regression_model %>%
  filter(!(Gene_ID %in% protein_mutation_model_gene)) %>%
  bind_rows(protein_mutation_model)


protein_regression_tidy <- protein_regression_model %>%
  mutate(model_info = map(model, broom::tidy)) %>%
  dplyr::select(-data, -model) %>%
  unnest

protein_cnv <- protein_regression_tidy %>%
  filter(term == "cnv") %>%
  rename_if(is.numeric, ~paste0(., "_cnv_protein")) %>%
  dplyr::select(-term)

protein_cnv_baf <- protein_regression_tidy %>%
  filter(term == "cnv_baf") %>%
  rename_if(is.numeric, ~paste0(., "_cnv_baf_protein")) %>%
  dplyr::select(-term)

protein_methy <- protein_regression_tidy %>%
  filter(term == "methy") %>%
  rename_if(is.numeric, ~paste0(., "_methy_protein")) %>%
  dplyr::select(-term)

protein_mutation_regression_tidy <- protein_mutation_model %>%
  mutate(model_info = map(model, broom::tidy)) %>%
  dplyr::select(Gene_ID, model_info) %>%
  unnest

protein_mutation <- protein_mutation_regression_tidy %>%
  filter(term == "mutation") %>%
  rename_if(is.numeric, ~paste0(., "_mutation_protein")) %>%
  dplyr::select(-term)

protein_info <- protein_cnv %>%
  inner_join(protein_cnv_baf) %>%
  inner_join(protein_methy)

protein_mutation_regression_glance <- protein_mutation_model %>%
  mutate(model_info = map(model, broom::glance)) %>%
  dplyr::select(Gene_ID, model_info) %>%
  unnest %>%
  rename_if(is.numeric, ~paste0(., "_protein"))

protein_regression_glance <- protein_regression_model %>%
  mutate(model_info = map(model, broom::glance)) %>%
  dplyr::select(-data, -model) %>%
  unnest %>%
  rename_if(is.numeric, ~paste0(., "_protein"))


# phospho -----------------------------------------------------------------

regression_long_with_phospho <- regression_long %>%
  inner_join(phospho_regression_long)

regression_long_with_mutation_phospho <- regression_long_with_phospho %>%
  inner_join(mutation_reg_111_long)

# regression_long_with_mutation_phospho %>%
#   janitor::get_dupes(Gene_ID, subject_id) %>%
#   arrange(desc(dupe_count))

phospho_mutation_model <- regression_long_with_mutation_phospho %>%
  group_by(index, peptide, Gene_ID) %>%
  nest() %>%
  mutate(data = map(data, ~na.omit(.))) %>%
  mutate(model_full = map(data, ~lm(phospho ~ cnv + cnv_baf + methy + purity + age + gender + age + gender + mutation, data = .)),
         model_small = map(data, ~lm(phospho ~ purity +age + gender + mutation, data = .)))

phospho_mutation_model_smallest_p <- phospho_mutation_model %>%
  mutate(anova = map2(model_small, model_full, anova)) %>%
  dplyr::select(-data, -model_full, -model_small) %>%
  unnest() %>%
  drop_na() %>%
  group_by(Gene_ID) %>%
  top_n(-1, `Pr(>F)`) %>%
  ungroup

phospho_mutation_model_final <- phospho_mutation_model %>%
  dplyr::select(-model_small, -data) %>%
  inner_join(phospho_mutation_model_smallest_p)

phospho_mutation_model_gene <- phospho_mutation_model_final %>% pull(1)

phospho_regression_model <- regression_long_with_phospho %>%
  group_by(index, peptide, Gene_ID) %>%
  nest() %>%
  mutate(data = map(data, ~na.omit(.))) %>%
  mutate(model_full = map(data, ~lm(phospho ~ cnv + cnv_baf + methy + purity + age + gender + age + gender, data = .)),
         model_small = map(data, ~lm(phospho ~ purity +age + gender, data = .)))

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
  inner_join(phospho_model_smallest_p)

phospho_regression_model_final <- phospho_model_final %>%
  filter(!(index %in% phospho_mutation_model_gene)) %>%
  bind_rows(phospho_mutation_model_final)

phospho_regression_tidy <- phospho_regression_model_final %>%
  mutate(model_info = map(model_full, broom::tidy)) %>%
  dplyr::select(index, peptide, Gene_ID, model_info) %>%
  unnest

phospho_cnv <- phospho_regression_tidy %>%
  filter(term == "cnv") %>%
  rename_if(is.numeric, ~paste0(., "_cnv_phospho")) %>%
  dplyr::select(-term)

phospho_cnv_baf <- phospho_regression_tidy %>%
  filter(term == "cnv_baf") %>%
  rename_if(is.numeric, ~paste0(., "_cnv_baf_phospho")) %>%
  dplyr::select(-term)

phospho_methy <- phospho_regression_tidy %>%
  filter(term == "methy") %>%
  rename_if(is.numeric, ~paste0(., "_methy_phospho")) %>%
  dplyr::select(-term)

phospho_mutation_regression_tidy <- phospho_mutation_model_final %>%
  mutate(model_info = map(model_full, broom::tidy)) %>%
  dplyr::select(index, peptide, Gene_ID, model_info) %>%
  unnest

phospho_mutation <- phospho_mutation_regression_tidy %>%
  filter(term == "mutation") %>%
  rename_if(is.numeric, ~paste0(., "_mutation_phospho")) %>%
  dplyr::select(-term)

phospho_info <- phospho_cnv %>%
  inner_join(phospho_cnv_baf) %>%
  inner_join(phospho_methy)

phospho_mutation_regression_glance <- phospho_mutation_model_final %>%
  mutate(model_info = map(model_full, broom::glance)) %>%
  dplyr::select(index, peptide, Gene_ID, model_info) %>%
  unnest %>%
  rename_if(is.numeric, ~paste0(., "_phospho"))

phospho_regression_glance <- phospho_regression_model_final %>%
  mutate(model_info = map(model_full, broom::glance)) %>%
  dplyr::select(index, peptide, Gene_ID, model_info) %>%
  unnest %>%
  rename_if(is.numeric, ~paste0(., "_phospho"))

# tibble(a = 1:10, b = 10:1, c = c(rep(1, 5), rep(2, 5))) %>%
#   group_by(c) %>%
#   top_n(1, a)


# combine results ---------------------------------------------------------
# rna_info %>%
#   dplyr::select(Gene_ID, estimate_cnv_rna) %>%
#   inner_join(protein_info %>%
#                dplyr::select(Gene_ID, estimate_cnv_protein)) %>%
#   inner_join(phospho_info %>%
#                dplyr::select(Gene_ID, estimate_cnv_phospho)) %>%
#   inner_join(cbind(multireg_lr$xName, multireg_lr$betas_J)) %>%
#   dplyr::select(estimate_cnv_rna, V1) %>%
#   mutate(diff = estimate_cnv_rna - V1) %>% pull(diff) %>% sum
#
# rna_info %>%
#   dplyr::select(Gene_ID, estimate_cnv_rna) %>%
#   inner_join(protein_info %>%
#                dplyr::select(Gene_ID, estimate_cnv_protein)) %>%
#   inner_join(phospho_info %>%
#                dplyr::select(Gene_ID, estimate_cnv_phospho)) %>%
#   inner_join(cbind(multireg_lr$xName, multireg_lr$betas_J)) %>%
#   dplyr::select(estimate_cnv_phospho, V3) %>%
#   mutate(diff = estimate_cnv_phospho - V3) %>% pull(diff) %>% sum


cnv_estimate <- rna_info %>%
  dplyr::select(Gene_ID, estimate_cnv_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, estimate_cnv_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, estimate_cnv_phospho))

cnv_std.error <- rna_info %>%
  dplyr::select(Gene_ID, std.error_cnv_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, std.error_cnv_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, std.error_cnv_phospho))

cnv_sigma2 <- rna_regression_glance %>%
  dplyr::select(Gene_ID, sigma2_rna = sigma_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, sigma2_protein = sigma_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, sigma2_phospho = sigma_phospho)) %>%
  mutate_if(is.numeric, ~(.^2))

# cbind(multireg_lr$xName,multireg_lr$sigma2_J) %>%
#   inner_join(cnv_sigma2) %>%
#   mutate(diff = V3-sigma2_phospho) %>% pull(diff) %>% sum

cnv_df.residual <- rna_regression_glance %>%
  dplyr::select(Gene_ID, df.residual_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, df.residual_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, df.residual_phospho))

# cbind(multireg_lr$xName, multireg_lr$dfs_J) %>%
#   inner_join(cnv_df.residual) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V3-df.residual_phospho) %>% pull(diff) %>% sum

cnv_r.squared <- rna_regression_glance %>%
  dplyr::select(Gene_ID, r.squared_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, r.squared_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, r.squared_phospho))

# cbind(multireg_lr$xName, multireg_lr$r_square_J) %>%
#   inner_join(cnv_r.squared) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-r.squared_rna) %>% pull(diff) %>% sum

cnv_v <- cnv_std.error %>%
  inner_join(rna_regression_glance %>%
               dplyr::select(Gene_ID, sigma_rna) %>%
               inner_join(protein_regression_glance %>%
                            dplyr::select(Gene_ID, sigma_protein)) %>%
               inner_join(phospho_regression_glance %>%
                            dplyr::select(peptide, Gene_ID, sigma_phospho))) %>%
  mutate(v_rna = (std.error_cnv_rna/sigma_rna) ^2,
         v_protein = (std.error_cnv_protein/sigma_protein) ^2,
         v_phospho = (std.error_cnv_phospho/sigma_phospho) ^2) %>%
  dplyr::select(peptide, Gene_ID, starts_with("v"))

# cbind(multireg_lr$xName, multireg_lr$v_g_J) %>%
#   inner_join(cnv_v) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-v_rna) %>%
#   arrange(diff)


cnv_statistic <- rna_info %>%
  dplyr::select(Gene_ID, statistic_cnv_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, statistic_cnv_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, statistic_cnv_phospho))

# cbind(multireg_lr$xName, multireg_lr$t_J) %>%
#   inner_join(cnv_statistic) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-statistic_cnv_rna) %>% pull(diff) %>% sum

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
  xName = cnv_full %>% dplyr::select(peptide, Gene_ID),
  yName = cnv_full %>% dplyr::select(Gene_ID),
  name = rep("cnv", 4021)
)

save(cnv_result, file = "CCRCC/RData/cnv_regression_unimputed_x_filterNA_50%_new_phospho.RData")

cnv_baf_estimate <- rna_info %>%
  dplyr::select(Gene_ID, estimate_cnv_baf_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, estimate_cnv_baf_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, estimate_cnv_baf_phospho))

cnv_baf_std.error <- rna_info %>%
  dplyr::select(Gene_ID, std.error_cnv_baf_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, std.error_cnv_baf_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, std.error_cnv_baf_phospho))

cnv_baf_sigma2 <- rna_regression_glance %>%
  dplyr::select(Gene_ID, sigma2_rna = sigma_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, sigma2_protein = sigma_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, sigma2_phospho = sigma_phospho)) %>%
  mutate_if(is.numeric, ~(.^2))

# cbind(multireg_lr$xName,multireg_lr$sigma2_J) %>%
#   inner_join(cnv_baf_sigma2) %>%
#   mutate(diff = V3-sigma2_phospho) %>% pull(diff) %>% sum

cnv_baf_df.residual <- rna_regression_glance %>%
  dplyr::select(Gene_ID, df.residual_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, df.residual_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, df.residual_phospho))

# cbind(multireg_lr$xName, multireg_lr$dfs_J) %>%
#   inner_join(cnv_baf_df.residual) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V3-df.residual_phospho) %>% pull(diff) %>% sum

cnv_baf_r.squared <- rna_regression_glance %>%
  dplyr::select(Gene_ID, r.squared_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, r.squared_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, r.squared_phospho))

# cbind(multireg_lr$xName, multireg_lr$r_square_J) %>%
#   inner_join(cnv_baf_r.squared) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-r.squared_rna) %>% pull(diff) %>% sum

cnv_baf_v <- cnv_baf_std.error %>%
  inner_join(rna_regression_glance %>%
               dplyr::select(Gene_ID, sigma_rna) %>%
               inner_join(protein_regression_glance %>%
                            dplyr::select(Gene_ID, sigma_protein)) %>%
               inner_join(phospho_regression_glance %>%
                            dplyr::select(peptide, Gene_ID, sigma_phospho))) %>%
  mutate(v_rna = (std.error_cnv_baf_rna/sigma_rna) ^2,
         v_protein = (std.error_cnv_baf_protein/sigma_protein) ^2,
         v_phospho = (std.error_cnv_baf_phospho/sigma_phospho) ^2) %>%
  dplyr::select(peptide, Gene_ID, starts_with("v"))

# cbind(multireg_lr$xName, multireg_lr$v_g_J) %>%
#   inner_join(cnv_baf_v) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-v_rna) %>%
#   arrange(diff)


cnv_baf_statistic <- rna_info %>%
  dplyr::select(Gene_ID, statistic_cnv_baf_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, statistic_cnv_baf_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, statistic_cnv_baf_phospho))

# cbind(multireg_lr$xName, multireg_lr$t_J) %>%
#   inner_join(cnv_baf_statistic) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-statistic_cnv_baf_rna) %>% pull(diff) %>% sum

cnv_baf_full <- cnv_baf_std.error %>%
  inner_join(cnv_baf_estimate) %>%
  inner_join(cnv_baf_sigma2) %>%
  inner_join(cnv_baf_df.residual) %>%
  inner_join(cnv_baf_v) %>%
  inner_join(cnv_baf_r.squared) %>%
  inner_join(cnv_baf_statistic)

cnv_baf_result <- list(
  betas_J = cnv_baf_full %>% dplyr::select(starts_with("estimate")) %>% as.matrix(),
  betas_se_J = cnv_baf_full %>% dplyr::select(starts_with("std.error")) %>% as.matrix(),
  sigma2_J = cnv_baf_full %>% dplyr::select(starts_with("sigma2")) %>% as.matrix(),
  dfs_J = cnv_baf_full %>% dplyr::select(starts_with("df")) %>% as.matrix(),
  v_g_J = cnv_baf_full %>% dplyr::select(starts_with("v")) %>% as.matrix(),
  r_square_J = cnv_baf_full %>% dplyr::select(starts_with("r.squared")) %>% as.matrix(),
  t_J = cnv_baf_full %>% dplyr::select(starts_with("statistic")) %>% as.matrix(),
  xName = cnv_baf_full %>% dplyr::select(peptide, Gene_ID),
  yName = cnv_baf_full %>% dplyr::select(Gene_ID),
  name = rep("cnv_baf", 4021)
)

save(cnv_baf_result, file = "CCRCC/RData/cnv_baf_regression_unimputed_x_filterNA_50%_new_phospho.RData")

methy_estimate <- rna_info %>%
  dplyr::select(Gene_ID, estimate_methy_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, estimate_methy_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, estimate_methy_phospho))

methy_std.error <- rna_info %>%
  dplyr::select(Gene_ID, std.error_methy_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, std.error_methy_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, std.error_methy_phospho))

methy_sigma2 <- rna_regression_glance %>%
  dplyr::select(Gene_ID, sigma2_rna = sigma_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, sigma2_protein = sigma_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, sigma2_phospho = sigma_phospho)) %>%
  mutate_if(is.numeric, ~(.^2))

# cbind(multireg_lr$xName,multireg_lr$sigma2_J) %>%
#   inner_join(methy_sigma2) %>%
#   mutate(diff = V3-sigma2_phospho) %>% pull(diff) %>% sum

methy_df.residual <- rna_regression_glance %>%
  dplyr::select(Gene_ID, df.residual_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, df.residual_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, df.residual_phospho))

# cbind(multireg_lr$xName, multireg_lr$dfs_J) %>%
#   inner_join(methy_df.residual) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V3-df.residual_phospho) %>% pull(diff) %>% sum

methy_r.squared <- rna_regression_glance %>%
  dplyr::select(Gene_ID, r.squared_rna) %>%
  inner_join(protein_regression_glance %>%
               dplyr::select(Gene_ID, r.squared_protein)) %>%
  inner_join(phospho_regression_glance %>%
               dplyr::select(peptide, Gene_ID, r.squared_phospho))

# cbind(multireg_lr$xName, multireg_lr$r_square_J) %>%
#   inner_join(methy_r.squared) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-r.squared_rna) %>% pull(diff) %>% sum

methy_v <- methy_std.error %>%
  inner_join(rna_regression_glance %>%
               dplyr::select(Gene_ID, sigma_rna) %>%
               inner_join(protein_regression_glance %>%
                            dplyr::select(Gene_ID, sigma_protein)) %>%
               inner_join(phospho_regression_glance %>%
                            dplyr::select(peptide, Gene_ID, sigma_phospho))) %>%
  mutate(v_rna = (std.error_methy_rna/sigma_rna) ^2,
         v_protein = (std.error_methy_protein/sigma_protein) ^2,
         v_phospho = (std.error_methy_phospho/sigma_phospho) ^2) %>%
  dplyr::select(peptide, Gene_ID, starts_with("v"))

# cbind(multireg_lr$xName, multireg_lr$v_g_J) %>%
#   inner_join(methy_v) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-v_rna) %>%
#   arrange(diff)


methy_statistic <- rna_info %>%
  dplyr::select(Gene_ID, statistic_methy_rna) %>%
  inner_join(protein_info %>%
               dplyr::select(Gene_ID, statistic_methy_protein)) %>%
  inner_join(phospho_info %>%
               dplyr::select(peptide, Gene_ID, statistic_methy_phospho))

# cbind(multireg_lr$xName, multireg_lr$t_J) %>%
#   inner_join(methy_statistic) %>%
#   as.tibble %>%
#   arrange(V1) %>%
#   mutate(diff = V1-statistic_methy_rna) %>% pull(diff) %>% sum

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
  xName = methy_full %>% dplyr::select(peptide, Gene_ID),
  yName = methy_full %>% dplyr::select(Gene_ID),
  name = rep("methy", 4021)
)

save(methy_result, file = "CCRCC/RData/methy_regression_unimputed_x_filterNA_50%_new_phospho.RData")

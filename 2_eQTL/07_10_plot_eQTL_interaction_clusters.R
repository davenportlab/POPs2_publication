# 07_10_plot_eQTL_interaction_clusters.R

################################################################################

# 7.10. Plot eQTL interaction clusters

################################################################################

# Aim: Plot eQTL interaction clusters

########################### Paths ##########################
path = "genotyping/analysis/"
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Load packages ###########################
library(tidyverse)

########################### Load data ###########################
# time interaction clusters
time_int_clust <- read.csv(paste0(path, "eQTL/output_data/interactions_flanders_2/time_interaction_clusters_2.csv")) %>%
  mutate(cluster_label = cluster_label + 1) %>%
  separate(gsp, sep = "_", into = c("Gene", "SNP"), remove = FALSE, extra = "merge")
# data to plot time interaction eQTL
load(paste0(path, "eQTL/input_data/eqtl_int_files_flanders.rda"))
# time interaction eQTL
time_interactions <- read.csv(paste0(path,"eQTL/output_data/interactions_flanders_2/time_interactions_flanders_2.csv"))
# gene id_name converter 
gtf <- rtracklayer::import(gtf_inpath)
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

########################### Analysis ###########################
n_genes_per_cluster <- time_int_clust %>%
  group_by(cluster_label) %>% 
  summarize(n = n()) 

cluster_ibeta_plt_df <- time_int_clust %>%
  left_join(n_genes_per_cluster, by = "cluster_label") %>%
  pivot_longer(cols = c(wk12, wk20, wk28, wk36), names_to = "timepoint", values_to = "scaled_beta_sum") %>%
  mutate(cluster_label = as.factor(paste0("Cluster t", cluster_label, ", n eQTL = ", n)),
         timepoint = gsub("wk", "", timepoint))

# mean variance per time-point per cluster
mean_per_clust <- cluster_ibeta_plt_df %>%
  # add column with cluster name
  group_by(cluster_label, timepoint) %>% 
  summarize(cluster_tp_mean = mean(scaled_beta_sum))

plot_main_effect <- function(n.peer, signal, id_name, interaction_results){
  
  gene <- interaction_results$Gene[signal]
  SNP <- interaction_results$SNP[signal]
  signals <- pairs.int[pairs.int$Gene == gene, ]
  other.snps <- signals$SNP[signals$SNP != SNP]

  if(nrow(signals) > 1) {
    # run the model
    model.test <- lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                         0 + # fake intercept - real intercept = timepoint
                         geno.int[, SNP] + # genotype
                         timepoint + # timepoint
                         #outcome + # outcome 
                         covariates + # covariates
                         peer_factors[, 1:n.peer] + # peer factors
                         #geno.int[, SNP]*outcome + # interaction term 
                         geno.int[, other.snps] +  # other snps 
                         (1|individual) + (1|batch), # individual and batch 
                       REML=FALSE)
  } else {
    # run the model
    model.test <- lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                         0 + # fake intercept - real intercept = timepoint
                         geno.int[, SNP] + # genotype
                         timepoint + # timepoint
                         #outcome + # outcome 
                         covariates + # covariates
                         peer_factors[, 1:n.peer] + # peer factors
                         #geno.int[, SNP]*outcome + # interaction term 
                         (1|individual) + (1|batch), # individual and batch 
                       REML=FALSE)
  }
  
  
  # Get row indices used in the model
  model_rows <- as.numeric(rownames(model.frame(model.test)))
  
  # Extract random intercepts for individual and batch (example)
  re_ind <- ranef(model.test)$individual
  re_batch <- ranef(model.test)$batch
  # make sure the order is correct
  re_ind_contrib <- re_ind[individual[model_rows], "(Intercept)"]
  re_batch_contrib <- re_batch[batch[model_rows], "(Intercept)"]
  # sum random effects
  random_effects <- re_ind_contrib + re_batch_contrib
  
  # fixed effects
  X <- model.matrix(model.test)
  betas <- fixef(model.test)
  # remove the terms we want to keep 
  keep_terms <- c("geno.int[, SNP]", 
                  "timepoint12_weeks", "timepoint20_weeks", "timepoint28_weeks", "timepoint36_weeks")#, "outcomeTRUE", "geno.int[, SNP]:outcomeTRUE")
  keep_idx <- which(names(betas) %in% keep_terms)
  non_interest_idx <- setdiff(seq_along(betas), keep_idx)
  
  X_non_interest <- X[, non_interest_idx, drop=FALSE]
  fitted_non_interest <- X_non_interest %*% betas[non_interest_idx]
  
  X_keep <- X[, keep_idx, drop = FALSE]                # Subset design matrix to terms of interest
  interest_effect <- X_keep %*% betas[keep_idx]        # Multiply by their fixed-effect coefficients
  
  observed <- getME(model.test, "y")
  partial_resid <- observed - fitted_non_interest - random_effects
  
  # check this is working as expected
  # partial_resid = manually computed partial residuals
  # interest_effect = fitted values from terms of interest
  # resid(model.test) = model residuals
  # all.equal(partial_resid, (interest_effect + resid(model.test)), tolerance = 1e-6)
  
  # Subset all vectors on the same index
  df <- data.frame(
    geno = geno.int[model_rows, SNP],
    #outcome = outcome[model_rows],
    #interaction = geno.int[idx, SNP] * sga[idx],  # recompute interaction here
    partial_resid = partial_resid
  )
  
  half_dodge <- 0.2
  
  gene_name <- id_name$gtf.gene_name[id_name$gtf.gene_id == gene]
  
  # pval = ifelse(interaction_results$BH_adj_pval[signal] < 0.001, 
  #               format(interaction_results$BH_adj_pval[signal], scientific = TRUE, digits = 2),
  #               round(interaction_results$BH_adj_pval[signal], digits = 2))
  
  plt <- df %>%
    mutate(geno = as.numeric(geno)#,  # ensure 0,1,2
           #outcome = as.factor(outcome),
           #group_label = interaction(snp_base, outcome),
           #x_pos = case_when(
           #   outcome == levels(outcome)[1] ~ snp_base - half_dodge,
           #   outcome == levels(outcome)[2] ~ snp_base + half_dodge
           #)) %>%
    ) %>%
    ggplot(aes(x = geno, y = partial_resid)) +
    geom_jitter(#aes(color = outcome),
      #color = "black",
      width = 0.1, height = 0, size = 0.5, alpha = 1) +
    geom_boxplot(aes(group = geno),#, color = outcome),
                 outlier.shape = NA, width = 0.3, alpha = 0.5) +
    geom_smooth(#aes(color = outcome), 
      color = "black",
      method = "lm", formula = y ~ x, se = FALSE) +
    scale_x_continuous(breaks = 0:2, labels = c("0", "1", "2")) + 
    # scale_color_manual(values = c("FALSE" = "#73A2C6", "TRUE" = "#F4777F"),
    #                    labels = c(paste0("No ", outcome_name), outcome_name)) + 
    xlab(paste("N", SNP, "minor alleles")) + 
    ylab(paste(gene_name, "expression")) +
    #labs(color = NULL) + 
    # annotate("text", x = -0.4, y = min(df$partial_resid), 
    #          label = paste0("βG: ", 
    #                         round(interaction_results$eQTL_beta[signal], digits = 2),
    #                         "\nβGxO: ",
    #                         round(interaction_results$Interaction_beta[signal], digits = 2),
    #                         "\npint: ",
    #                         pval), hjust = 0, vjust = 0) + 
    ggtitle("Main effect eQTL") + 
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1))
  
  return(plt)
}

plot_interaction_time <- function(n.peer, signal, id_name, interaction_results){
  
  gene <- interaction_results$Gene[signal]
  SNP <- interaction_results$SNP[signal]
  signals <- pairs.int[pairs.int$Gene == gene, ]
  other.snps <- signals$SNP[signals$SNP != SNP]

  if(nrow(signals) > 1) {
    # run the model
    model.test <- lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                         0 + # fake intercept - real intercept = timepoint
                         geno.int[, SNP] + # genotype
                         timepoint + # timepoint
                         covariates + # covariates
                         peer_factors[, 1:n.peer] + # peer factors
                         geno.int[, SNP]*timepoint + # interaction term 
                         geno.int[, other.snps] +  # other snps 
                         (1|individual) + (1|batch), # individual and batch 
                       REML=FALSE)
  } else {
    # run the model
    model.test <- lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                         0 + # fake intercept - real intercept = timepoint
                         geno.int[, SNP] + # genotype
                         timepoint + # timepoint
                         covariates + # covariates
                         peer_factors[, 1:n.peer] + # peer factors
                         geno.int[, SNP]*timepoint + # interaction term 
                         (1|individual) + (1|batch), # individual and batch 
                       REML=FALSE)
  }
  
  
  # Get row indices used in the model
  model_rows <- as.numeric(rownames(model.frame(model.test)))
  
  # Extract random intercepts for individual and batch (example)
  re_ind <- ranef(model.test)$individual
  re_batch <- ranef(model.test)$batch
  # make sure the order is correct
  re_ind_contrib <- re_ind[individual[model_rows], "(Intercept)"]
  re_batch_contrib <- re_batch[batch[model_rows], "(Intercept)"]
  # sum random effects
  random_effects <- re_ind_contrib + re_batch_contrib
  
  # fixed effects
  X <- model.matrix(model.test)
  betas <- fixef(model.test)
  # remove the terms we want to keep 
  keep_terms <- c("geno.int[, SNP]", 
                  "timepoint12_weeks", "timepoint20_weeks", "timepoint28_weeks", "timepoint36_weeks", 
                  "geno.int[, SNP]:timepoint20_weeks", "geno.int[, SNP]:timepoint28_weeks", "geno.int[, SNP]:timepoint36_weeks")
  keep_idx <- which(names(betas) %in% keep_terms)
  non_interest_idx <- setdiff(seq_along(betas), keep_idx)
  
  X_non_interest <- X[, non_interest_idx, drop=FALSE]
  fitted_non_interest <- X_non_interest %*% betas[non_interest_idx]
  
  X_keep <- X[, keep_idx, drop = FALSE]                # Subset design matrix to terms of interest
  interest_effect <- X_keep %*% betas[keep_idx]        # Multiply by their fixed-effect coefficients
  
  observed <- getME(model.test, "y")
  partial_resid <- observed - fitted_non_interest - random_effects
  
  # check this is working as expected
  # partial_resid = manually computed partial residuals
  # interest_effect = fitted values from terms of interest
  # resid(model.test) = model residuals
  # all.equal(partial_resid, (interest_effect + resid(model.test)), tolerance = 1e-6)
  
  # Subset all vectors on the same index
  df <- data.frame(
    geno = geno.int[model_rows, SNP],
    timepoint = timepoint[model_rows],
    #outcome = outcome[model_rows],
    #interaction = geno.int[idx, SNP] * sga[idx],  # recompute interaction here
    partial_resid = partial_resid
  )
  
  gene_name <- id_name$gtf.gene_name[id_name$gtf.gene_id == gene]
  
  pval = ifelse(interaction_results$BH_adj_pval[signal] < 0.005, 
                format(interaction_results$BH_adj_pval[signal], scientific = TRUE, digits = 2),
                round(interaction_results$BH_adj_pval[signal], digits = 2))
  
  half_dodge <- 0.06
  
  plt <- df %>%
    mutate(snp_base = as.numeric(geno),  # ensure 0,1,2
           #outcome = as.factor(outcome),
           timepoint = as.factor(timepoint),
           group_label = interaction(snp_base, timepoint),
           x_pos = case_when(
             timepoint == levels(timepoint)[1] ~ snp_base - 1.5 * half_dodge,
             timepoint == levels(timepoint)[2] ~ snp_base - 0.5 * half_dodge,
             timepoint == levels(timepoint)[3] ~ snp_base + 0.5 * half_dodge,
             timepoint == levels(timepoint)[4] ~ snp_base + 1.5 * half_dodge
           )) %>%
    ggplot(aes(x = x_pos, y = partial_resid)) +
    geom_jitter(aes(color = timepoint), width = 0.01, height = 0, size = 0.5, alpha = 1) +
    geom_boxplot(aes(x = x_pos, group = interaction(snp_base, timepoint), color = timepoint),
                 outlier.shape = NA, width = 0.25, alpha = 0.5) +
    geom_smooth(aes(color = timepoint), method = "lm", formula = y ~ x, se = FALSE) +
    scale_x_continuous(breaks = 0:2, labels = c("0", "1", "2")) + 
    scale_color_manual(values=c("#BCE4D8", "#83C4CB", "#439FB7", "#32769B")) + 
    xlab(paste("N", SNP, "minor alleles")) + 
    ylab(paste(gene_name, "expression")) +
    labs(color = NULL) + 
    annotate("text", x = -0.4, y = max(df$partial_resid), 
             label = paste0("βG: ", 
                            round(interaction_results$eQTL_beta_baseline[signal], digits = 2),#interaction_results$eQTL_beta[signal], digits = 2),
                            "\nβGx20: ",
                            round(interaction_results$eqtl_diff_wk20_interaction[signal], digits = 2),#interaction_results$Interaction_wk20_beta[signal], digits = 2),
                            "\nβGx28: ",
                            round(interaction_results$eqtl_diff_wk28_interaction[signal], digits = 2),#interaction_results$Interaction_wk28_beta[signal], digits = 2),
                            "\nβGx36: ",
                            round(interaction_results$eqtl_diff_wk36_interaction[signal], digits = 2),#interaction_results$Interaction_wk36_beta[signal], digits = 2),
                            "\npint: ",
                            pval), 
             hjust = 0, vjust = 1) + 
    ggtitle("Time-point interaction") + 
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1))
  
  return(plt)
}

# get most significant interactions per group 
top_gene_snp_pairs_per_clust <- time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)) %>%
  left_join(time_int_clust %>% select(gsp, cluster_label)) %>%
  left_join(id_name, by = c("Gene" = "gtf.gene_id")) %>% 
  filter(eQTL_beta_baseline > 0) %>%
  group_by(cluster_label) %>%
  top_n(n = -6, BH_adj_pval) %>%
  select(Gene, SNP, BH_adj_pval, cluster_label) %>%
  mutate(rank = rank(BH_adj_pval)) %>%
  arrange(cluster_label, BH_adj_pval) %>%
  filter(ifelse(cluster_label %in% c(2, 5), rank == 1, 
                ifelse(cluster_label == 1, rank == 2, 
                       ifelse(cluster_label == 3, rank == 6, 
                              ifelse(cluster_label == 4, rank == 3)))))

make_time_int_plts <- function(cluster_ibeta_plt_df, cluster_label_in, mean_per_clust, 
                               n.peer = 30, signal,
                               id_name = id_name, 
                               interaction_results = time_interactions){
  # clust 1 
  cluster_traj <- cluster_ibeta_plt_df %>%
    filter(cluster_label == cluster_label_in) %>%
    ggplot(aes(x = timepoint, y = scaled_beta_sum, group = gsp)) + 
    geom_line(color = "#598eca", alpha = 1, linewidth = 0.05) + 
    geom_line(
      data = mean_per_clust %>% filter(cluster_label == cluster_label_in), 
      aes(x = timepoint, y = cluster_tp_mean, group = cluster_label), 
      color = "black", 
      #linetype = "dashed",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) + 
    geom_point(
      data = mean_per_clust %>% filter(cluster_label == cluster_label_in), 
      aes(x = timepoint, y = cluster_tp_mean, group = cluster_label), 
      color = "black"
    ) + 
    facet_wrap(~cluster_label, nrow = 2) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 13)) + 
    labs(#title = "Scaled interaction eQTL betas",
      x = "Sample time-point (weeks)",
      y = "Scaled magnitude of effect") + 
    scale_x_discrete(expand = c(0.05, 0)) + 
    scale_y_continuous(breaks = c(1, 0, -1))
  
  time_main_clust_plt <- plot_main_effect(n.peer = 30, 
                                          signal = signal, 
                                          id_name = id_name, 
                                          interaction_results = time_interactions)
  
  time_int_clust_plt <- plot_interaction_time(n.peer = 30, 
                                              signal = signal, 
                                              id_name = id_name, 
                                              interaction_results = time_interactions)
  return(list(cluster_traj, time_main_clust_plt, time_int_clust_plt))
}

plt_list <- list()
for(i in 1:5){
  plts_i <- make_time_int_plts(cluster_ibeta_plt_df, cluster_label_in = levels(cluster_ibeta_plt_df$cluster_label)[i], 
                               mean_per_clust, 
                               n.peer = 30, 
                               signal =  which(
                                 time_interactions$Gene == top_gene_snp_pairs_per_clust$Gene[top_gene_snp_pairs_per_clust$cluster_label == i] & 
                                   time_interactions$SNP == top_gene_snp_pairs_per_clust$SNP[top_gene_snp_pairs_per_clust$cluster_label == i]
                               ),
                               id_name = id_name, 
                               interaction_results = time_interactions)
  
  plt_list[[i]] <- plts_i
}

# Main figures 6A, B
all_time_int_clusters <- gridExtra::grid.arrange(plt_list[[1]][[1]] + labs(tag = "A") + xlab(""), 
                                                 plt_list[[2]][[1]] + labs(tag = " ") + xlab(""), 
                                                 plt_list[[3]][[1]] + labs(tag = " "), 
                                                 plt_list[[4]][[1]] + labs(tag = " ") + xlab(""), 
                                                 plt_list[[5]][[1]] + labs(tag = " ") + xlab(""), 
                                                 plt_list[[1]][[3]] + labs(tag = "B") + theme(legend.position="none"),
                                                 plt_list[[2]][[3]] + labs(tag = " ") + theme(legend.position="none"),
                                                 plt_list[[3]][[3]] + labs(tag = " ") + theme(legend.position="none"),
                                                 plt_list[[4]][[3]] + labs(tag = " ") + theme(legend.position="none"),
                                                 plt_list[[5]][[3]] + labs(tag = " ") + theme(legend.position="none"),
                                                 ncol = 5, nrow = 2, 
                                                 #widths = c(1, 1, 1.2), 
                                                 heights = c(1, 1),
                                                 layout_matrix = rbind(c(1, 2, 3, 4, 5),
                                                                       c(6, 7, 8, 9, 10)))

                                                 
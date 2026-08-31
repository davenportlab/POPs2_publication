# 07_7_get_signif_interactions.R

################################################################################

# 7.7. Run multiple testing correction on interaction eQTL 

################################################################################

# Aim: Identify significant eQTL interactions 

########################### Paths ##########################
path <- "genotyping/analysis/eQTL/output_data/"

########################### Load packages ###########################
library(tidyverse)
library(ggupset)

########################### Load data ###########################
# read in time interactions
time_interactions <- readRDS(paste0(path, "interactions_flanders_2/tp_int_results_flanders_2.rds")) %>%
  mutate(eQTL_beta_baseline = as.numeric(eQTL_beta_baseline), 
         eqtl_diff_wk20_interaction = as.numeric(eqtl_diff_wk20_interaction),
         eqtl_diff_wk28_interaction = as.numeric(eqtl_diff_wk28_interaction),
         eqtl_diff_wk36_interaction = as.numeric(eqtl_diff_wk36_interaction), 
         tp1_nminoralleles = as.numeric(tp1_nminoralleles), 
         tp2_nminoralleles = as.numeric(tp2_nminoralleles), 
         tp3_nminoralleles = as.numeric(tp3_nminoralleles), 
         tp4_nminoralleles = as.numeric(tp4_nminoralleles))
# read in neutrophil interactions
neutrophil_interactions <- readRDS(paste0(path, "interactions_flanders_2/neutrophil_int_results_flanders_2.rds")) %>%
  mutate(bottom_half_nminoralleles = as.numeric(bottom_half_nminoralleles), 
         top_half_nminoralleles = as.numeric(top_half_nminoralleles)) 
# read in monocyte interactions
monocyte_interactions <- readRDS(paste0(path, "interactions_flanders_2/monocyte_int_results_flanders_2.rds"))
# read in neutrophil time interactions
neutrophil_time_interactions <- readRDS(paste0(path, "interactions_flanders_2/neutrophil_int_and_time_results_flanders_2.rds")) %>%
  mutate(bottom_half_nminoralleles = as.numeric(bottom_half_nminoralleles), 
         top_half_nminoralleles = as.numeric(top_half_nminoralleles),
         tp1_nminoralleles = as.numeric(tp1_nminoralleles),
         tp2_nminoralleles = as.numeric(tp2_nminoralleles),
         tp3_nminoralleles = as.numeric(tp3_nminoralleles),
         tp4_nminoralleles = as.numeric(tp4_nminoralleles))
# read in monocyte time interactions 
monocyte_time_interactions <- readRDS(paste0(path, "interactions_flanders_2/monocyte_int_and_time_results_flanders_2.rds")) %>%
  mutate(bottom_half_nminoralleles = as.numeric(bottom_half_nminoralleles), 
         top_half_nminoralleles = as.numeric(top_half_nminoralleles),
         tp1_nminoralleles = as.numeric(tp1_nminoralleles),
         tp2_nminoralleles = as.numeric(tp2_nminoralleles),
         tp3_nminoralleles = as.numeric(tp3_nminoralleles),
         tp4_nminoralleles = as.numeric(tp4_nminoralleles)) 

# read in conditional results with sign flips for time clustering 
all_conditional_signals_aligned <- readRDS(paste0(path, "flanders_conditional_eQTL/all_conditional_signals_sign_aligned_pops.rds"))

########################### Analysis ###########################

################################################################################
############## multiple testing correction: time interactions ##################
################################################################################
# what proportion were tested? 
time_interactions %>% filter(reason_not_tested == "tested", 
                             tp1_nminoralleles > 1, 
                             tp2_nminoralleles > 1, 
                             tp3_nminoralleles > 1, 
                             tp4_nminoralleles > 1) %>% nrow() / time_interactions %>% nrow()
# 0.8304204

time_interactions <- time_interactions %>% 
  # subset only to gene SNP pairs that were tested for interaction (for multiple testing correction)
  filter(reason_not_tested == "tested", 
         tp1_nminoralleles > 1, 
         tp2_nminoralleles > 1, 
         tp3_nminoralleles > 1, 
         tp4_nminoralleles > 1) %>%
  # make a bh adjusted pvalue 
  mutate(BH_adj_pval = p.adjust(Interaction_pval, method = "BH")) %>%
  # filter to significant bh adjusted pvalues
  filter(BH_adj_pval < 0.05) 

time_interactions %>%
  nrow()
# 1193

# write this out 
time_interactions %>% write.csv(paste0(path,"interactions_flanders_2/time_interactions_flanders_2.csv"), row.names = FALSE, quote = FALSE)

# make df for input into time interaction clustering 
time_interactions_forclust <- time_interactions %>%
  # calculate the beta values at each time-point
  mutate(beta_wk12 = as.numeric(eQTL_beta_baseline), 
         beta_wk20 = beta_wk12 + as.numeric(eqtl_diff_wk20_interaction),
         beta_wk28 = beta_wk12 + as.numeric(eqtl_diff_wk28_interaction),
         beta_wk36 = beta_wk12 + as.numeric(eqtl_diff_wk36_interaction)) %>%
  left_join(all_conditional_signals_aligned %>%
              select(snps, phenotype_id, bC_aligned), 
            by = c("SNP" = "snps", "Gene" = "phenotype_id")) %>%
  select(SNP, Gene, beta_wk12, beta_wk20, beta_wk28, beta_wk36, bC_aligned) %>%
  mutate(sign_beta = ifelse(as.numeric(bC_aligned) < 0, -1, 1)) %>%
  rowwise() %>%
  mutate(
    n_negative = sum(c_across(
      c(beta_wk12,
        beta_wk20,
        beta_wk28,
        beta_wk36)
    ) < 0)
  ) %>%
  ungroup()

time_interactions_forclust %>%
  mutate(gsp = paste0(Gene, "_", SNP), 
         wk12 = beta_wk12*sign_beta, 
         wk20 = beta_wk20*sign_beta,
         wk28 = beta_wk28*sign_beta, 
         wk36 = beta_wk36*sign_beta) %>%
  select(gsp, wk12, wk20, wk28, wk36) %>%
  write.csv(paste0(path,"interactions_flanders_2/time_interactions_flanders_forclust_2.csv"), row.names = FALSE, quote = FALSE)

################################################################################
########### multiple testing correction: neutrophil interactions ###############
################################################################################
# what proportion were tested? 
neutrophil_interactions %>% filter(reason_not_tested == "tested", 
                                   bottom_half_nminoralleles > 1, 
                                   top_half_nminoralleles > 1) %>% nrow() / neutrophil_interactions %>% nrow()
# 0.8395155 

neutrophil_interactions <- neutrophil_interactions %>% 
  filter(reason_not_tested == "tested", 
         bottom_half_nminoralleles > 1, 
         top_half_nminoralleles > 1) %>%
  mutate(BH_adj_pval = p.adjust(Interaction_pval, method = "BH")) %>%
  filter(BH_adj_pval < 0.05) 

neutrophil_interactions %>%
  nrow()
# [1] 2154 

# write this out 
neutrophil_interactions %>% write.csv(paste0(path,"interactions_flanders_2/neutrophil_interactions_flanders_2.csv"), row.names = FALSE, quote = FALSE)

################################################################################
########### multiple testing correction: monocyte interactions #################
################################################################################
# what proportion were tested? 
monocyte_interactions %>% filter(reason_not_tested == "tested", 
                                 bottom_half_nminoralleles > 1, 
                                 top_half_nminoralleles > 1) %>% nrow() / monocyte_interactions %>% nrow()
# 0.8399346

monocyte_interactions <- monocyte_interactions %>% 
  filter(reason_not_tested == "tested", 
         bottom_half_nminoralleles > 1, 
         top_half_nminoralleles > 1) %>%
  mutate(BH_adj_pval = p.adjust(Interaction_pval, method = "BH")) %>%
  filter(BH_adj_pval < 0.05) 

monocyte_interactions %>%
  nrow()
# 953

# write this out 
monocyte_interactions %>% write.csv(paste0(path,"interactions_flanders_2/monocyte_interactions_flanders_2.csv"), row.names = FALSE, quote = FALSE)

################################################################################
########### multiple testing correction: neutrophil time interactions ##########
################################################################################
# what proportion were tested? 
neutrophil_time_interactions %>% filter(!is.na(reason_not_tested), 
                                        bottom_half_nminoralleles > 1, 
                                        top_half_nminoralleles > 1, 
                                        tp1_nminoralleles > 1, 
                                        tp2_nminoralleles > 1, 
                                        tp3_nminoralleles > 1, 
                                        tp4_nminoralleles > 1) %>% nrow() / neutrophil_time_interactions %>% nrow()
# 0.8232952

neutrophil_time_interactions <- neutrophil_time_interactions %>% 
  filter(!is.na(reason_not_tested), 
         bottom_half_nminoralleles > 1, 
         top_half_nminoralleles > 1, 
         tp1_nminoralleles > 1, 
         tp2_nminoralleles > 1, 
         tp3_nminoralleles > 1, 
         tp4_nminoralleles > 1) %>%
  mutate(BH_adj_pval = p.adjust(Interaction_pval, method = "BH")) %>%
  filter(BH_adj_pval < 0.05) 

neutrophil_time_interactions %>%
  nrow()
# [1] 774

# write this out 
neutrophil_time_interactions %>% write.csv(paste0(path,"interactions_flanders_2/neutrophil_time_interactions_flanders_2.csv"), row.names = FALSE, quote = FALSE)

################################################################################
########### multiple testing correction: monocyte time interactions ############
################################################################################
# what proportion were tested? 
monocyte_time_interactions %>% filter(!is.na(reason_not_tested), 
                                      bottom_half_nminoralleles > 1, 
                                      top_half_nminoralleles > 1, 
                                      tp1_nminoralleles > 1, 
                                      tp2_nminoralleles > 1, 
                                      tp3_nminoralleles > 1, 
                                      tp4_nminoralleles > 1) %>% nrow() / monocyte_time_interactions %>% nrow()
# 0.8228761

monocyte_time_interactions <- monocyte_time_interactions %>% 
  filter(!is.na(reason_not_tested), 
         bottom_half_nminoralleles > 1, 
         top_half_nminoralleles > 1, 
         tp1_nminoralleles > 1, 
         tp2_nminoralleles > 1, 
         tp3_nminoralleles > 1, 
         tp4_nminoralleles > 1) %>%
  mutate(BH_adj_pval = p.adjust(Interaction_pval, method = "BH")) %>%
  filter(BH_adj_pval < 0.05) 

monocyte_time_interactions %>%
  nrow()
# [1] 785

# write this out 
monocyte_time_interactions %>% write.csv(paste0(path,"interactions_flanders_2/monocyte_time_interactions_flanders_2.csv"), row.names = FALSE, quote = FALSE)

################################################################################
########### compare overlapping time and cellprop interactions #################
################################################################################

# load significant interactions 
time_interactions <- read.csv(paste0(path,"interactions_flanders_2/time_interactions_flanders_2.csv"))
neutrophil_interactions <- read.csv(paste0(path,"interactions_flanders_2/neutrophil_interactions_flanders_2.csv"))
monocyte_interactions <- read.csv(paste0(path,"interactions_flanders_2/monocyte_interactions_flanders_2.csv"))
monocyte_time_interactions <- read.csv(paste0(path,"interactions_flanders_2/monocyte_time_interactions_flanders_2.csv"))
neutrophil_time_interactions <- read.csv(paste0(path,"interactions_flanders_2/neutrophil_time_interactions_flanders_2.csv"))


listInput_gsp <- list(
  `Time-point` = (time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Neutrophils = (neutrophil_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp,
  Monocytes = (monocyte_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp,
  Monocytes_time = (monocyte_time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp,
  Neutrophils_time = (neutrophil_time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp)

# Convert list to long tibble
gsp_long <- enframe(listInput_gsp, name = "interaction", value = "gene_snp_pair") %>%
  unnest(gene_snp_pair)

# Combine by gene to form sets per gene
gsp_sets <- gsp_long %>%
  group_by(gene_snp_pair) %>%
  summarise(interactions = list(interaction), .groups = "drop")

# number of time interactions that overlap neutrophil and monocyte interactions
gsp_long %>%
  filter(!(interaction %in% c("Neutrophils_time", "Monocytes_time"))) %>%
  distinct(interaction, gene_snp_pair) %>%   # avoid duplicates
  group_by(gene_snp_pair) %>%
  summarise(
    combination = paste(sort(interaction), collapse = " & "),
    .groups = "drop"
  ) %>%
  count(combination, sort = TRUE) %>%
  filter(grepl("Time-point", combination)) %>%
  mutate(
    percent = n / sum(n),
    percent_label = scales::percent(percent, accuracy = 0.1)
  )
# combination                              n percent percent_label
# <chr>                                   <int>   <dbl> <chr>        
# 1 Time-point                             483  0.405  40.5%        
# 2 Monocytes & Neutrophils & Time-point   424  0.355  35.5%        
# 3 Neutrophils & Time-point               223  0.187  18.7%        
# 4 Monocytes & Time-point                  63  0.0528 5.3%     

# Time-point interactions that don't overlap neutrophils or monocytes: 40.5% (n = 483)
# Time-point interactions that DO overlap neutrophils or monocytes: 59.5% (n = 710)

# overlap of neutrophil and time interactions
gsp_long %>%
  filter(!(interaction %in% c("Monocytes", "Monocytes_time"))) %>%
  distinct(interaction, gene_snp_pair) %>%   # avoid duplicates
  group_by(gene_snp_pair) %>%
  summarise(
    combination = paste(sort(interaction), collapse = " & "),
    .groups = "drop"
  ) %>%
  count(combination, sort = TRUE) %>%
  filter(grepl("Time-point", combination)) %>%
  mutate(
    percent = n / sum(n),
    percent_label = scales::percent(percent, accuracy = 0.1)
  )
# # A tibble: 4 × 4
# combination                                     n percent percent_label
# <chr>                                           <int>   <dbl> <chr>        
# 1 Neutrophils_time & Time-point                 405   0.339 33.9%        
# 2 Neutrophils & Neutrophils_time & Time-point   328   0.275 27.5%        
# 3 Neutrophils & Time-point                      319   0.267 26.7%        
# 4 Time-point                                    141   0.118 11.8%   

# of eQTL with both neutrophil and time interactions, 
# ~50% have evidence of additional time effects (27.5%) and ~50% don't (26.7%)

# overlap of neutrophil and time interactions
gsp_long %>%
  filter(!(interaction %in% c("Neutrophils", "Neutrophils_time"))) %>%
  distinct(interaction, gene_snp_pair) %>%   # avoid duplicates
  group_by(gene_snp_pair) %>%
  summarise(
    combination = paste(sort(interaction), collapse = " & "),
    .groups = "drop"
  ) %>%
  count(combination, sort = TRUE) %>%
  filter(grepl("Time-point", combination)) %>%
  mutate(
    percent = n / sum(n),
    percent_label = scales::percent(percent, accuracy = 0.1)
  )

# # A tibble: 4 × 4
# combination                                 n percent percent_label
# <chr>                                     <int>   <dbl> <chr>        
# 1 Monocytes_time & Time-point               460   0.386 38.6%     # only time effects   
# 2 Monocytes & Monocytes_time & Time-point   280   0.235 23.5%     # time and additional monocyte effects
# 3 Time-point                                246   0.206 20.6%     # only time effects   
# 4 Monocytes & Time-point                    207   0.174 17.4%     # separate time and monocyte effects   

# of eQTL with both monocyte and time interactions, 
# ~60% have evidence of additional time effects (23.5%) and ~40% don't (17.4%)


########################### Plot ###########################
# time and cell prop sep
listInput_gsp_time_int_sep <- list(
  `Time-point` = (time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Neutrophils = (neutrophil_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp,
  Monocytes = (monocyte_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp
)


# Convert list to long tibble
gsp_long_time_int_sep <- enframe(listInput_gsp_time_int_sep, name = "interaction", value = "gene_snp_pair") %>%
  unnest(gene_snp_pair)

# Combine by gene to form sets per gene
gsp_sets_time_int_sep <- gsp_long_time_int_sep %>%
  group_by(gene_snp_pair) %>%
  summarise(interactions = list(interaction), .groups = "drop")

# Count frequencies for labeling
counts_time_int_sep <- gsp_sets_time_int_sep %>%
  count(interactions)

int_gsp_upset_plt_time_int_sep <- ggplot(gsp_sets_time_int_sep, aes(x = interactions)) +
  geom_bar(color = "black", fill = "black") +
  geom_text(
    data = counts_time_int_sep,
    aes(x = interactions, y = n, label = n),
    vjust = -0.3, size = 3
  ) +
  scale_x_upset(order_by = "freq") +
  labs(
    x = "Interaction intersections",
    y = "n SNP-gene pairs"
  ) + 
  ggtitle("Overlap of eQTL interactions with time and cell proportions") + theme(
    panel.border = element_blank(),       # remove panel border
    panel.background = element_blank(),   # remove panel background
    panel.grid = element_blank(),         # remove grid lines
    axis.line = element_line(color = "white")  # optional: add axis lines
  )

# monocyte only
listInput_gsp_mon_only <- list(
  `Time-point` = (time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Monocytes = (monocyte_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Monocytes_time = (monocyte_time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp)

# Convert list to long tibble
gsp_long_mon_only <- enframe(listInput_gsp_mon_only, name = "interaction", value = "gene_snp_pair") %>%
  unnest(gene_snp_pair)

# Combine by gene to form sets per gene
gsp_sets_mon_only <- gsp_long_mon_only %>%
  group_by(gene_snp_pair) %>%
  summarise(interactions = list(interaction), .groups = "drop")

# Count frequencies for labeling
counts_mon_only <- gsp_sets_mon_only %>%
  count(interactions)

# Plot overlap between time-point, neutrophil and monocyte interactions
# (Supplementary Figure 19A)
int_gsp_upset_plt_mon_only <- ggplot(gsp_sets_mon_only, aes(x = interactions)) +
  geom_bar(color = "black", fill = "black") +
  geom_text(
    data = counts_mon_only,
    aes(x = interactions, y = n, label = n),
    vjust = -0.3, size = 3
  ) +
  scale_x_upset(order_by = "freq") +
  labs(
    x = "Interaction intersections",
    y = "n SNP-gene pairs"
  ) + 
  ggtitle("Overlap of eQTL interactions with time and monocytes") + theme(
    panel.border = element_blank(),       # remove panel border
    panel.background = element_blank(),   # remove panel background
    panel.grid = element_blank(),         # remove grid lines
    axis.line = element_line(color = "white")  # optional: add axis lines
  )

# monocyte only
listInput_gsp_mon_only <- list(
  `Time-point` = (time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Monocytes = (monocyte_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Monocytes_time = (monocyte_time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp)

# Convert list to long tibble
gsp_long_mon_only <- enframe(listInput_gsp_mon_only, name = "interaction", value = "gene_snp_pair") %>%
  unnest(gene_snp_pair)

# Combine by gene to form sets per gene
gsp_sets_mon_only <- gsp_long_mon_only %>%
  group_by(gene_snp_pair) %>%
  summarise(interactions = list(interaction), .groups = "drop")

# Count frequencies for labeling
counts_mon_only <- gsp_sets_mon_only %>%
  count(interactions)

# Plot overlaps between time-point, monocyte and monocyte-time interactions
# (Supplementary Figure 19C)
int_gsp_upset_plt_mon_only <- ggplot(gsp_sets_mon_only, aes(x = interactions)) +
  geom_bar(color = "black", fill = "black") +
  geom_text(
    data = counts_mon_only,
    aes(x = interactions, y = n, label = n),
    vjust = -0.3, size = 3
  ) +
  scale_x_upset(order_by = "freq") +
  labs(
    x = "Interaction intersections",
    y = "n SNP-gene pairs"
  ) + 
  ggtitle("Overlap of eQTL interactions with time and monocytes") + theme(
    panel.border = element_blank(),       # remove panel border
    panel.background = element_blank(),   # remove panel background
    panel.grid = element_blank(),         # remove grid lines
    axis.line = element_line(color = "white")  # optional: add axis lines
  )

# neutrophil only
listInput_gsp_neu_only <- list(
  `Time-point` = (time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Neutrophils = (neutrophil_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp, 
  Neutrophils_time = (neutrophil_time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)))$gsp)

# Convert list to long tibble
gsp_long_neu_only <- enframe(listInput_gsp_neu_only, name = "interaction", value = "gene_snp_pair") %>%
  unnest(gene_snp_pair)

# Combine by gene to form sets per gene
gsp_sets_neu_only <- gsp_long_neu_only %>%
  group_by(gene_snp_pair) %>%
  summarise(interactions = list(interaction), .groups = "drop")

# Count frequencies for labeling
counts_neu_only <- gsp_sets_neu_only %>%
  count(interactions)

# Plot overlaps between time-point, neutrophil and neutrophil-time interactions
# (Supplementary Figure 19B)
int_gsp_upset_plt_neu_only <- ggplot(gsp_sets_neu_only, aes(x = interactions)) +
  geom_bar(color = "black", fill = "black") +
  geom_text(
    data = counts_neu_only,
    aes(x = interactions, y = n, label = n),
    vjust = -0.3, size = 3
  ) +
  scale_x_upset(order_by = "freq") +
  labs(
    x = "Interaction intersections",
    y = "n SNP-gene pairs"
  ) + 
  ggtitle("Overlap of eQTL interactions with time and neutrophils") + theme(
    panel.border = element_blank(),       # remove panel border
    panel.background = element_blank(),   # remove panel background
    panel.grid = element_blank(),         # remove grid lines
    axis.line = element_line(color = "white")  # optional: add axis lines
  )
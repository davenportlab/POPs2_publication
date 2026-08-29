# 03_3_compile_conditional_flanders_signals.R

################################################################################

# 3.3. Compile conditional results

################################################################################

# Aim: generate log2cpm count matrix for following analyses (including gene filtering)

########################### Output paths ##########################
dir <- "genotyping/analysis/eQTL/"

########################### Input paths ###########################

########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)

########################### Load data ###########################

########################### Analysis ###########################

full_cond_results <- readRDS(paste0(dir, "output_data/flanders_conditional_eQTL/conditional_eQTL_results.rds")) %>%
  # Filter to unique rows (to get rid of one problem gene-snp pair)
  distinct()

signals_for_conditional_mapping <- readRDS(paste0(dir, "input_data/conditional_eQTL/gene_snp_pairs_for_cond_mapping.rds")) %>%
  arrange(chr, phenotype_id, signal_rank)

# investigate - how often am I losing the most significant signal 
lost_primary <- signals_for_conditional_mapping %>% 
  left_join(full_cond_results %>% select(Gene, SNP, conditionally_independent), 
            by = c("snps" = "SNP", "phenotype_id" = "Gene")) %>% 
  filter(conditionally_independent == FALSE, 
         signal_rank == 1) 
lost_primary %>%
  nrow()
# [1] 124

# check if there are any genes where no signals are conditionally independent 
no_cond_ind <- full_cond_results %>%
  group_by(Gene) %>%
  summarise(
    any_ci = any(conditionally_independent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!any_ci) 

no_cond_ind %>%
  nrow()
# [1] 100

# how many of these overlap? 
length(intersect(lost_primary$phenotype_id, no_cond_ind$Gene))
# [1] 100
# 100 of the 124 will be rescued

# select original flanders signals based on the following
# IF at least one signal in a gene is conditionally independent: select all conditionally independent signals
# ELSE select the most significant signal in the gene from the flanders output 
signals_selected <- signals_for_conditional_mapping %>%
  # filter out the problem gene snp pair second signal
  filter(!((phenotype_id == "ENSG00000143851" & 
              snps == "rs2924110") & 
             signal_rank != 1)) %>%
  # bring in conditional independence info
  left_join(
    full_cond_results %>%
      select(Gene, SNP, conditionally_independent),
    by = c("phenotype_id" = "Gene", "snps" = "SNP")
  ) %>%
  group_by(phenotype_id) %>%
  # decision tree
  filter(
    if (any(conditionally_independent, na.rm = TRUE)) {
      # in cases of complete LD between eSNPs, the model controling for one and testing the other is rank deficient so one of the SNPs is dropped. In this case the anova testing for a difference between models does not give a pvalue, so it is NA
      conditionally_independent
    } else {
      top_pvalue == min(top_pvalue)
    }
  ) %>%
  # re-rank within gene
  mutate(
    signal_rank = rank(top_pvalue, ties.method = "first")
  ) %>%
  ungroup() %>%
  arrange(chr, phenotype_id, signal_rank)

dim(signals_for_conditional_mapping)
#[1] 5638   17
dim(signals_selected)
# [1] 3718   18

# next, join this to the list of signals that don't need conditional mapping 
# read in signals that didn't need conditional mapping
sigs_no_cond_mapping <- readRDS(paste0(dir, "input_data/conditional_eQTL/gene_snp_pairs_single_locus_no_cond_mapping.rds"))

all_conditional_signals <- rbind(sigs_no_cond_mapping, 
                                 signals_selected %>% select(-conditionally_independent)) %>%
  select(-i.chr) %>%
  arrange(chr, phenotype_id, signal_rank)

# how many genes have conditional signals (should be all that were finemapped)
all_conditional_signals %>% select(phenotype_id) %>% distinct() %>% nrow()
# [1] 14268

# write this out as the final list of conditional signals 
saveRDS(all_conditional_signals, 
        paste0(dir, "output_data/flanders_conditional_eQTL/all_conditional_signals_pops.rds"))

# then prepare for interaction mapping
# swap signs for the beta values based on the minor/major alleles + write this out 
all_conditional_signals_aligned = all_conditional_signals %>% 
  mutate(bC_aligned = ifelse(A1 != a1, -bC, bC))

# write this out 
saveRDS(all_conditional_signals_aligned, 
        paste0(dir, "output_data/flanders_conditional_eQTL/all_conditional_signals_sign_aligned_pops.rds"))
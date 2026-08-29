# 03_1_prep_finemapped_sigs_for_conditional_mapping.R

################################################################################

# 3.1. Identification of flanders signals for conditional analysis

################################################################################

# Aim: Identify flanders signals where more than one signal was discovered in 
#     more than one locus per gene

########################### Output paths ##########################
path <- "genotyping/analysis/eQTL/"

########################### Input paths ###########################

########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)
library(data.table)

########################### Load data ###########################
finemapped.loci <- readRDS(paste0(path, "output_data/flanders_finemap/finemapped_loci_pops.rds"))

########################### Analysis ###########################


gene_locus_summary <- finemapped.loci %>% 
  # group by gene
  group_by(phenotype_id) %>%
  summarise(n_unique_loci = n_distinct(path_rds))

##############################################################
# Signals where conditional independence needs testing 
##############################################################

# how many genes have more than one locus where signals were discovered? 
gene_locus_summary %>% filter(n_unique_loci > 1) %>% nrow()

# get subset of genes where conditional independence will need to be tested
genes_w_multi_loci <- gene_locus_summary %>% 
  filter(n_unique_loci > 1)
dim(genes_w_multi_loci)

# get a list of the signals for each gene that does require conditional mapping 
genes_for_cond_mapping <- finemapped.loci %>% 
  filter(phenotype_id %in% genes_w_multi_loci$phenotype_id) %>%
  mutate(locus = path_rds) %>% #gsub(paste0(path, "output_data/flanders_4/flanders_output/results/finemap/"), "", path_rds)) %>%
  select(phenotype_id, chr, start, end, top_pvalue, locus, snp, a1, a0, bC, bC_se) %>%
  # add a signal rank based on pvalue per gene 
  group_by(phenotype_id) %>%
  mutate(signal_rank = rank(top_pvalue, ties.method = "first")) %>%
  ungroup()

# write out list of genes
unique(genes_for_cond_mapping$phenotype_id) %>% write.csv(paste0(path, "input_data/conditional_eQTL/genes_for_cond_mapping.csv"), row.names = FALSE, quote = FALSE)

# read in original SNP IDs
snp_id_key <- read.table(paste0(path, "input_data/flanders/POPS2_main_effect_input/POPS2_flanders_main_effects.txt"),
                         header = TRUE)

# add a flanders format snp id
setDT(snp_id_key)  # make sure it’s a data.table

# reorder A1 and A2 alphabetically and remake flanders_snp
snp_id_key[, `:=`(
  # reorder alleles
  allele1 = pmin(A1, A2),
  allele2 = pmax(A1, A2),
  # remake flanders_snp with new order
  flanders_snp = paste0("chr", chr, ":", SNPpos, ":", pmin(A1, A2), ":", pmax(A1, A2))
)]

# write this out
snp_id_key %>% select(snps, gene, chr, SNPpos, A1, A2, flanders_snp) %>%
  saveRDS(paste0(path, "input_data/conditional_eQTL/flanders_snp_rsid_key.rds"))

snp_id_key_in <- readRDS(paste0(path, "input_data/conditional_eQTL/flanders_snp_rsid_key.rds"))

# convert to data.table if not already
setDT(genes_for_cond_mapping)
setDT(snp_id_key_in)

# remove 'gene' and take unique rows
# select only the columns you need and unique
snp_key_unique <- unique(snp_id_key_in[, .(snps, chr, SNPpos, A1, A2, flanders_snp)])

# perform the left join
cond_mapping_signals_out <- snp_key_unique[genes_for_cond_mapping, 
                                           on = c("flanders_snp" = "snp")]

saveRDS(cond_mapping_signals_out, 
        paste0(path, "input_data/conditional_eQTL/gene_snp_pairs_for_cond_mapping.rds"))

##############################################################
# Signals where conditional independence does not need testing 
##############################################################

# get subset of genes where conditional independence will not need to be tested 
# because all signals were detected in the same locus 
genes_w_single_locus <- gene_locus_summary %>%
  filter(n_unique_loci == 1)
dim(genes_w_single_locus)

# get a list of signals for each genes that don't need conditional mapping
genes_no_cond_mapping <- finemapped.loci %>% 
  filter(phenotype_id %in% genes_w_single_locus$phenotype_id) %>%
  mutate(locus = path_rds) %>% 
  select(phenotype_id, chr, start, end, top_pvalue, locus, snp, a1, a0, bC, bC_se) %>%
  # add a signal rank based on pvalue per gene 
  group_by(phenotype_id) %>%
  mutate(signal_rank = rank(top_pvalue, ties.method = "first")) %>%
  ungroup() %>%
  arrange(chr, phenotype_id, signal_rank)

setDT(genes_no_cond_mapping)

# perform the left join
no_cond_mapping_signals_out <- snp_key_unique[genes_no_cond_mapping, 
                                              on = c("flanders_snp" = "snp")]

saveRDS(no_cond_mapping_signals_out, paste0(path, "input_data/conditional_eQTL/gene_snp_pairs_single_locus_no_cond_mapping.rds"))

# 06_1_direct_interval_comparison.R

################################################################################

# 6.1. Comparison of POPs2 TP eQTL to Interval

################################################################################

# Aim: format data for comparison between the five eQTL datasets

########################### Output paths ##########################

########################### Input paths ###########################
dir <-  "genotyping/analysis/eQTL/"
interval_dir <- "interval_rna/POPS2_comparison/"

########################### Parameters ############################
options(stringsAsFactors = FALSE)

########################### Load packages ###########################
library(data.table)
library(ggplot2)

########################### Load data ###########################

########################### Analysis ###########################

pops.bim <- data.frame(fread(paste0(dir, "input_data/single_tp_eQTL/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps.bim"),
                             header=FALSE))
# POPS lead SNP-gene pairs
pops_12wk <- read.delim(paste0(dir, "output_data/single_tp_analysis/12_weeks/eigenMT/12_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  mutate(snps = gsub("\\.", ":", snps)) %>%
  left_join(pops.bim %>% rename("snps" = "V2", "minor" = "V5", "major" = "V6") %>% select(snps, minor, major), by = "snps") %>%
  mutate(Pairs = paste0(snps, "_", gene), Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene)) 
pops_20wk <- read.delim(paste0(dir, "output_data/single_tp_analysis/20_weeks/eigenMT/20_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  mutate(snps = gsub("\\.", ":", snps)) %>%
  left_join(pops.bim %>% rename("snps" = "V2", "minor" = "V5", "major" = "V6") %>% select(snps, minor, major), by = "snps") %>%
  mutate(Pairs = paste0(snps, "_", gene), Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene))
pops_28wk <- read.delim(paste0(dir, "output_data/single_tp_analysis/28_weeks/eigenMT/28_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  mutate(snps = gsub("\\.", ":", snps)) %>%
  left_join(pops.bim %>% rename("snps" = "V2", "minor" = "V5", "major" = "V6") %>% select(snps, minor, major), by = "snps") %>%
  mutate(Pairs = paste0(snps, "_", gene), Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene))
pops_36wk <- read.delim(paste0(dir, "output_data/single_tp_analysis/36_weeks/eigenMT/36_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  mutate(snps = gsub("\\.", ":", snps)) %>%
  left_join(pops.bim %>% rename("snps" = "V2", "minor" = "V5", "major" = "V6") %>% select(snps, minor, major), by = "snps") %>%
  mutate(Pairs = paste0(snps, "_", gene), Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene))

# make a list of all SNP-gene pairs that are lead/unique
all_lead_pairs_notrsID <- unique(c(pops_12wk$Pairs_notrsID, pops_20wk$Pairs_notrsID, pops_28wk$Pairs_notrsID, pops_36wk$Pairs_notrsID))

#lead_snp_bim <- pops.bim %>% filter(V2 %in% all_lead_snps)
setDT(pops.bim)
# Rename and select relevant columns (faster than dplyr)
pops.bim_small <- pops.bim[, .(snps = V2, minor = V5, major = V6)]

# read in full results and subset to pairs that are lead in at least one time-point
# 12 weeks
pops_12wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/12_weeks/12_weeks_ciseqtl_all.rds")) 
setDT(pops_12wk_full)
pops_12wk_full <- pops.bim_small[pops_12wk_full, on = "snps"]
pops_12wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
pops_12wk_res_all <- pops_12wk_full[Pairs_notrsID %in% all_lead_pairs_notrsID]
# 20 weeks
pops_20wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/20_weeks/20_weeks_ciseqtl_all.rds")) 
setDT(pops_20wk_full)
pops_20wk_full <- pops.bim_small[pops_20wk_full, on = "snps"]
pops_20wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
pops_20wk_res_all <- pops_20wk_full[Pairs_notrsID %in% all_lead_pairs_notrsID]
# 28 weeks
pops_28wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/28_weeks/28_weeks_ciseqtl_all.rds")) 
setDT(pops_28wk_full)
pops_28wk_full <- pops.bim_small[pops_28wk_full, on = "snps"]
pops_28wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
pops_28wk_res_all <- pops_28wk_full[Pairs_notrsID %in% all_lead_pairs_notrsID]
# 36 weeks
pops_36wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/36_weeks/36_weeks_ciseqtl_all.rds")) 
setDT(pops_36wk_full)
pops_36wk_full <- pops.bim_small[pops_36wk_full, on = "snps"]
pops_36wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
pops_36wk_res_all <- pops_36wk_full[Pairs_notrsID %in% all_lead_pairs_notrsID]

# Load interval results
interval.results <- readRDS(paste0(interval_dir, "eQTL_main_effects/interval_ciseqtl_all.rds"))
setDT(interval.results)
# bim
interval.bim <- data.frame(fread(paste0(interval_dir, "data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_maf0.05_snps.bim"),
                                 header=FALSE))
setDT(interval.bim)
interval.bim_small <- interval.bim[, .(snps = V2, minor = V5, major = V6)]
interval.results <- interval.bim_small[interval.results, on = "snps"]
interval.results[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
interval.results_all <- interval.results[Pairs_notrsID %in% all_lead_pairs_notrsID]

write.table(interval.results_all, paste0(dir, "input_data/interval_comp/interval_snp_gene_pairs_in_pops_lead_byPos.txt"), sep="\t", quote=F)
interval.results_all <- interval.results_all %>% rename("interval_statistic" = "statistic", "interval_pvalue" = "pvalue", "interval_beta" = "beta", "interval_se" = "se", 
                                                        "interval_threshold" = "threshold", "interval_minor" = "minor", "interval_major" = "major", "interval_snps" = "snps")

# now merge the interval results with the single tp results
# merge these
# each will now have a different n pairs because each time-point could have a different lead SNP which may or may not be present in Interval...
comp_12wk <- merge(interval.results_all, 
                   pops_12wk_res_all %>% 
                     rename("pops_statistic" = "statistic", "pops_pvalue" = "pvalue", "pops_beta" = "beta", "pops_se" = "se", 
                            "pops_threshold" = "threshold", "pops_minor" = "minor", "pops_major" = "major", "pops_snps" = "snps"), 
                   by=c("Pairs_notrsID", "gene", "chr", "SNPpos", "TSS"))
comp_20wk <- merge(interval.results_all, pops_20wk_res_all %>%
                     rename("pops_statistic" = "statistic", "pops_pvalue" = "pvalue", "pops_beta" = "beta", "pops_se" = "se", 
                            "pops_threshold" = "threshold", "pops_minor" = "minor", "pops_major" = "major", "pops_snps" = "snps"), 
                   by=c("Pairs_notrsID", "gene", "chr", "SNPpos", "TSS"))
comp_28wk <- merge(interval.results_all, pops_28wk_res_all %>%
                     rename("pops_statistic" = "statistic", "pops_pvalue" = "pvalue", "pops_beta" = "beta", "pops_se" = "se", 
                            "pops_threshold" = "threshold", "pops_minor" = "minor", "pops_major" = "major", "pops_snps" = "snps"), 
                   by=c("Pairs_notrsID", "gene", "chr", "SNPpos", "TSS"))
comp_36wk <- merge(interval.results_all, pops_36wk_res_all %>%
                     rename("pops_statistic" = "statistic", "pops_pvalue" = "pvalue", "pops_beta" = "beta", "pops_se" = "se", 
                            "pops_threshold" = "threshold", "pops_minor" = "minor", "pops_major" = "major", "pops_snps" = "snps"), 
                   by=c("Pairs_notrsID", "gene", "chr", "SNPpos", "TSS"))


compare_eQTL <- function(comp, file_name){
  # check alleles/direction of effect
  comp$action <- "NA"
  comp[which(comp$pops_minor == comp$interval_minor), "action"] <- "same"
  
  # Alleles can't distinguish because of strand https://wanggroup.org/compbio_tutorial/allele_qc.html
  comp[which(comp$interval_minor=="A" & comp$interval_major=="T" |
               comp$interval_minor=="T" & comp$interval_major=="A" |
               comp$interval_minor=="C" & comp$interval_major=="G" |
               comp$interval_minor=="G" & comp$interval_major=="C"), "action"] <-"strand"
  
  # Same but different strand
  comp[which(comp$interval_minor=="A" & comp$pops_minor=="T" & comp$action != "strand" |
               comp$interval_minor=="T" & comp$pops_minor=="A" & comp$action != "strand" |
               comp$interval_minor=="C" & comp$pops_minor=="G" & comp$action != "strand" |
               comp$interval_minor=="G" & comp$pops_minor=="C" & comp$action != "strand"),
       "action"] <-"same"
  
  # Flip effects for those that are opposite
  comp[which(comp$action == "NA"), "action"] <- "flip"
  table(comp$action, useNA = "always")
  comp[which(comp$action == "flip"), "interval_beta"] <- comp[which(comp$action == "flip"), "interval_beta"] * (-1)
  
  # Remove SNPs which can't be distinguished because of strand
  comp$interval_beta[which(comp$action == "strand")] <- NA
  comp <- comp[complete.cases(comp), ]
  
  comp$interval.sig <- comp$interval_pvalue <= comp$interval_threshold
  comp$interval.sig[is.na(comp$interval.sig)] <- "FALSE"
  comp$pops.sig <- comp$pops_pvalue <= comp$pops_threshold
  comp_sig_tbl <- table(comp$interval.sig, comp$pops.sig)
  print(comp_sig_tbl)
  #        FALSE TRUE        
  # FALSE  2101 1086  
  # TRUE    483 9377  
  # out of 
  # TRUE TRUE 9377 genes are significant in both interval and pops
  # TRUE FALSE 483 genes are significant in interval but not significant in pops
  # FALSE TRUE 1086 are significant in POPS but not significant in interval
  # FALSE FALSE 2101 are not significant in either dataset 
  # Out of all genes that were significant in POPS (in this subsetted context: 9377 + 1086) = 10463,  1086 or 1086/10463=0.1037943=10.38% were pregnancy-specific
  print("writing out results")
  write.table(comp, paste0(dir, "input_data/interval_comp/", file_name, "_interval_pops_comp_lead_for_mashr_byPos_allLeads.txt"), sep="\t", quote=F)
  return(comp)
  
}

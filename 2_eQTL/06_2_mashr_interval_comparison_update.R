# 06_2_mashr_interval_comparison_update.R

################################################################################

# 6.2. Run mashR

################################################################################

# Aim: Use mashR to compare interval FRA and POPs2 TP eQTL effects

########################### Output paths ##########################

########################### Input paths ###########################
dir <-  "genotyping/analysis/eQTL/"
interval_dir <- "interval_rna/POPS2_comparison/"

########################### Parameters ############################
set.seed(14022022)

########################### Load packages ###########################
library(mashr)
library(tidyverse)
library(data.table)

########################### Load data and Analysis ###########################

comp_12wk_final <- read.table(paste0(dir, "input_data/interval_comp/12_weeks_interval_pops_comp_lead_for_mashr_byPos_allLeads.txt"), sep="\t")
comp_20wk_final <- read.table(paste0(dir, "input_data/interval_comp/20_weeks_interval_pops_comp_lead_for_mashr_byPos_allLeads.txt"), sep="\t")
comp_28wk_final <- read.table(paste0(dir, "input_data/interval_comp/28_weeks_interval_pops_comp_lead_for_mashr_byPos_allLeads.txt"), sep="\t")
comp_36wk_final <- read.table(paste0(dir, "input_data/interval_comp/36_weeks_interval_pops_comp_lead_for_mashr_byPos_allLeads.txt"), sep="\t")
setDT(comp_12wk_final)
setDT(comp_20wk_final)
setDT(comp_28wk_final)
setDT(comp_36wk_final)

# make a list of leads that have a significant effect
sig_leads <- unique(c(comp_12wk_final %>% filter(pops.sig == TRUE) %>% pull(Pairs_notrsID), 
                      comp_20wk_final %>% filter(pops.sig == TRUE) %>% pull(Pairs_notrsID), 
                      comp_28wk_final %>% filter(pops.sig == TRUE) %>% pull(Pairs_notrsID), 
                      comp_36wk_final %>% filter(pops.sig == TRUE) %>% pull(Pairs_notrsID)))

# Make strong tests df - Keep only the relevant columns and rename
tp12 <- comp_12wk_final[, .(Pairs_notrsID, gene, chr, SNPpos,
                            beta_12 = pops_beta, se_12 = pops_se,
                            interval_beta, interval_se)]
tp20 <- comp_20wk_final[, .(Pairs_notrsID, beta_20 = pops_beta, se_20 = pops_se)]
tp28 <- comp_28wk_final[, .(Pairs_notrsID, beta_28 = pops_beta, se_28 = pops_se)]
tp36 <- comp_36wk_final[, .(Pairs_notrsID, beta_36 = pops_beta, se_36 = pops_se)]

# Merge across timepoints
lead_snps <- Reduce(function(x, y) merge(x, y, by = "Pairs_notrsID", all = TRUE),
                    list(tp12, tp20, tp28, tp36))

# Reorder columns for clarity
setcolorder(lead_snps, c("Pairs_notrsID", "gene", "chr", "SNPpos",
                         "beta_12", "se_12", "beta_20", "se_20",
                         "beta_28", "se_28", "beta_36", "se_36",
                         "interval_beta", "interval_se"))

# filter lead SNPs to those that were significant at at least one time-point
lead_sig_snps <- lead_snps[Pairs_notrsID %in% sig_leads]

# Save sig tests for mashR training
write.table(lead_sig_snps, paste0(dir, "input_data/interval_comp/pops_leads_sig_interval_merged_for_mashr.txt"), sep = "\t")
lead_sig_snps <- read.table(paste0(dir, "input_data/interval_comp/pops_leads_sig_interval_merged_for_mashr.txt"), sep = "\t")
# Save all leads for final model fit
write.table(lead_snps, paste0(dir, "input_data/interval_comp/pops_leads_interval_merged_for_mashr_modelfit.txt"), sep = "\t")

pops.strong_12wk <- data.frame("beta"=lead_sig_snps$beta_12, "SE"=lead_sig_snps$se_12)
pops.strong_20wk <- data.frame("beta"=lead_sig_snps$beta_20, "SE"=lead_sig_snps$se_20)
pops.strong_28wk <- data.frame("beta"=lead_sig_snps$beta_28, "SE"=lead_sig_snps$se_28)
pops.strong_36wk <- data.frame("beta"=lead_sig_snps$beta_36, "SE"=lead_sig_snps$se_36)
interval.strong <- data.frame("beta"=lead_sig_snps$interval_beta, "SE"=lead_sig_snps$interval_se)

# random results - compare all results (initial pass) and sub-sample
# read in full results for each time-point and interval 
# read in full results and subset to pairs that are lead in at least one time-point
pops.bim <- data.frame(fread(paste0(dir, "input_data/single_tp_eQTL/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps.bim"),
                             header=FALSE))
setDT(pops.bim)
# Rename and select relevant columns (faster than dplyr)
pops.bim_small <- pops.bim[, .(snps = V2, minor = V5, major = V6)]
# 12 weeks
pops_12wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/12_weeks/12_weeks_ciseqtl_all.rds")) 
setDT(pops_12wk_full)
pops_12wk_full <- pops.bim_small[pops_12wk_full, on = "snps"]
pops_12wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
# 20 weeks
pops_20wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/20_weeks/20_weeks_ciseqtl_all.rds")) 
setDT(pops_20wk_full)
pops_20wk_full <- pops.bim_small[pops_20wk_full, on = "snps"]
pops_20wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
# 28 weeks
pops_28wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/28_weeks/28_weeks_ciseqtl_all.rds")) 
setDT(pops_28wk_full)
pops_28wk_full <- pops.bim_small[pops_28wk_full, on = "snps"]
pops_28wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]
# 36 weeks
pops_36wk_full <- readRDS(paste0(dir, "output_data/single_tp_analysis/36_weeks/36_weeks_ciseqtl_all.rds")) 
setDT(pops_36wk_full)
pops_36wk_full <- pops.bim_small[pops_36wk_full, on = "snps"]
pops_36wk_full[, Pairs_notrsID := paste0(chr, "_", SNPpos, "_", gene)]

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

# merge interval and 12 week results to harmonize alleles across datasets
# This works because the same snps were tested from the same dataset across time-points
comp_12wk <- merge(interval.results %>% rename("interval_statistic" = "statistic", "interval_pvalue" = "pvalue", "interval_beta" = "beta", "interval_se" = "se", 
                                               "interval_threshold" = "threshold", "interval_minor" = "minor", "interval_major" = "major", "interval_snps" = "snps"), 
                   pops_12wk_full %>% 
                     rename("pops_statistic" = "statistic", "pops_pvalue" = "pvalue", "pops_beta" = "beta", "pops_se" = "se", 
                            "pops_threshold" = "threshold", "pops_minor" = "minor", "pops_major" = "major", "pops_snps" = "snps"), 
                   by=c("Pairs_notrsID", "gene", "chr", "SNPpos", "TSS"))

# check alleles/direction of effect
comp_12wk$action <- "NA"
comp_12wk[which(comp_12wk$pops_minor == comp_12wk$interval_minor), "action"] <- "same"

# Alleles can't distinguish because of strand https://wanggroup.org/compbio_tutorial/allele_qc.html
comp_12wk[which(comp_12wk$interval_minor=="A" & comp_12wk$interval_major=="T" |
                  comp_12wk$interval_minor=="T" & comp_12wk$interval_major=="A" |
                  comp_12wk$interval_minor=="C" & comp_12wk$interval_major=="G" |
                  comp_12wk$interval_minor=="G" & comp_12wk$interval_major=="C"), "action"] <-"strand"

# Same but different strand
comp_12wk[which(comp_12wk$interval_minor=="A" & comp_12wk$pops_minor=="T" & comp_12wk$action != "strand" |
                  comp_12wk$interval_minor=="T" & comp_12wk$pops_minor=="A" & comp_12wk$action != "strand" |
                  comp_12wk$interval_minor=="C" & comp_12wk$pops_minor=="G" & comp_12wk$action != "strand" |
                  comp_12wk$interval_minor=="G" & comp_12wk$pops_minor=="C" & comp_12wk$action != "strand"),
          "action"] <-"same"

# Flip effects for those that are opposite
comp_12wk[which(comp_12wk$action == "NA"), "action"] <- "flip"
table(comp_12wk$action, useNA = "always")
comp_12wk[which(comp_12wk$action == "flip"), "interval_beta"] <- comp_12wk[which(comp_12wk$action == "flip"), "interval_beta"] * (-1)

# Remove SNPs which can't be distinguished because of strand
comp_12wk$interval_beta[which(comp_12wk$action == "strand")] <- NA
comp_12wk <- comp_12wk[complete.cases(comp_12wk), ]

saveRDS(comp_12wk, paste0(dir, "input_data/interval_comp/12wk_interval_pops_comp_for_mashr_byPos_allPairs.rds"))

# now, join this to the remaining eQTL from 20, 28 and 36 weeks
# Make sure all datasets are data.tables
setDT(comp_12wk)       
setDT(pops_20wk_full)       
setDT(pops_28wk_full)       
setDT(pops_36wk_full)       

# 1. Rename columns in tp datasets so that after join, they reflect the time-point
setnames(comp_12wk, old = c("pops_statistic", "pops_pvalue", "pops_beta", "pops_se"),
         new = c("pops_statistic_12wk", "pops_pvalue_12wk", "pops_beta_12wk", "pops_se_12wk"))

setnames(pops_20wk_full, old = c("statistic", "pvalue", "beta", "se", "snps"),
         new = c("pops_statistic_20wk", "pops_pvalue_20wk", "pops_beta_20wk", "pops_se_20wk", "pops_snps"))

setnames(pops_28wk_full, old = c("statistic", "pvalue", "beta", "se", "snps"),
         new = c("pops_statistic_28wk", "pops_pvalue_28wk", "pops_beta_28wk", "pops_se_28wk", "pops_snps"))

setnames(pops_36wk_full, old = c("statistic", "pvalue", "beta", "se", "snps"),
         new = c("pops_statistic_36wk", "pops_pvalue_36wk", "pops_beta_36wk", "pops_se_36wk", "pops_snps"))

# Remove redundant columns (minor, major, threshold) from single TP datasets
pops_20wk_full[, c("minor", "major", "threshold") := NULL]
pops_28wk_full[, c("minor", "major", "threshold") := NULL]
pops_36wk_full[, c("minor", "major", "threshold") := NULL]

# Merge sequentially by key columns
key_cols <- c("Pairs_notrsID", "pops_snps", "gene", "chr", "SNPpos", "TSS")

merged_df <- merge(comp_12wk, pops_20wk_full, by = key_cols, all.x = TRUE)
merged_df <- merge(merged_df, pops_28wk_full, by = key_cols, all.x = TRUE)
merged_df <- merge(merged_df, pops_36wk_full, by = key_cols, all.x = TRUE)

saveRDS(merged_df, paste0(dir, "input_data/interval_comp/alltp_interval_pops_comp_for_mashr_byPos_allPairs.rds"))
merged_df <- readRDS(paste0(dir, "input_data/interval_comp/alltp_interval_pops_comp_for_mashr_byPos_allPairs.rds"))

# sub-sample random results from 
merged_df_sub <- merged_df[sample(1:nrow(merged_df), size=nrow(merged_df)*0.01), ]
# write this out to be read back in 
saveRDS(merged_df_sub, paste0(dir, "input_data/interval_comp/sub_sampled_random_results_for_mashr.rds"))

pops.random_12wk <- data.frame("beta"=merged_df_sub$pops_beta_12wk, "SE"=merged_df_sub$pops_se_12wk)
pops.random_20wk <- data.frame("beta"=merged_df_sub$pops_beta_20wk, "SE"=merged_df_sub$pops_se_20wk)
pops.random_28wk <- data.frame("beta"=merged_df_sub$pops_beta_28wk, "SE"=merged_df_sub$pops_se_28wk)
pops.random_36wk <- data.frame("beta"=merged_df_sub$pops_beta_36wk, "SE"=merged_df_sub$pops_se_36wk)
interval.random <- data.frame("beta"=merged_df_sub$interval_beta, "SE"=merged_df_sub$interval_se)

# set up data for mash
effects.strong <- as.matrix(cbind(pops.strong_12wk$beta, pops.strong_20wk$beta, pops.strong_28wk$beta, pops.strong_36wk$beta, interval.strong$beta))
ses.strong <- as.matrix(cbind(pops.strong_12wk$SE, pops.strong_20wk$SE, pops.strong_28wk$SE, pops.strong_36wk$SE, interval.strong$SE))
colnames(effects.strong) <- c("POPS2_12wk", "POPS2_20wk", "POPS2_28wk", "POPS2_36wk", "Interval")
colnames(ses.strong) <- c("POPS2_12wk", "POPS2_20wk", "POPS2_28wk", "POPS2_36wk", "Interval")
data.strong <- mash_set_data(effects.strong, ses.strong)

effects.random <- as.matrix(cbind(pops.random_12wk$beta, pops.random_20wk$beta, pops.random_28wk$beta, pops.random_36wk$beta, interval.random$beta))
ses.random <- as.matrix(cbind(pops.random_12wk$SE, pops.random_20wk$SE, pops.random_28wk$SE, pops.random_36wk$SE, interval.random$SE))
colnames(effects.random) <- c("POPS2_12wk", "POPS2_20wk", "POPS2_28wk", "POPS2_36wk", "Interval")
colnames(ses.random) <- c("POPS2_12wk", "POPS2_20wk", "POPS2_28wk", "POPS2_36wk", "Interval")
data.random <- mash_set_data(effects.random, ses.random)

# account for correlation structure among measurements (i.e. across time-points)
# find correlation structure from random subset
# as described in this tutorial https://stephenslab.github.io/mashr/articles/intro_correlations.html
Vhat = estimate_null_correlation_simple(data.random)

# Set up our main data objects with this correlation structure in place:
data.random = mash_set_data(effects.random, ses.random, V=Vhat)
data.strong = mash_set_data(effects.strong, ses.strong, V=Vhat)

# Use the strong tests to set up data-driven covariances
# as described in this tutorial https://stephenslab.github.io/mashr/articles/intro_mash_dd.html
# here, I choose to use 5 PCs because I am mashing 5 sets of summary stats so it reflects the possible axes of variation in the data
# (Compared to Katie's choice of 2 because she was mashing 2 datasets)
U.pca = cov_pca(data.strong, 5)
print(names(U.pca))
U.ed = cov_ed(data.strong, U.pca)

# Now we fit mash to the random tests using both data-driven and canonical covariances. 
# (Remember the Crucial Rule! We have to fit using a random set of tests, and not a 
# dataset that is enriched for strong tests.) The outputlevel=1 option means that it 
# will not compute posterior summaries for these tests (which saves time).
# It is recommend to run mash with both data-driven and canonical covariances (U.c and U.ed)
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
# - Computing 414277 x 497 likelihood matrix.
# - Likelihood calculations took 459.82 seconds.
# - Fitting model with 497 mixture components.
# - Model fitting took 2503.31 seconds. (41 mins)

# Now we can compute posterior summaries etc for any subset of tests using the 
# above mash fit. Here we do this for the strong tests. We do this using the same 
# mash function as above, but we specify to use the fit from the previous run of 
# mash by specifying g=get_fitted_g(m), fixg=TRUE. (In mash the parameter g is 
# used to denote the mixture model which we learned above.)

# TO DO decidie if I need to re-make data.strong but with all of the strong tests, not just the significant ones 
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
# - Computing 17704 x 497 likelihood matrix.
# - Likelihood calculations took 18.62 seconds.
# - Computing posterior matrices.
# - Computation allocated took 2.05 seconds.

print(get_pairwise_sharing(m2)) # the same sign and within a factor 0.5 of each other
#             POPS2_12wk POPS2_20wk POPS2_28wk POPS2_36wk  Interval
# POPS2_12wk  1.0000000  0.9924515  0.9906818  0.9866148 0.8557067
# POPS2_20wk  0.9924515  1.0000000  0.9947536  0.9921023 0.8408876
# POPS2_28wk  0.9906818  0.9947536  1.0000000  0.9937496 0.8322904
# POPS2_36wk  0.9866148  0.9921023  0.9937496  1.0000000 0.8264848
# Interval    0.8557067  0.8408876  0.8322904  0.8264848 1.0000000
# list of different ones:
post.mean <- data.frame(get_pm(m2))
rownames(post.mean) <- lead_sig_snps$Pairs_notrsID

# column names
cols <- colnames(post.mean)
# Create empty lists to store results
ratio_list <- list()
diff_list <- list()

# Loop over all pairwise combinations of columns
for (pair in combn(cols, 2, simplify = FALSE)) {
  col1 <- pair[1]
  col2 <- pair[2]
  
  # Column names for ratio and diff
  ratio_name <- paste0("ratio_", col1, "_vs_", col2)
  diff_name <- paste0("diff_", col1, "_vs_", col2)
  
  # Calculate ratio
  post.mean[[ratio_name]] <- post.mean[[col1]] / post.mean[[col2]]
  
  # Flag extreme differences (<0.5 or >2)
  post.mean[[diff_name]] <- post.mean[[ratio_name]] < 0.5 | post.mean[[ratio_name]] > 2
}

# Quick summary: number of extreme differences per pair
diff_cols <- grep("^diff_", colnames(post.mean), value = TRUE)
sapply(post.mean[, diff_cols], table)

#       diff_POPS2_12wk_vs_POPS2_20wk diff_POPS2_12wk_vs_POPS2_28wk
# FALSE                         17283                         17257
# TRUE                            421                           447
#       diff_POPS2_12wk_vs_POPS2_36wk diff_POPS2_12wk_vs_Interval
# FALSE                         17171                       14787
# TRUE                            533                        2917
#       diff_POPS2_20wk_vs_POPS2_28wk diff_POPS2_20wk_vs_POPS2_36wk
# FALSE                         17356                         17275
# TRUE                            348                           429
#       diff_POPS2_20wk_vs_Interval diff_POPS2_28wk_vs_POPS2_36wk
# FALSE                       14545                         17325
# TRUE                         3159                           379
#       diff_POPS2_28wk_vs_Interval diff_POPS2_36wk_vs_Interval
# FALSE                       14409                       14295
# TRUE                         3295                        3409

post.sig <- get_lfsr(m2)
colnames(post.sig) <- paste0(colnames(post.sig), ".lfsr")
post.mean <- data.frame(post.mean, post.sig)

theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

effect_cols <- c("POPS2_12wk","POPS2_20wk","POPS2_28wk","POPS2_36wk","Interval")

# All pairwise combinations
pairs <- t(combn(effect_cols, 2))
pairs_df <- data.frame(xvar = pairs[,1], yvar = pairs[,2], stringsAsFactors = FALSE)

# Build long dataframe for plotting
plot_data <- do.call(rbind, lapply(1:nrow(pairs_df), function(i){
  x <- pairs_df$xvar[i]
  y <- pairs_df$yvar[i]
  data.frame(
    xval = post.mean[[x]],
    yval = post.mean[[y]],
    diff = abs(post.mean[[x]] / post.mean[[y]]) < 0.5 | abs(post.mean[[x]] / post.mean[[y]]) > 2,
    xvar = x,
    yvar = y
  )
}))

# Keep only upper-triangular pairs: xvar index < yvar index
plot_data <- plot_data %>%
  filter(match(xvar, effect_cols) < match(yvar, effect_cols))

# Drop unused factor levels
plot_data$xvar <- factor(plot_data$xvar, levels = effect_cols)
plot_data$yvar <- factor(plot_data$yvar, levels = effect_cols)

# Plot
ggplot(plot_data, aes(x = xval, y = yval, color = diff)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(size = 0.01) +
  scale_color_manual(values = c("lightgrey", "darkblue")) +
  facet_grid(yvar ~ xvar, scales = "free", drop = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) + ylab(NULL) 

write.table(post.mean, paste0(dir, "output_data/interval_comp/pops_alltps_intervalFRA_mashr_results.txt"), sep="\t", quote=F, row.names = TRUE)

lead_sig_snps$diff_POPS2_12wk_vs_POPS2_20wk_mashr <- post.mean$diff_POPS2_12wk_vs_POPS2_20wk
lead_sig_snps$diff_POPS2_12wk_vs_POPS2_28wk_mashr <- post.mean$diff_POPS2_12wk_vs_POPS2_28wk
lead_sig_snps$diff_POPS2_12wk_vs_POPS2_36wk_mashr <- post.mean$diff_POPS2_12wk_vs_POPS2_36wk
lead_sig_snps$diff_POPS2_12wk_vs_Interval_mashr <- post.mean$diff_POPS2_12wk_vs_Interval
lead_sig_snps$diff_POPS2_20wk_vs_POPS2_28wk_mashr <- post.mean$diff_POPS2_20wk_vs_POPS2_28wk
lead_sig_snps$diff_POPS2_20wk_vs_POPS2_36wk_mashr <- post.mean$diff_POPS2_20wk_vs_POPS2_36wk
lead_sig_snps$diff_POPS2_20wk_vs_Interval_mashr <- post.mean$diff_POPS2_20wk_vs_Interval
lead_sig_snps$diff_POPS2_28wk_vs_POPS2_36wk_mashr <- post.mean$diff_POPS2_28wk_vs_POPS2_36wk
lead_sig_snps$diff_POPS2_28wk_vs_Interval_mashr <- post.mean$diff_POPS2_28wk_vs_Interval
lead_sig_snps$diff_POPS2_36wk_vs_Interval_mashr <- post.mean$diff_POPS2_36wk_vs_Interval

interval_diff_cols <- grep("_vs_Interval_mashr$", colnames(lead_sig_snps), value = TRUE)

lead_sig_snps_interval_diff <- lead_sig_snps[rowSums(lead_sig_snps[, interval_diff_cols]) > 0, ]
write.table(lead_sig_snps_interval_diff, paste0(dir, "output_data/interval_comp/pregnancy_specific_pops2_lead_mashr.txt"), 
            sep="\t", quote=F)

print(get_pairwise_sharing(m2, factor=0)) # same sign
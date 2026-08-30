# 06_3_defining_pregnancy_specificity.R

################################################################################

# 6.3. Define context-specific eQTL signals 

################################################################################

# Aim: Using mashR output, categorize eQTL context-specificity 

########################### Paths ##########################
dir <- "genotyping/analysis/eQTL/"

########################### Load packages ###########################
library(tidyverse)


theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Analysis ###########################

# read in mashR results (across four time-points and interval FRA)
mashr.results <- read.table(paste0(dir, "output_data/interval_comp/pops_alltps_intervalFRA_mashr_results.txt"), sep="\t")
# subset mashR results to only comparisons between Interval and POPs
mashr.results_vInterval <- mashr.results %>% 
  select(POPS2_12wk, POPS2_20wk, POPS2_28wk, POPS2_36wk, Interval, 
         ratio_POPS2_12wk_vs_Interval, ratio_POPS2_20wk_vs_Interval, ratio_POPS2_28wk_vs_Interval, ratio_POPS2_36wk_vs_Interval, 
         diff_POPS2_12wk_vs_Interval, diff_POPS2_20wk_vs_Interval, diff_POPS2_28wk_vs_Interval, diff_POPS2_36wk_vs_Interval, 
         POPS2_12wk.lfsr, POPS2_20wk.lfsr, POPS2_28wk.lfsr, POPS2_36wk.lfsr, Interval.lfsr)

# define sharing at each time-point vs. interval comparison into the following categories: 
# 1. Not significant in POPs after shrinkage 
# 2. Significant in POPs but not Interval 
# 3. Signal enhanced in pops
# 4. Signal dampened in pops
# 5. Opposite direction of effect in pops and interval 
# 6. Shared effect in pops and interval 

################################################################################
# 12 weeks
################################################################################
# add sharing column
mashr.results_vInterval$shared_POPS2_12wk_vs_Interval <- NA
# Where pops signal isn't significant, mark "not significant"
mashr.results_vInterval$shared_POPS2_12wk_vs_Interval[which(mashr.results_vInterval$POPS2_12wk.lfsr > 0.05)] <- "not significant"
# If pops signal is significant, use diff column to define whether the signal is shared or not shared
mashr.results_vInterval$shared_POPS2_12wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_12wk_vs_Interval)] <- ifelse(
  mashr.results_vInterval$diff_POPS2_12wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_12wk_vs_Interval)] == "FALSE", "shared", "not shared")
# If both pops and interval signals are significant but the ratio is less than 0 (indicating different signs across datasets), mark opposite
mashr.results_vInterval$shared_POPS2_12wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_12wk_vs_Interval < 0 &
                                                              mashr.results_vInterval$Interval.lfsr < 0.05 &
                                                              mashr.results_vInterval$POPS2_12wk.lfsr < 0.05)] <- "Opposite direction of effect"
# If the ratio is less than 0 (indicating opposite direction of effect) and only the POPs signal is significant, indicate only significant in POPS
mashr.results_vInterval$shared_POPS2_12wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_12wk_vs_Interval < 0 & 
                                                              mashr.results_vInterval$Interval.lfsr >= 0.05 &
                                                              mashr.results_vInterval$POPS2_12wk.lfsr < 0.05)] <- "Only significant in POPS2"
# Divide the not-shared results by which dataset has a bigger effect
mashr.results_vInterval$shared_POPS2_12wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_12wk_vs_Interval == "not shared")] <- ifelse(
  mashr.results_vInterval$ratio_POPS2_12wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_12wk_vs_Interval == "not shared")] >2, 
  "Bigger effect in POPS2", 
  "Bigger effect in Interval")

################################################################################
# 20 weeks
################################################################################
# add sharing column
mashr.results_vInterval$shared_POPS2_20wk_vs_Interval <- NA
# Where pops signal isn't significant, mark "not significant"
mashr.results_vInterval$shared_POPS2_20wk_vs_Interval[which(mashr.results_vInterval$POPS2_20wk.lfsr > 0.05)] <- "not significant"
# If pops signal is significant, use diff column to define whether the signal is shared or not shared
mashr.results_vInterval$shared_POPS2_20wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_20wk_vs_Interval)] <- ifelse(
  mashr.results_vInterval$diff_POPS2_20wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_20wk_vs_Interval)] == "FALSE", "shared", "not shared")
# If both pops and interval signals are significant but the ratio is less than 0 (indicating different signs across datasets), mark opposite
mashr.results_vInterval$shared_POPS2_20wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_20wk_vs_Interval < 0 &
                                                              mashr.results_vInterval$Interval.lfsr < 0.05 &
                                                              mashr.results_vInterval$POPS2_20wk.lfsr < 0.05)] <- "Opposite direction of effect"
# If the ratio is less than 0 (indicating opposite direction of effect) and only the POPs signal is significant, indicate only significant in POPS
mashr.results_vInterval$shared_POPS2_20wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_20wk_vs_Interval < 0 & 
                                                              mashr.results_vInterval$Interval.lfsr >= 0.05 &
                                                              mashr.results_vInterval$POPS2_20wk.lfsr < 0.05)] <- "Only significant in POPS2"
# Divide the not-shared results by which dataset has a bigger effect
mashr.results_vInterval$shared_POPS2_20wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_20wk_vs_Interval == "not shared")] <- ifelse(
  mashr.results_vInterval$ratio_POPS2_20wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_20wk_vs_Interval == "not shared")] > 2, 
  "Bigger effect in POPS2", 
  "Bigger effect in Interval")

################################################################################
# 28 weeks
################################################################################
# add sharing column
mashr.results_vInterval$shared_POPS2_28wk_vs_Interval <- NA
# Where pops signal isn't significant, mark "not significant"
mashr.results_vInterval$shared_POPS2_28wk_vs_Interval[which(mashr.results_vInterval$POPS2_28wk.lfsr > 0.05)] <- "not significant"
# If pops signal is significant, use diff column to define whether the signal is shared or not shared
mashr.results_vInterval$shared_POPS2_28wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_28wk_vs_Interval)] <- ifelse(
  mashr.results_vInterval$diff_POPS2_28wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_28wk_vs_Interval)] == "FALSE", "shared", "not shared")
# If both pops and interval signals are significant but the ratio is less than 0 (indicating different signs across datasets), mark opposite
mashr.results_vInterval$shared_POPS2_28wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_28wk_vs_Interval < 0 &
                                                              mashr.results_vInterval$Interval.lfsr < 0.05 &
                                                              mashr.results_vInterval$POPS2_28wk.lfsr < 0.05)] <- "Opposite direction of effect"
# If the ratio is less than 0 (indicating opposite direction of effect) and only the POPs signal is significant, indicate only significant in POPS
mashr.results_vInterval$shared_POPS2_28wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_28wk_vs_Interval < 0 & 
                                                              mashr.results_vInterval$Interval.lfsr >= 0.05 &
                                                              mashr.results_vInterval$POPS2_28wk.lfsr < 0.05)] <- "Only significant in POPS2"
# Divide the not-shared results by which dataset has a bigger effect
mashr.results_vInterval$shared_POPS2_28wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_28wk_vs_Interval == "not shared")] <- ifelse(
  mashr.results_vInterval$ratio_POPS2_28wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_28wk_vs_Interval == "not shared")] > 2, 
  "Bigger effect in POPS2", 
  "Bigger effect in Interval")
################################################################################
# 36 weeks
################################################################################
# add sharing column
mashr.results_vInterval$shared_POPS2_36wk_vs_Interval <- NA
# Where pops signal isn't significant, mark "not significant"
mashr.results_vInterval$shared_POPS2_36wk_vs_Interval[which(mashr.results_vInterval$POPS2_36wk.lfsr > 0.05)] <- "not significant"
# If pops signal is significant, use diff column to define whether the signal is shared or not shared
mashr.results_vInterval$shared_POPS2_36wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_36wk_vs_Interval)] <- ifelse(
  mashr.results_vInterval$diff_POPS2_36wk_vs_Interval[is.na(mashr.results_vInterval$shared_POPS2_36wk_vs_Interval)] == "FALSE", "shared", "not shared")
# If both pops and interval signals are significant but the ratio is less than 0 (indicating different signs across datasets), mark opposite
mashr.results_vInterval$shared_POPS2_36wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_36wk_vs_Interval < 0 &
                                                              mashr.results_vInterval$Interval.lfsr < 0.05 &
                                                              mashr.results_vInterval$POPS2_36wk.lfsr < 0.05)] <- "Opposite direction of effect"
# If the ratio is less than 0 (indicating opposite direction of effect) and only the POPs signal is significant, indicate only significant in POPS
mashr.results_vInterval$shared_POPS2_36wk_vs_Interval[which(mashr.results_vInterval$ratio_POPS2_36wk_vs_Interval < 0 & 
                                                              mashr.results_vInterval$Interval.lfsr >= 0.05 &
                                                              mashr.results_vInterval$POPS2_36wk.lfsr < 0.05)] <- "Only significant in POPS2"
# Divide the not-shared results by which dataset has a bigger effect
mashr.results_vInterval$shared_POPS2_36wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_36wk_vs_Interval == "not shared")] <- ifelse(
  mashr.results_vInterval$ratio_POPS2_36wk_vs_Interval[which(mashr.results_vInterval$shared_POPS2_36wk_vs_Interval == "not shared")] > 2, 
  "Bigger effect in POPS2", 
  "Bigger effect in Interval")

################################################################################
# Processing across all time-points
################################################################################
# write out the mashR results with sharing defined in Interval 
write.table(mashr.results_vInterval, paste0(dir, "output_data/interval_comp/pops_alltps_intervalFRA_mashr_results_sharing.txt"), sep="\t", quote = FALSE)

# read in lead SNP gene pairs per time-point for filtering
# POPS lead SNP-gene pairs WHERE tested in interval and no strand issues in interval
pops_12wk_sig_leads <- read.delim(paste0(dir, "output_data/single_tp_analysis/12_weeks/eigenMT/12_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  filter(Sig == TRUE) %>% mutate(Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene), lead_POPS2_12wk = "POPS2_12wk") %>% select(Pairs_notrsID, lead_POPS2_12wk)
pops_20wk_sig_leads <- read.delim(paste0(dir, "output_data/single_tp_analysis/20_weeks/eigenMT/20_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  filter(Sig == TRUE) %>% mutate(Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene), lead_POPS2_20wk = "POPS2_20wk") %>% select(Pairs_notrsID, lead_POPS2_20wk)
pops_28wk_sig_leads <- read.delim(paste0(dir, "output_data/single_tp_analysis/28_weeks/eigenMT/28_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  filter(Sig == TRUE) %>% mutate(Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene), lead_POPS2_28wk = "POPS2_28wk") %>% select(Pairs_notrsID, lead_POPS2_28wk)
pops_36wk_sig_leads <- read.delim(paste0(dir, "output_data/single_tp_analysis/36_weeks/eigenMT/36_weeks_ciseqtl_eigenMT_corrected.txt")) %>% 
  filter(Sig == TRUE) %>% mutate(Pairs_notrsID = paste0(chr, "_", SNPpos, "_", gene), lead_POPS2_36wk = "POPS2_36wk") %>% select(Pairs_notrsID, lead_POPS2_36wk)

# now, make the dataset longer and subset to the lead pair per time-point
mashr.results_vInterval_long <- mashr.results_vInterval %>%
  rownames_to_column(var = "Pairs_notrsID") %>%
  select(
    Pairs_notrsID, POPS2_12wk, POPS2_20wk, POPS2_28wk, POPS2_36wk, Interval,
    shared_POPS2_12wk_vs_Interval, shared_POPS2_20wk_vs_Interval,
    shared_POPS2_28wk_vs_Interval, shared_POPS2_36wk_vs_Interval
  ) %>%
  pivot_longer(
    cols = starts_with("POPS2_"),   # the POPS2 effect size columns
    names_to = "POPS2_timepoint",
    values_to = "POPS2_effect"
  ) %>%
  mutate(
    # Map the correct shared column to each POPS2 timepoint
    shared = case_when(
      POPS2_timepoint == "POPS2_12wk" ~ shared_POPS2_12wk_vs_Interval,
      POPS2_timepoint == "POPS2_20wk" ~ shared_POPS2_20wk_vs_Interval,
      POPS2_timepoint == "POPS2_28wk" ~ shared_POPS2_28wk_vs_Interval,
      POPS2_timepoint == "POPS2_36wk" ~ shared_POPS2_36wk_vs_Interval
    )
  ) %>% 
  select(-c(shared_POPS2_12wk_vs_Interval, shared_POPS2_20wk_vs_Interval, shared_POPS2_28wk_vs_Interval, shared_POPS2_36wk_vs_Interval)) %>%
  mutate(shared = factor(shared, levels = c("Bigger effect in POPS2", "Only significant in POPS2", "Bigger effect in Interval", "Opposite direction of effect", "shared", "not significant"))) %>% 
  arrange(shared) %>%
  left_join(pops_12wk_sig_leads, by = "Pairs_notrsID") %>%
  left_join(pops_20wk_sig_leads, by = "Pairs_notrsID") %>%
  left_join(pops_28wk_sig_leads, by = "Pairs_notrsID") %>%
  left_join(pops_36wk_sig_leads, by = "Pairs_notrsID") %>%
  filter(lead_POPS2_12wk == POPS2_timepoint | lead_POPS2_20wk == POPS2_timepoint | lead_POPS2_28wk == POPS2_timepoint | lead_POPS2_36wk == POPS2_timepoint)

# make sure the numbers work 
mashr.results_vInterval_long %>% separate(col = Pairs_notrsID,
                                          into = c("chr", "pos", "gene"),
                                          sep = "_",
                                          remove = FALSE) %>% group_by(gene) %>% summarize(n = n()) %>% filter(n > 4)
# how many signals remain per time-point
mashr.results_vInterval_long %>% group_by(POPS2_timepoint) %>% summarize(n = n())

# write this out for enrichment testing 
mashr.results_vInterval_long %>%
  select(-c(lead_POPS2_12wk, lead_POPS2_20wk, lead_POPS2_28wk, lead_POPS2_36wk)) %>%
  write.table(paste0(dir, "output_data/interval_comp/pops_alltps_intervalFRA_mashr_results_sharing_filtered_leads.txt"), sep="\t", quote = FALSE)

########################### Plot ###########################

# read in mashR results
mashr.results_Interval <- read.table(paste0(dir, "output_data/interval_comp/pops_alltps_intervalFRA_mashr_results_sharing_filtered_leads.txt"), sep="\t") %>%
  mutate(shared = factor(shared, levels = c("Bigger effect in POPS2", "Only significant in POPS2", "Bigger effect in Interval", "Opposite direction of effect", "shared", "not significant")))

# Scatterplot of context categories (Main Figure 5A)
mashr.results_Interval %>%
  filter(shared != "not significant") %>%
  mutate(POPS2_timepoint = paste0(gsub("POPS2_", "", POPS2_timepoint), "GA")) %>%
  arrange(desc(shared == "shared")) %>%
  ggplot(aes(Interval, POPS2_effect)) +
  geom_point(aes(colour=shared, fill = shared, alpha = shared), pch = 21) +
  scale_colour_manual(values=c("#598eca", "#96b3d6", "#cb2f43", "goldenrod", "gray", "black"), labels = c("Pregnancy-magnified", "Pregnancy-magnified\n(not significant in INTERVAL)", "Pregnancy-dampened", "Opposite direction of effect", "Shared", "Not significant")) +
  scale_fill_manual(values=c("#598eca", "#96b3d6", "#cb2f43", "goldenrod", "gray", "transparent"), labels = c("Pregnancy-magnified", "Pregnancy-magnified\n(not significant in INTERVAL)", "Pregnancy-dampened", "Opposite direction of effect", "Shared", "Not significant")) +
  scale_alpha_manual(values = c(1, 1, 1, 1, 1, 0.6), labels = c("Pregnancy-magnified", "Pregnancy-magnified\n(not significant in Interval)", "Pregnancy-dampened", "Opposite direction of effect", "Shared", "Not significant")) +
  #ggtitle("Pregnancy vs non-pregnancy effect sizes") +
  ylab("POPs2 posterior effect size") +
  xlab("Interval posterior effect size") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_abline(intercept = 0, slope=1, lty=2, linewidth = 0.3) +
  geom_abline(intercept = 0, slope=0.5, lty=2, linewidth = 0.3) +
  geom_abline(intercept = 0, slope=1/0.5, lty=2, linewidth = 0.3) +
  # theme(
  #   legend.position = c(0.01, 0.95),      # x, y coordinates (relative to plot area)
  #   legend.justification = c("left", "top")
  # ) +
  facet_wrap(~ POPS2_timepoint, scales = "free_x", ncol = 4) +
  labs(color = "Context dependency",
       fill = "Context dependency") +
  guides(color = guide_legend(override.aes = list(size = 4)),
         alpha = "none")  +
  theme(strip.text = element_text(size = 10)) 
#legend.position = "top")  # Increase the font size here)

# plot number of signals in each category per time-point
df_perc <- mashr.results_Interval %>%
  filter(shared != "not significant") %>%
  count(shared, POPS2_timepoint) %>%                    # count signals per group
  group_by(POPS2_timepoint) %>%                         # group by timepoint
  mutate(percent = 100 * n / sum(n)) %>%                # convert to %
  ungroup()

# Plot n signals in each category across time-points (Main Figure 5B)
df_perc %>%
  mutate(POPS2_timepoint = gsub("wk", " weeks", gsub("POPS2_", "", POPS2_timepoint)), 
         shared = ifelse(shared == "Bigger effect in POPS2", "Pregnancy-magnified", 
                         ifelse(shared == "Only significant in POPS2", "Pregnancy-magnified\n(not significant in INTERVAL)", 
                                ifelse(shared == "Bigger effect in Interval", "Pregnancy-dampened", 
                                       ifelse(shared == "shared", "Shared", 
                                              ifelse(shared == "not significant", "Not significant", "Opposite direction of effect"))))),
         shared = factor(shared, levels = c("Shared", "Pregnancy-magnified", "Not significant", "Pregnancy-magnified\n(not significant in INTERVAL)", "Opposite direction of effect", "Pregnancy-dampened"))) %>% 
  ggplot(aes(x = shared, y = percent, fill = POPS2_timepoint)) +
  geom_col(position = "dodge") +                        # side-by-side bars
  geom_text(
    aes(label = paste0(sprintf("%.1f", percent), "%")),             # label format
    position = position_dodge(width = 0.9), 
    vjust = -0.3,                                        # above bars
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "12 weeks" = "#BCE4D8",
      "20 weeks" = "#83C4CB",
      "28 weeks" = "#439FB7",
      "36 weeks" = "#32769B"
    ), 
    labels = c(
      "12 weeks" = "12wkGA",
      "20 weeks" = "20wkGA",
      "28 weeks" = "28wkGA",
      "36 weeks" = "36wkGA"
    )
  ) +
  labs(
    x = "Context dependency category",
    y = "Percent of signals",
    #title = "Quantification of eQTL context specificity per time-point",
    fill = "POPS2 timepoint"
  ) + 
  theme(
    legend.position = c(0.98, 0.95),      # x, y coordinates (relative to plot area)
    legend.justification = c("right", "top"),
    axis.text.x = element_text(color="black")
  ) 

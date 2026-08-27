# 02_5_cibersort_continuous_trajectories.R 

################################################################################

# 2.5. Continuous trajectories of imputed cell proportions

################################################################################

# Aim: Plot continuous trajectories of cibersort outputs and validate them compared to Bar et al 2025

########################### Input paths ###########################
cibersort_inpath <- "rna-seq/analysis/cibersort/outputs/CIBERSORTx_Adjusted.txt"
sample_info_inpath <- "rna-seq/data/sample_covariates_clin_tech.csv"
bar_cellcount_inpath <- "rna-seq/analysis/outcome_pred/Ref_data/Bar_2025_cell_props/"

########################### Load packages ###########################
library(tidyverse)
library(variancePartition)
library(ggeffects)
library(interactions)

# set ggplot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Load data ###########################

# read in cibersort output imputed cell proportions
cell_prop_only <- read.delim(cibersort_inpath)
# read in sample info 
sample.info <- read.csv(sample_info_inpath) 
# bar cell count data for validation 
Bar_eos <- read.csv(paste0(bar_cellcount_inpath, "EOS_perc.csv"))
Bar_lym <- read.csv(paste0(bar_cellcount_inpath, "LYMperc.csv"))
Bar_mon <- read.csv(paste0(bar_cellcount_inpath, "MONperc.csv"))
Bar_neu <- read.csv(paste0(bar_cellcount_inpath, "NEUTperc.csv"))

########################### Analysis ###########################
cell_prop <- cell_prop_only %>% left_join(sample.info,  by=c("Mixture" ="RNA_sanger_sample_id"))


# Variance parition - determine which covariates belong in model
# transpose cell proportions for input into variance partition
var_part_cell_prop <- cell_prop %>% 
  column_to_rownames(var = "Mixture") %>% 
  select(contains("cells") | contains("Macrophage") | contains("phils") | contains("Monocytes"),
         # remove cells that showed no differential abundance over time
         -c(Dendritic.cells.resting, Mast.cells.activated, T.cells.follicular.helper, T.cells.gamma.delta)) %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(.)))

# make variance partition formula
variance_partition_formula <- ~(1|Sample_taken_at)+(1|id_run_position)+(1|ANON_ID)+(1|an_alcohol_status)+an_est_age+an_scan1_bmi+(1|an_level_edu_str)+(1|an_diab_tf)+pn_monthdel+(1|an_eth_cat_str)+(1|an_smokstat)

# make sample info df with clinical covariates 
sample_info_vp <- sample.info %>%
  mutate(an_diab_tf = as.logical(an_diab),
         id_run = as.factor(id_run)) %>% 
  arrange(factor(RNA_sanger_sample_id, levels = colnames(var_part_cell_prop))) %>%
  column_to_rownames("RNA_sanger_sample_id")

# run variance partition on cell proportions (takes ~one minute)
varPart <- fitExtractVarPartModel(var_part_cell_prop, variance_partition_formula, sample_info_vp)

# plot variance parititon (Supplementary Figure 22)
plotVarPart(sortCols(varPart), , col = c(rep("white", 11), "grey85")) + 
  scale_x_discrete(labels=c("ANON_ID" = "Individual", 
                            "id_run_position" = "Sequencing\nbatch", 
                            #"id_run" = "Sequencing\nflowcell",
                            "Sample_taken_at" = "Sample\ntime-point", 
                            "an_eth_cat_str" = "Self-reported\nethnicity", 
                            "an_scan1_bmi" = "BMI", 
                            "pn_monthdel" = "Month of\ndelivery", 
                            "an_est_age" = "Age", 
                            "an_level_edu_str" = "Education\nlevel", 
                            "an_alcohol_status" = "Alcohol\nstatus", 
                            "an_smokstat" = "Smoking\nstatus", 
                            "an_diab_tf" = "Diabetes\nstatus")) + 
  xlab("Covariate") + 
  ggtitle("Variance partition plot for cell proportions") + 
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), 
        text=element_text(size=15), 
        plot.title = element_text(hjust=0))

# Will go forward with individual, sequencing batch (modeling sample GA)
# I will fit the following model based on the variancePartition I ran above
# mixed.lmer <- lmer(cell_prop ~ SampleGA + (1|ANON_ID) + (1|id_run_position), data = cellprop.time.data)

# Model changes in cell prop over time 

# prepare df
cellprop.time.data <- cell_prop %>% 
  select(Mixture, SampleGA, id_run_position, ANON_ID, Sample_taken_at,
         contains("cells") | contains("Macrophage") | contains("phils") | contains("Monocytes"),
         # remove cells that showed no differential abundance over time
         -c(Dendritic.cells.resting, Mast.cells.activated, T.cells.follicular.helper, T.cells.gamma.delta)) %>%
  mutate(ANON_ID = as.factor(ANON_ID),
         Sample_taken_at = as.factor(Sample_taken_at),
         id_run_position = as.factor(id_run_position))

#somewhere to store results
lme.results_cell <- list()
lme_prediction_df_cell <- c()
lme.partial.residuals_cell <- c()

for (cellprop in rownames(var_part_cell_prop)) {
  print(cellprop)
  res <- lmerTest::lmer(as.formula(get(cellprop) ~ splines::bs(x = SampleGA, degree = 3) + (1|ANON_ID) + (1|id_run_position)), 
                        data = cellprop.time.data, 
                        REML = FALSE)
  
  lme.results_cell[[cellprop]] <- res
  
  # based on the model, make a prediction 
  # use ggemmeans to estimate confidence intervals/fits
  prediction <- as.data.frame(ggemmeans(lme.results_cell[[cellprop]], c("SampleGA [all]"))) %>% mutate(Cellprop = cellprop)
  colnames(prediction) <- c("SampleGA", "pred_cellprop_val", "std.err", "conf.low", "conf.high", "group", "Cellprop") 
  lme_prediction_df_cell <- rbind(lme_prediction_df_cell, prediction)
  
  # extract partial residuals 
  # Using the argument partial.residuals = TRUE, what is plotted instead is the observed data with the effects of all the control variables accounted for
  plt <- interact_plot(model = lme.results_cell[[cellprop]], pred = SampleGA, modx = id_run_position, plot.points = TRUE, partial.residuals = TRUE, data = cellprop.time.data, colors = rep("black", length(unique(cellprop.time.data$id_run_position))))
  partial_residuals <- plt$layers[[2]]$data %>% select(`get(cellprop)`)
  colnames(partial_residuals) <- "partial_residuals"
  SampleGA <- plt$layers[[2]]$data$SampleGA
  RNA_sanger_sample_id <- plt$layers[[2]]$data$Mixture
  partial_residuals_plt_dat <- as.data.frame(cbind(cbind(partial_residuals, SampleGA), RNA_sanger_sample_id)) %>% mutate(Cellprop = cellprop)
  lme.partial.residuals_cell <- rbind(lme.partial.residuals_cell, partial_residuals_plt_dat)
  
}

# Based on these models, make time series predictions for each module (Supplementary figure 4)
lme_prediction_df_cell %>% 
  mutate(Cellprop = as.factor(gsub("\\.", " ",Cellprop))) %>% 
  ggplot(aes(x = SampleGA)) + 
  geom_line(aes(SampleGA, pred_cellprop_val), size = 1.7, color = "#32769B") +
  geom_ribbon(aes(SampleGA, pred_cellprop_val, ymin = conf.low, ymax = conf.high), alpha = 0.055, fill = "#32769B") + 
  ylab("Imputed cell proportion") + xlab("Sample gestational age (weeks)") + 
  facet_wrap(~Cellprop, scales = "free", ncol = 3) + 
  ggtitle("Imputed cell proportion values modeled over time") + 
  theme(legend.position = "none", strip.background = element_blank())

# Compare cellprops to Bar cellprops
Bar_all_cellprops <- as.data.frame(rbind(Bar_eos %>% mutate(celltype = "Eosinophils"),
                                         Bar_lym %>% mutate(celltype = "Lymphocytes"),
                                         Bar_mon %>% mutate(celltype = "Monocytes"),
                                         Bar_neu %>% mutate(celltype = "Neutrophils")))
# Plot trajectories in Bar data 
# week: Week of the test relative to delivery (=0).
box_width <- 1
# plot bar cellprops
Bar_cellprop_plot <- Bar_all_cellprops %>%
  select(week, val_5, val_10, val_25, val_50, val_75, val_90, val_95, celltype) %>%
  mutate(
    week_start = as.numeric(str_extract(week, "-?\\d+")),
    week_end = as.numeric(str_extract(week, "(?<=,)-?\\d+")),
    week_mid = (week_start + week_end) / 2 + 0.5,
    gest_age = 40 + week_mid
  ) %>%
  filter(gest_age >= min(sample.info$SampleGA), 
         gest_age <= max(sample.info$SampleGA)) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = gest_age - box_width / 2,
      xmax = gest_age + box_width / 2,
      ymin = val_25/100,
      ymax = val_75/100
    ),
    fill = "white",
    color = "black",
    alpha = 0.6
  ) + 
  # Median line
  geom_segment(
    aes(
      x = gest_age - box_width / 2,
      xend = gest_age + box_width / 2,
      y = val_50/100,
      yend = val_50/100
    ),
    color = "black",
    size = 0.8
  ) +
  labs(
    x = "Gestational age (weeks)",
    y = "Measured cell fraction",
    title = "Bar measured cell fraction across pregnancy"
  ) + geom_smooth(aes(x = gest_age, y = val_50/100), 
                  color = "red",
                  se = FALSE) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  facet_wrap(~celltype, scales = "free", nrow = 1) + 
  theme(plot.title = element_text(size = 13),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 13))

# Re-create cell proportion trajectories for relevant celltypes
pops_cellprop_plts <- lme_prediction_df_cell %>% 
  mutate(Cellprop = ifelse(grepl("B\\.cells\\.|T\\.cells\\.|NK\\.cells", Cellprop), "Lymphocytes", Cellprop)) %>%
  filter(Cellprop %in% c("Eosinophils", "Lymphocytes", "Monocytes", "Neutrophils")) %>%
  group_by(SampleGA, Cellprop) %>%
  summarize(pred_cellprop_val = sum(pred_cellprop_val), 
            std.err = sum(std.err),
            conf.low = sum(conf.low), 
            conf.high = sum(conf.high)) %>%
  ggplot(aes(x = SampleGA)) + 
  geom_line(aes(SampleGA, pred_cellprop_val), size = 1.7, color = "#32769B") +
  geom_ribbon(aes(SampleGA, pred_cellprop_val, ymin = conf.low, ymax = conf.high), alpha = 0.055, fill = "#32769B") + 
  facet_wrap(~Cellprop, scales = "free_y", nrow = 1) + 
  xlab("Gestational age (weeks)") +
  ggtitle("POPs2 imputed cell proportion values modeled over time") +
  ylab("Imputed cell proportion") + 
  theme(plot.title = element_text(size = 13),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 13))


# Update cell proportions collapsing lymphocytes to match bar 
cellprop_for_comp <- cell_prop %>%
  mutate(Lymphocytes = rowSums(
    select(
      .,
      matches("B\\.cells\\.|T\\.cells\\.|NK\\.cells")
    ),
    na.rm = TRUE
  )
  ) %>%
  select(SampleGA, Eosinophils, Lymphocytes, Monocytes, Neutrophils) %>%
  pivot_longer(cols = -SampleGA, names_to = "celltype", values_to = "prop")

# Plot correlation of weekly mean values 
cell_corr_plt <- Bar_all_cellprops %>%
  mutate(
    week_start = as.numeric(str_extract(week, "-?\\d+")),
    week_end = as.numeric(str_extract(week, "(?<=,)-?\\d+")),
    week_mid = (week_start + week_end) / 2 + 0.5,
    gest_age = 40 + week_mid
  ) %>%
  select(gest_age, val_mean, celltype) %>%
  filter(gest_age >= min(sample.info$SampleGA), 
         gest_age <= max(sample.info$SampleGA)) %>%
  rename("val_mean" = "Bar_mean") %>%
  left_join(cellprop_for_comp %>%
              mutate(gest_age = floor(SampleGA)) %>%
              group_by(gest_age, celltype) %>%
              summarize(POPS_mean = mean(prop)) , by = c("gest_age", "celltype")) %>%
  ggplot(aes(x = POPS_mean, y = Bar_mean/100, color = gest_age)) + 
  #geom_jitter(alpha = 0.05, size = 0.1, height = 0.001) + 
  geom_point() +
  facet_wrap(~celltype, nrow = 1, scales = "free") + 
  scale_color_gradient(low = "#BCE4D8", high = "#32769B",
                       breaks = c(12, 20, 28, 36),
                       name = "Gestational\nage (weeks)") + 
  ggpubr::stat_cor(method = "pearson", size = 3) + 
  labs(x = "POPs2 mean imputed cell proportions per week",
       y = "Bar mean measured\ncell proportions per week",
       title = "Correlation of mean cell proportions per week between POPs2 and Bar") + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),     # smaller legend boxes
        legend.text = element_text(size = 7),  # smaller legend labels
        legend.title = element_text(size = 8),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),    # reduce outer plot margin
        legend.spacing.x = unit(0.2, 'cm'),                  # tighten spacing between keys
        legend.spacing.y = unit(0.2, 'cm'),
        plot.title = element_text(size = 13),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 13))

# Plot full validation panel (Supplementary figure 5)
Bar_comp_sup <- gridExtra::grid.arrange(cell_corr_plt + labs(tag = "A"), 
                                        pops_cellprop_plts + labs(tag = "B"), 
                                        Bar_cellprop_plot + labs(tag = "C"),
                                        ncol = 1, nrow = 3, 
                                        widths = c(1), 
                                        heights = c(1.2, 1, 1),
                                        layout_matrix = rbind(c(1),
                                                              c(2),
                                                              c(3)))


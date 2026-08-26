# 02_5_cibersort_continuous_trajectories.R 

################################################################################

# 2.5. Continuous trajectories of imputed cell proportions

################################################################################

# Aim: Plot continuous trajectories of cibersort outputs 

########################### Input paths ###########################
cibersort_inpath <- "rna-seq/analysis/cibersort/outputs/CIBERSORTx_Adjusted.txt"
sample_info_inpath <- "rna-seq/data/sample_covariates_clin_tech.csv"

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

# plot variance parititon 
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

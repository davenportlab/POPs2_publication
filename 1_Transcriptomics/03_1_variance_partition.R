# 03_1_variance_partition.R 

################################################################################

# 3.1. Variance partition 

################################################################################

# Aim: Determine which clinical and technical variables contribute to variation in gene expression 

########################### Output paths ##########################
out_obj_path <- "rna-seq/analysis/dream/outputs/all_covariates_vp.rds"

########################### Input paths ###########################
voom_obj_path <- "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"

########################### Parameters ############################
variance_partition_formula <- "~(1|Sample_taken_at)+(1|id_run_position)+(1|ANON_ID)+(1|an_alcohol_status)+an_est_age+an_scan1_bmi+(1|an_level_edu_str)+(1|an_diab_tf)+pn_monthdel+(1|an_eth_cat_str)+(1|an_smokstat)"

########################### Load packages ###########################
print("loading packages")
library("edgeR")
library("BiocParallel")
library("variancePartition")
suppressPackageStartupMessages(library(tidyverse))
library(stats)
library(lme4)
print("packages loaded")

########################### Load data ###########################
print("Loading data...")

# load voom object
vobjDream <- readRDS(voom_obj_path)

# read covariate data
sample.info <- read.csv(sample_info_inpath)

print("Data loaded.")
Sys.time()

########################### Analysis ###########################
set.seed(1)

print(c("starting at"))
Sys.time()

print(paste("Voom object:", voom_obj_path))
print(paste("Variance partition formula:", variance_partition_formula))
print(paste("Output object:", out_obj_path))

# 1. Format data
sample.info <- sample.info %>%
  # make diabetes t/f
  mutate(an_diab_tf = as.logical(an_diab)) %>% 
  select(RNA_sanger_sample_id, everything()) %>%
  # order by count matrix 
  arrange(factor(RNA_sanger_sample_id, levels = colnames(vobjDream$E)))

rownames(sample.info) <- sample.info$RNA_sanger_sample_id

# check we have 3779 samples and 60678 genes
dim(sample.info)
sample.info[1:5, 1:5]
dim(vobjDream$E)
vobjDream$E[1:5, 1:5]
all(colnames(vobjDream$E) == sample.info$RNA_sanger_sample_id)


# 4. Look at variance partition of all covariates together (takes ~18hrs)
print("start variance partition")
vp <- fitExtractVarPartModel(vobjDream, variance_partition_formula, sample.info)
print("Finished variance partition")
Sys.time()

#write data 
print("Writing out variance partition data")
saveRDS(vp, file = out_obj_path)

print("finished all")
Sys.time()

########################### Plot ###########################
# read in variance partition data
vp <- readRDS(out_obj_path)

# Plot variance partition (Supplementary Figure 23)
plotVarPart(sortCols(vp), , col = c(rep("white", 11), "grey85")) + 
  scale_x_discrete(labels=c("ANON_ID" = "Individual", "id_run_position" = "Sequencing\nbatch", "Sample_taken_at" = "Time-point", "an_eth_cat_str" = "Maternal\nself-reported\nethnicity", "an_scan1_bmi" = "Maternal\nBMI", "pn_monthdel" = "Delivery\nmonth", "an_est_age" = "Maternal\nage", "an_level_edu_str" = "Maternal\neducation\nlevel", "an_alcohol_status" = "Maternal\nalcohol\nstatus", "an_smokstat" = "Maternal\nsmoking\nstatus", "an_diab_tf" = "Maternal\ndiabetes\nstatus", "Residuals" = "Residuals")) + 
  xlab("Covariate") + 
  #theme_minimal() +
  labs(title="Variance partition plot",
       subtitle = paste("Covariates tested for dream differential gene expression model")) + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks = element_blank()) 
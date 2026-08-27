# 03_6_run_outcome_vp.R 

################################################################################

# 3.6. Run variance partition of pregnancy outcome groups

################################################################################

# Aim: Identify the contribution of pregnancy outcomes to gene expression variance

########################### Output paths ##########################
out_obj_path <- "rna-seq/analysis/outcome_pred/outputs/de_all_tps/all_clinical_outcome_covariates"

########################### Input paths ###########################
voom_obj_path <- "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
# sample info input 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"

########################### Parameters ############################
variance_partition_formula <- "~(1|Sample_taken_at)+(1|id_run_position)+(1|ANON_ID)+(1|pn_gdm)+(1|sga)+(1|pn_petACOG13)+(1|pn_ptd_sp)"

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
# load voom object
vobjDream <- readRDS(voom_obj_path)
# read in metadata 
sample.info <- read.csv(sample_info_inpath) 

########################### Analysis ###########################

set.seed(1)

print(c("starting at"))
Sys.time()

print(paste("Voom object:", voom_obj_path))
print(paste("Variance partition formula:", variance_partition_formula))
print(paste("Output object:", out_obj_path))

# 1. load data 
print("Loading data...")

# format covariate data
sample.info <- sample.info %>%
  mutate(pn_gdm = as.factor(pn_diabetes_3cat == 2),
         sga = as.factor(BW_Centile_Br1990 < 10),
         pn_petACOG13 = as.factor(pn_petACOG13), 
         pn_ptd_sp = as.factor(pn_ptd_sp)) %>% 
  dplyr::select(RNA_sanger_sample_id, everything()) %>%
  # order by count matrix 
  arrange(factor(RNA_sanger_sample_id, levels = colnames(vobjDream$E)))

rownames(sample.info) <- sample.info$RNA_sanger_sample_id

# check 
dim(sample.info)
sample.info[1:5, 1:5]
dim(vobjDream$E)
vobjDream$E[1:5, 1:5]
all(colnames(vobjDream$E) == sample.info$RNA_sanger_sample_id)
print("Data loaded.")
Sys.time()

# 4. Look at variance partition of all covariates together 
print("start variance partition")
vp <- fitExtractVarPartModel(vobjDream, variance_partition_formula, sample.info)
print("Finished variance partition")
Sys.time()

#write data 
print("Writing out variance partition data")
saveRDS(vp, file = paste0(out_obj_path, "_vp.rds"))

print("finished all")
Sys.time()

########################### Plot ###########################
vp <- readRDS(paste0(out_obj_path, "_vp.rds"))

plotVarPart(sortCols(vp), , col = c(rep("white", 11), "grey85")) + 
  scale_x_discrete(labels=c("ANON_ID" = "Individual", 
                            "id_run_position" = "Sequencing\nbatch", 
                            "Sample_taken_at" = "Time-point", 
                            "an_eth_cat_str" = "Maternal\nself-reported\nethnicity", 
                            "pn_gdm" = "Gestational\ndiabetes", 
                            "pn_petACOG13" = "Pre-eclampsia", 
                            "sga" = "Small for\ngestational age", 
                            "pn_ptd_sp" = "Spontaneous\npreterm\nbirth", 
                            "Residuals" = "Residuals")) + 
  xlab("Covariate") + 
  labs(title="Variance partition plot",
       subtitle = paste("Covariates used in dream differential gene expression model and clinical covariates")) + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks = element_blank()) 


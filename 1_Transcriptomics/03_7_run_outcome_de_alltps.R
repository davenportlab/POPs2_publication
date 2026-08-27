# 03_7_run_outcome_de_alltps.R 

################################################################################

# 3.7. Run differential expression across pregnancy outcome groups

################################################################################

# Aim: Identify genes differentially expressed across pregnancy outcome groups 
# Note: This script runs DE on samples from all time-points. 
#       In this example script, we are only running DE on the GDM outcome. 
#       However, the models for PE, SGA and sPTB can all be run in the same way,
#       so are exemplified in comments in this script.

########################### Output paths ##########################
out_obj_path_gdm <- "rna-seq/analysis/outcome_pred/outputs/de_all_tps/GDM"
# out_obj_path_pe <-  "rna-seq/analysis/outcome_pred/outputs/de_all_tps/PE"
# out_obj_path_sga <- "rna-seq/analysis/outcome_pred/outputs/de_all_tps/SGA"
# out_obj_path_ptd <- "rna-seq/analysis/outcome_pred/outputs/de_all_tps/PTD"

########################### Input paths ###########################
voom_obj_path <-  "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
# sample info input 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"

########################### Parameters ############################
# DE formula
form_gdm <- "~0+pn_gdm+Sample_taken_at+(1|ANON_ID)+(1|id_run_position)"
# form_pe <- "~0+pn_petACOG13+Sample_taken_at+(1|ANON_ID)+(1|id_run_position)"
# form_sga <- "~0+sga+Sample_taken_at+(1|ANON_ID)+(1|id_run_position)"
# form_ptd <- "~0+pn_ptd_sp+Sample_taken_at+(1|ANON_ID)+(1|id_run_position)"

# VP formula
variance_partition_formula_gdm <- "~(1|Sample_taken_at)+(1|ANON_ID)+(1|id_run_position)+(1|pn_gdm)"
# variance_partition_formula_pe <- "~(1|Sample_taken_at)+(1|ANON_ID)+(1|id_run_position)+(1|pn_petACOG13)"
# variance_partition_formula_sga <- "~(1|Sample_taken_at)+(1|ANON_ID)+(1|id_run_position)+(1|sga)"
# variance_partition_formula_ptd <- "~(1|Sample_taken_at)+(1|ANON_ID)+(1|id_run_position)+(1|pn_ptd_sp)"

########################### Load packages ###########################
print("loading packages")
library(edgeR)
library(BiocParallel)
library(variancePartition)
suppressPackageStartupMessages(library(tidyverse))
library(stats)
library(lme4)
print("packages loaded")

########################### Load data ###########################
print("Loading data...")
# load voom object
vobjDream <- readRDS(voom_obj_path)
# read in metadata 
sample.info <- read.csv(sample_info_inpath)
print("Data loaded.")
Sys.time()

########################### Analysis ###########################

set.seed(1)

print(c("starting at"))
Sys.time()

sample.info <- sample.info %>%
  mutate(pn_gdm = as.character(as.factor(pn_diabetes_3cat == 2)),
         sga = as.character(as.factor(BW_Centile_Br1990 < 10)),
         pn_petACOG13 = as.character(as.factor(pn_petACOG13 == 1)), 
         pn_ptd_sp = as.character(as.factor(pn_ptd_sp == 1))) %>% 
  dplyr::select(RNA_sanger_sample_id, Sample_taken_at, ANON_ID, id_run_position,
         pn_gdm, sga, pn_petACOG13, pn_ptd_sp) %>%
  arrange(factor(RNA_sanger_sample_id, levels = colnames(vobjDream$E))) 

# add covariate 
rownames(sample.info) <- sample.info$RNA_sanger_sample_id
# check we have 3779 samples and 60678 genes
dim(sample.info)
dim(vobjDream$E)
# check things are in correct order
all(rownames(sample.info) == colnames(vobjDream$E))


# 2. make contrasts 
# gdm
L <- makeContrastsDream(formula = form_gdm, 
                            data = sample.info,
                            contrasts = c(compare_GDM_vs_noGDM = "pn_gdmFALSE - pn_gdmTRUE"))
# # sga
# L <- makeContrastsDream(formula = form_sga, 
#                             data = clin_tech_var,
#                             contrasts = c(compare_SGA_vs_noSGA = "sgaFALSE - sgaTRUE"))
# # pe
# L <- makeContrastsDream(formula = form_pe, 
#                            data = clin_tech_var,
#                            contrasts = c(compare_PE_vs_noPE = "pn_petACOG13FALSE - pn_petACOG13TRUE"))
# # pn_ptd_sp
# L <- makeContrastsDream(formula = form_ptd, 
#                             data = clin_tech_var,
#                             contrasts = c(compare_PTD_vs_noPTD = "pn_ptd_spFALSE - pn_ptd_spTRUE"))

table(sample.info$pn_gdm)


# 3. fit dream model with contrasts
print("start fit outcome dream model with contrasts")
Sys.time()
fit <- dream(vobjDream, form_gdm, sample.info, L)
fit <- eBayes(fit)
print("end fit dream model with contrasts")
Sys.time()

#write data 
print("Writing out model fit data")
saveRDS(fit, file = paste0(out_obj_path_gdm, "_dream_fit.rds"))

# 4. Look at variance partition of modeled covariates 
print("start variance partition")
vp <- fitExtractVarPartModel(vobjDream, variance_partition_formula_gdm, sample.info)
print("Finished variance partition")
Sys.time()

#write data 
print("Writing out variance partition data")
saveRDS(vp, file = paste0(out_obj_path_gdm, "_vp.rds"))

print("finished all modeling tasks")
Sys.time()

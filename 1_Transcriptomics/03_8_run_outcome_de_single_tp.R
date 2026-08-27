# 03_8_run_outcome_de_single_tp.R 

################################################################################

# 3.8. Run differential expression across pregnancy outcome groups

################################################################################

# Aim: Identify genes differentially expressed across pregnancy outcome groups 
# Note: This script runs DE on samples from single time-points. 
#       In this example script, we are only running DE on samples from the 12 week time-point. 
#       However, the models for 20, 28 and 36 week samples can all be run in the same way,
#       so are exemplified in comments in this script.

########################### Output paths ##########################
# Need to update for other time-points
out_obj_path_gdm <- "./rna-seq/analysis/outcome_pred/outputs/de_12wk/GDM_12wk"
out_obj_path_sga <- "./rna-seq/analysis/outcome_pred/outputs/de_12wk/SGA_12wk"
out_obj_path_pe <- "./rna-seq/analysis/outcome_pred/outputs/de_12wk/PE_12wk"
out_obj_path_ptd <- "./rna-seq/analysis/outcome_pred/outputs/de_12wk/PTD_12wk"

########################### Input paths ###########################
voom_obj_path <-  "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
# sample info input 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Parameters ############################
# time-point 
tp <- "12_weeks"
# gdm
form_gdm <- "~0+pn_gdm+(1|id_run_position)"
variance_partition_formula_gdm <- "~(1|id_run_position)+(1|pn_gdm)"
# sga
form_sga <- "~0+sga+(1|id_run_position)"
variance_partition_formula_sga <- "~(1|id_run_position)+(1|sga)"
# pn_petACOG13 (pe)
form_pe <- "~0+pn_petACOG13+(1|id_run_position)"
variance_partition_formula_pe <- "~(1|id_run_position)+(1|pn_petACOG13)"
# pn_ptd_sp
form_ptd <- "~0+pn_ptd_sp+(1|id_run_position)"
variance_partition_formula_ptd <- "~(1|id_run_position)+(1|pn_ptd_sp)"

# define deg cutoffs
p_val = 0.01
log_FC = log2(1.5)

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
# read covariate data
sample.info <- read.csv(sample_info_inpath)
print("Data loaded.")
Sys.time()

########################### Analysis ###########################

set.seed(1)

print(c("starting at"))
Sys.time()

miscarried_samples_to_remove <- c("pops212970284", "pops212971619")

# format data
sample.info <- sample.info %>%
  mutate(pn_gdm = as.character(as.factor(pn_diabetes_3cat == 2)),
         sga = as.character(as.factor(BW_Centile_Br1990 < 10)),
         pn_petACOG13 = as.character(as.factor(pn_petACOG13 == 1)), 
         pn_ptd_sp = as.character(as.factor(pn_ptd_sp == 1)),
         # set the outcome variables to NA for all samples that are not at specified time point (tp)
         # and misscarried samples as a way to filter all non-tp samples 
         # from analysis without subsetting the vobj
         pn_gdm = ifelse((Sample_taken_at != tp | 
                            RNA_sanger_sample_id %in% miscarried_samples_to_remove), 
                         NA, pn_gdm),
         sga = ifelse((Sample_taken_at != tp | 
                         RNA_sanger_sample_id %in% miscarried_samples_to_remove), 
                      NA, sga),
         pn_petACOG13 = ifelse((Sample_taken_at != tp | 
                                  RNA_sanger_sample_id %in% miscarried_samples_to_remove), 
                               NA, pn_petACOG13),
         pn_ptd_sp = ifelse((Sample_taken_at != tp | 
                               RNA_sanger_sample_id %in% miscarried_samples_to_remove), 
                            NA, pn_ptd_sp)) %>% 
  dplyr::select(RNA_sanger_sample_id, id_run_position,
         pn_gdm, sga, pn_petACOG13, pn_ptd_sp) %>%
  arrange(factor(RNA_sanger_sample_id, levels = colnames(vobjDream$E))) 


# add covariate 
rownames(sample.info) <- sample.info$RNA_sanger_sample_id
dim(sample.info)
dim(vobjDream$E)
# check things are in correct order
all(rownames(sample.info) == colnames(vobjDream$E))

# 2. make contrasts 
# gdm
L_gdm <- makeContrastsDream(formula = form_gdm, 
                            data = sample.info,
                            contrasts = c(compare_GDM_vs_noGDM = "pn_gdmFALSE - pn_gdmTRUE"))
table(sample.info$pn_gdm)
plotContrasts(L_gdm)
# sga
L_sga <- makeContrastsDream(formula = form_sga, 
                            data = sample.info,
                            contrasts = c(compare_SGA_vs_noSGA = "sgaFALSE - sgaTRUE"))
table(sample.info$sga)
plotContrasts(L_sga)
# pn_petACOG13
L_pe <- makeContrastsDream(formula = form_pe, 
                           data = sample.info,
                           contrasts = c(compare_PE_vs_noPE = "pn_petACOG13FALSE - pn_petACOG13TRUE"))
table(sample.info$pn_petACOG13)
plotContrasts(L_pe)
# pn_ptd_sp
L_ptd <- makeContrastsDream(formula = form_ptd, 
                            data = sample.info,
                            contrasts = c(compare_PTD_vs_noPTD = "pn_ptd_spFALSE - pn_ptd_spTRUE"))
table(sample.info$pn_ptd_sp)
plotContrasts(L_ptd)


# 3. fit dream model with contrasts
print("start fit gdm dream model with contrasts")
Sys.time()
fit_gdm <- dream(vobjDream, form_gdm, sample.info, L_gdm)
fit_gdm <- eBayes(fit_gdm)
print("end fit gdm dream model with contrasts")
print("start fit sga dream model with contrasts")
Sys.time()
fit_sga <- dream(vobjDream, form_sga, sample.info, L_sga)
fit_sga <- eBayes(fit_sga)
print("end fit sga dream model with contrasts")
print("start fit pe dream model with contrasts")
Sys.time()
fit_pe <- dream(vobjDream, form_pe, sample.info, L_pe)
fit_pe <- eBayes(fit_pe)
print("end fit pe dream model with contrasts")
print("start fit ptd dream model with contrasts")
Sys.time()
fit_ptd <- dream(vobjDream, form_ptd, sample.info, L_ptd)
fit_ptd <- eBayes(fit_ptd)
print("end fit ptd dream model with contrasts")
Sys.time()
# get names of available coefficients and contrasts for testing
colnames(fit_gdm)
colnames(fit_sga)
colnames(fit_pe)
colnames(fit_ptd)
print("Finished dream model")
Sys.time()

#write data 
print("Writing out model fit data")
saveRDS(fit_gdm, file = paste0(out_obj_path_gdm, "_dream_fit.rds"))
saveRDS(fit_sga, file = paste0(out_obj_path_sga, "_dream_fit.rds"))
saveRDS(fit_pe, file = paste0(out_obj_path_pe, "_dream_fit.rds"))
saveRDS(fit_ptd, file = paste0(out_obj_path_ptd, "_dream_fit.rds"))

# 4. Look at variance partition of all covariates together 
print("start variance partition gdm")
vp_gdm <- fitExtractVarPartModel(vobjDream, variance_partition_formula_gdm, sample.info)
print("Finished variance partition gdm")
Sys.time()
print("start variance partition sga")
vp_sga <- fitExtractVarPartModel(vobjDream, variance_partition_formula_sga, sample.info)
print("Finished variance partition sga")
Sys.time()
print("start variance partition pe")
vp_pe <- fitExtractVarPartModel(vobjDream, variance_partition_formula_pe, sample.info)
print("Finished variance partition pe")
Sys.time()
print("start variance partition ptd")
vp_ptd <- fitExtractVarPartModel(vobjDream, variance_partition_formula_ptd, sample.info)
print("Finished variance partition ptd")
Sys.time()
print("Finished variance partition")

#write data 
print("Writing out variance partition data")
saveRDS(vp_gdm, file = paste0(out_obj_path_gdm, "_vp.rds"))
saveRDS(vp_sga, file = paste0(out_obj_path_sga, "_vp.rds"))
saveRDS(vp_pe, file = paste0(out_obj_path_pe, "_vp.rds"))
saveRDS(vp_ptd, file = paste0(out_obj_path_ptd, "_vp.rds"))

print("finished all")
Sys.time()



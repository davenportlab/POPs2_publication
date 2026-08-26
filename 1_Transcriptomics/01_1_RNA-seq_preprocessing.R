# 01.1_RNA-seq_preprocessing.R 

################################################################################

# 1.1. RNA-seq_preprocessing

################################################################################

# Aim: generate log2cpm count matrix for following analyses (including gene filtering)

########################### Output paths ##########################
# voom output file (.rds)
voom_outpath = "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
# cpm count matrix output file (.txt)
cpm_counts_outpath = "rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt"

########################### Input paths ###########################
# raw count matrix as input
raw_count_matrix_inpath = "rna-seq/data/star-fc-genecounts_3_filt_samples_HLApm.txt"
# sample info input 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"

########################### Parameters ############################
# voom formula (we use ~1 to make this count matrix generalizeable for future analyses)
form <- ~ 1

set.seed(1)

print(c("starting at"))
Sys.time()

print("Loading packages...")
########################### Load packages ###########################
library(tidyverse)
library(edgeR)
library(variancePartition)
print("Packages loaded.")

print("Loading data...")
# load data 
counts <- read.table(raw_count_matrix_inpath, row.names = 1, header = TRUE) 
sample.info <- read.csv(sample_info_inpath) %>%
  # order by count matrix 
  arrange(factor(RNA_sanger_sample_id, levels = colnames(counts)))

print("Data loaded.")

print("Filtering genes...")
##### filter genes based on cpm ##### 
# filter to genes with 10 reads per library size in at least 5% of the samples 
# calculate average library size
lib_sizes <- counts %>% colSums()
mean_lib_size <- mean(lib_sizes) # 24100407
# define cpm threshold: Count of 10 from 24 million reads 10/24.100407 = 0.4149307
cpm_threshold <- 10/(mean_lib_size/1000000)
# calculate cpm 
cpm_counts <- cpm(counts)
# use cpm threshold of 0.4149307 for avg of 24 million reads
keep <- filterByExpr(cpm_counts, min.count = cpm_threshold, min.prop = 0.05) # could increase min prop to 0.1 to be more stringent
cpm_counts_filtered <- cpm_counts[keep,]
# how many genes to keep
dim(cpm_counts_filtered)

# now we have generated the list of genes we will keep, but will go back and filter on the raw count matrix 
# make list of genes to keep 
genes_to_keep <- rownames(cpm_counts_filtered)
# make original/raw count matrix into DGElist object 
d0 <- DGEList(counts)
# actually filter to the genes we will keep, as determined above
d1 <- d0[genes_to_keep,]
dim(d1)
print("Genes filtered")

print("Calculating normalization factors...")
# calculate TMM normalization factors
d1 <- calcNormFactors(d1)
print("TMM normalization factors calculated.")

print("Running voom...")
param <- SnowParam(4, "SOCK", progressbar = TRUE)
# run voom WITHOUT formula
print(paste(c("Formula:", form)))
vobjDream_noform <- voomWithDreamWeights(d1, form, sample.info, BPPARAM = param, plot=TRUE, save.plot=TRUE)
print("Finished running voom.")
# write out voom object for use in dream 
saveRDS(vobjDream_noform, file = voom_outpath)
print("Finished saving voom output object.")

# write out cpm count matrix for use in plotting and cibersort
print("Writing out cpm count matrix.")
cpm_counts <- vobjDream_noform$E
# check this is what you expect
cpm_counts[1:5, 1:5]
dim(cpm_counts)
cpm_counts %>% write.table(file = cpm_counts_outpath, 
                           quote = FALSE, sep = "\t")
print("Finished saving cpm count matrix output object.")

# check dimensions match the count matrix it was made from
cpm_counts <- read.table(cpm_counts_outpath, header = TRUE)
cpm_counts[1:5, 1:5]
dim(cpm_counts)

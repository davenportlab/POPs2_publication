# 02_1_cibersort_prep.R 

################################################################################

# 2.1. Prepare data for CIBERSORTx

################################################################################

# Aim: prepare the data in the correct format for cibersort input 

########################### Output paths ##########################
# cibersort input matrix 
cibersort_input_count_matrix_outpath <- "rna-seq/analysis/cibersort/outputs/cibersort_input_counts.tsv"

########################### Input paths ###########################
# cpm count matrix input (made in RNA-seq_preprocessing.R)
cpm_counts_inpath <- "rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt"
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Load packages ###########################
library(edgeR)
library(tidyverse)
library(statmod)
library(rtracklayer)

# read in log2cpm counts (made with voom)
counts <- read.table(cpm_counts_inpath, header = TRUE) %>%
  rownames_to_column(var = "ENSEMBL_ID")
dim(counts) # [1] 18826  3777
counts[1:5, 1:5]

# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
# edit since we added two genes that are not in the reference 
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

# Change gene IDs to gene names
counts_named <- counts %>% left_join(id_name, by = c("ENSEMBL_ID" = "gtf.gene_id")) %>%
  # do something about the HLAs 
  mutate(gtf.gene_name = ifelse(is.na(gtf.gene_name), ENSEMBL_ID, gtf.gene_name))

# Sum non-unique gene names to collapse into one 
counts_named_compact <- counts_named %>% 
  select(-ENSEMBL_ID) %>%
  group_by(gtf.gene_name) %>%
  summarize(across(everything(), sum), 
            .groups = 'drop') 
# check this looks as expected
counts_named_compact[1:5, 1:5]
dim(counts_named_compact)

# make the matrix collapsed by gene names into a df in the correct format
counts_named_compact <- data.frame(counts_named_compact)
rownames(counts_named_compact) <- counts_named_compact$gtf.gene_name
counts_named_compact[,1] <- NULL
counts_named_compact <- cbind(Genes = rownames(counts_named_compact), counts_named_compact)
counts_named_compact <- as.tibble(counts_named_compact)
# check this looks as expected
counts_named_compact[1:5, 1:5]
dim(counts_named_compact)

# write out this matrix for use as cibersort input 
write.table(counts_named_compact, cibersort_input_count_matrix_outpath,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

# test reading back in 
counts_named_compact <- read.table(cibersort_input_count_matrix_outpath, header = TRUE)
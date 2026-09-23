# 04_2_prepare_subsetted_interval_cohort.R 

################################################################################

# 4.2. Prepare subsetted interval cohort for eQTL mapping

################################################################################

# Aim: Subset interval cohort to females of reproductive age and preprocess data

########################### Input paths ###########################
interval_path <- "projects_v2/interval_rna/" # 
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Load packages ###########################
library(tidyverse)
library(edgeR)
library(variancePartition)

# set ggplot theme 
theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Analysis ###########################

# interval metadata
sexcheck <- read.table(paste0(interval_path, "recovery_intermediate_files/INTERVAL_RNAseq/genotype/chrX/INTERVAL_chrX_sexcheck.sexcheck"), 
                       header = TRUE)
interval_metadata <- read.delim(paste0(interval_path, "recovery_intermediate_files/INTERVAL_RNAseq/covariate/INTERVAL_RNAseq_phase1-2_fullcovariates.txt"), header = TRUE) %>%
  select(affymetrix_ID, 
         sample_id, 
         paste0("sequencingBatch_", 1:15), 
         age_RNA,
         sex,
         paste0("PC", 1:10)
  ) %>%
  pivot_longer(
    cols = starts_with("sequencingBatch_"),
    names_to = "sequencingBatch",
    values_to = "inBatch"
  ) %>%
  filter(inBatch == 1) %>%
  mutate(
    sequencingBatch = readr::parse_number(sequencingBatch) # extract 1–15
  ) %>%
  select(-inBatch) %>%
  mutate(sex = ifelse(sex == 0, 2, 1)) %>%
  left_join(sexcheck %>% # join to sexcheck to confirm 
              select(IID, SNPSEX), by = c("affymetrix_ID" = "IID")) # SNPSEX	Imputed sex code (1 = male, 2 = female, 0 = unknown)
# interval count matrix (will re-read in raw count matrix when we have access)
interval_counts_raw <- read.table(paste0(interval_path, 
                                         "recovery_intermediate_files/INTERVAL_RNAseq/expression_data/raw/5591-star-fc-genecounts.txt"), 
                                  sep = "\t", header = TRUE)
# gtf
gtf <- as.data.frame(rtracklayer::import(gtf_inpath))
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

################################################################################
# Subset to females < 45yo
################################################################################

# check that all sex in metadata match plink sexcheck
nrow(interval_metadata %>% filter(sex != SNPSEX)) == 0

# plot distribution of age by sex 
interval_metadata %>%
  mutate(sex = as.factor(sex)) %>%
  ggplot(aes(x = age_RNA, group = sex, color = sex)) + 
  geom_density() + 
  geom_vline(xintercept = 45) + 
  ggtitle("Age distribution by sex in Interval")

# filter interval metadata to females less than 45 
interval_metadata_using <- interval_metadata %>%
  filter(sex == 2, 
         age_RNA < 45)

# Plot the age distribution of the samples in the females of reproductive age subset 
# (Supplementary figure 33A)
interval_metadata_using %>%
  ggplot(aes(x = age_RNA)) + 
  geom_histogram(fill = "black") + 
  ggtitle("Age distribution") + 
  xlab("Age")
summary(interval_metadata$age_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.30   50.50   58.10   55.83   64.90   79.10 

dim(interval_metadata_using) # n = 445

# subset the gene expression matrix to these samples 
dim(interval_counts_raw)
# interval_counts_using <- interval_counts %>% 
#   filter(sample_id %in% interval_metadata_using$sample_id) %>%
#   column_to_rownames(var = "sample_id") %>%
#   t()
interval_counts_using <- interval_counts_raw %>% 
  select(ENSEMBL_ID, interval_metadata_using$sample_id) %>% 
  column_to_rownames(var = "ENSEMBL_ID")
dim(interval_counts_using) 

# write out new count matrix
interval_counts_using %>% write.table(file = paste0(interval_path, 
                                                    "POPS2_comparison/data_subset/Interval_counts_female_less45_raw.txt"), 
                                      quote = FALSE, sep = "\t")
# write out new metadata file
interval_metadata_using %>% select(-SNPSEX) %>%
  write.csv(file = paste0(interval_path, 
                          "POPS2_comparison/data_subset/Interval_metadata_female_less45.csv"), 
            row.names = FALSE, quote = FALSE)

# write out list of samples we will use so we can subset the interval files 
cbind(interval_metadata_using$affymetrix_ID, interval_metadata_using$affymetrix_ID) %>%
  write.table(file = paste0(interval_path,
                            "POPS2_comparison/data_subset/Interval_female_less45_geno_sampleIDs.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)

# ################################################################################
# # RNA-seq pre-processing
# ################################################################################
# 
# Re-process the RNA-seq samples as I did in the POPS cohort 
# I will proceed with them as they are
# calculate average library size
lib_sizes <- interval_counts_using %>% colSums()
mean_lib_size <- mean(lib_sizes) # 28391464 
# define cpm threshold: Count of 10 from 28 million reads 10/28.39146 = 0.3522186
cpm_threshold <- 10/(mean_lib_size/1000000)
# calculate cpm
interval_cpm_counts <- cpm(interval_counts_using)
# use cpm threshold of 0.3522186 for avg of 28 million reads
keep <- filterByExpr(interval_cpm_counts, min.count = cpm_threshold, min.prop = 0.05)
interval_cpm_counts_filtered <- interval_cpm_counts[keep,]
# how many genes to keep
dim(interval_cpm_counts_filtered) # 20132 - more than POPs2...

# now we have generated the list of genes we will keep, but will go back and filter on the raw count matrix
# make list of genes to keep
genes_to_keep <- rownames(interval_cpm_counts_filtered)
# make original/raw count matrix into DGElist object
d0 <- DGEList(interval_counts_using)
# actually filter to the genes we will keep, as determined above
d1 <- d0[genes_to_keep,]
dim(d1)
print("Genes filtered")

print("Calculating normalization factors...")
# calculate TMM normalization factors
d1 <- calcNormFactors(d1)
print("TMM normalization factors calculated.")

interval_metadata_using <- interval_metadata_using %>%
  # order by count matrix
  arrange(factor(sample_id, levels = colnames(interval_counts_using)))

form = ~ 1
print("Running voom...")
param <- SnowParam(4, "SOCK", progressbar = TRUE)
# run voom WITHOUT formula
print(paste(c("Formula:", form)))
vobjDream_noform <- voomWithDreamWeights(d1, form, interval_metadata_using %>% column_to_rownames(var = "sample_id"), BPPARAM = param, plot=TRUE, save.plot=TRUE)
print("Finished running voom.")

# write out cpm count matrix for use in plotting and cibersort
print("Writing out cpm count matrix.")
cpm_counts <- vobjDream_noform$E
# check this is what you expect
cpm_counts[1:5, 1:5]
dim(cpm_counts)
cpm_counts %>% write.table(file = paste0(interval_path, 
                                         "POPS2_comparison/data_subset/Interval_counts_female_less45_log2cpm.txt"),
                           quote = FALSE, sep = "\t")
print("Finished saving cpm count matrix output object.")

# check all counts are between 0-15
max(cpm_counts) # 16.05587
min(cpm_counts) # -7.692074

################################################################################
# Cibersort 
################################################################################
# read in the cpm count matrix
interval_cpm_counts_using <- read.table(paste0(interval_path, 
                                               "POPS2_comparison/data_subset/Interval_counts_female_less45_log2cpm.txt"), 
                                        sep = "\t", header = TRUE)

# Change gene IDs to gene names
counts_named <- interval_cpm_counts_using %>% 
  as.data.frame() %>%
  # make gene names a column
  rownames_to_column(var = "gtf.gene_id") %>%
  # join to gene names
  left_join(id_name, by = c("gtf.gene_id")) %>%
  # do something about the HLAs 
  mutate(gtf.gene_name = ifelse(is.na(gtf.gene_name), gtf.gene_id, gtf.gene_name))

# Sum non-unique gene names to collapse into one 
counts_named_compact <- counts_named %>% 
  select(-gtf.gene_id) %>%
  group_by(gtf.gene_name) %>%
  summarize(across(everything(), sum), 
            .groups = 'drop') 
# check this looks as expected
counts_named_compact[1:5, 1:5]
dim(counts_named_compact)

# make the matrix collapsed by gene names into a df in the correct format
counts_named_compact <- counts_named_compact %>% rename("Genes" = "gtf.gene_name")
# check this looks as expected
counts_named_compact[1:5, 1:5]
dim(counts_named_compact)

# write out this matrix for use as cibersort input 
write.table(counts_named_compact, paste0(interval_path, "POPS2_comparison/cibersort/Input/Interval_cibersort_input_counts_myprocess.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

################################################################################
#                     RUN IN BASH
################################################################################
# run cibersort 
# # start interactive job
# fash 20
# # load module
# module load HGI/common/Cibersortx/2020-04-04
# # cd to correct dir 
# cd interval_rna/POPS2_comparison/cibersort/
# # run cibersort on docker 
# run_cibersortx_fractions --input-dir ./Input --output-dir ./Output_myprocess --username sh50@sanger.ac.uk --token <your token> --mixture Interval_cibersort_input_counts_myprocess.tsv --sigmatrix LM22.txt --rmbatchBmode TRUE --perm 100 --verbose TRUE --QN FALSE
################################################################################

# read in the results 
interval_cibersort <- read.delim(paste0(interval_path, "POPS2_comparison/cibersort/Output_myprocess/CIBERSORTx_Adjusted.txt"))

# Plot neutrophil and monocyte distributions centered at zero for eQTL model
# (Supplementary Figure 33B)
interval_cibersort %>% select(Mixture, Neutrophils, Monocytes) %>%
  pivot_longer(-Mixture, names_to = "Cell_type", values_to = "Cell_prop") %>%
  group_by(Cell_type) %>%
  mutate(scaled_cell_prop = scale(Cell_prop, center = TRUE, scale = FALSE)) %>%
  ggplot(aes(x = scaled_cell_prop)) + 
  geom_histogram(fill = "black") + 
  facet_wrap(~Cell_type, scales = "free", ncol = 1) + 
  ggtitle("Centered imputed cell proportions") + 
  xlab("Centered imputed cell proportion")
scaled_cellprop_plt

################################################################################
# Subset genotyping data to females under 45
################################################################################

################################################################################
#                     RUN IN BASH
################################################################################
# # start interactive job
# fash 20
# # load module with plink 
# module load HGI/softpack/groups/trynka/popgen/1.0
# # cd to correct dir
# cd interval_rna/recovery_intermediate_files/INTERVAL_RNAseq/genotype
# 
# # subset autosomes to FRA samples
# plink \
# --bfile INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_AllAutosomes \
# --keep ../../../POPS2_comparison/data_subset/Interval_female_less45_geno_sampleIDs.txt \
# --make-bed \
# --out ../../../POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllAutosomes_female_u45
# 
# # subset x chromosome to FRA samples
# plink \
# --bfile INTERVAL_RNAseq_imputed_b38_biallelic_MAF0.005_chrX \
# --keep ../../../POPS2_comparison/data_subset/Interval_female_less45_geno_sampleIDs.txt \
# --make-bed \
# --out ../../../POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_chrX_female_u45
# 
# cd interval_rna/POPS2_comparison/data_subset/genotyping
# 
# # change the FID in autosomes from 0 to the same as the IID to match the X chr
# awk '{print $2, $2, $3, $4, $5, $6}' INTERVAL_RNAseq_imputed_b38_AllAutosomes_female_u45.fam > tmp.fam
# mv tmp.fam INTERVAL_RNAseq_imputed_b38_AllAutosomes_female_u45.fam
# 
# # Merge the autosomes and X chromosome
# plink \
# --bfile INTERVAL_RNAseq_imputed_b38_AllAutosomes_female_u45 \
# --bmerge INTERVAL_RNAseq_imputed_b38_chrX_female_u45 \
# --make-bed \
# --out INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45
# 
# # Filter variants in all chromosomes (to maf 0.01)
# plink \
# --bfile INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45 \
# --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed \
# --out INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt
#
################################################################################

# Use PCs from full cohort (to capture the most axes of genetic variation)
PC_plt1 <- interval_metadata_using %>%
  as.data.frame() %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 0.5) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01)  # fixed decimals
  )
PC_plt2 <- interval_metadata_using %>%
  as.data.frame() %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(size = 0.5) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01)  # fixed decimals
  )
PC_plt3 <- interval_metadata_using %>%
  as.data.frame() %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point(size = 0.5) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01)  # fixed decimals
  )

# Plot all Interval genotype PCs
# (Supplementary Figure 33C)
PCA_plt <- gridExtra::grid.arrange(PC_plt1 + labs(tag = "C"),
                                   PC_plt2 + labs(tag = " "),
                                   PC_plt3 + labs(tag = " "),
                                   ncol = 3, nrow = 1, 
                                   widths = c(1, 1, 1), 
                                   heights = c(1),
                                   layout_matrix = rbind(c(1, 2, 3)),
                                   top = "Top 6 genotyping PCs")

################################################################################
# Calculate PEER factors
################################################################################

interval_cibersort <- read.delim(paste0(interval_path, "POPS2_comparison/cibersort/Output_myprocess/CIBERSORTx_Adjusted.txt"))

# Check we have genotyping sample for all RNA-seq samples 
samples <- read.table(paste0(interval_path, "POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt.fam"))
setequal(samples$V1, interval_metadata_using$affymetrix_ID)

# make covariates combined matrix
# centered neutrophils, centered monocytes, genotype PCs 1-4
covariates <- interval_metadata_using %>% 
  select(affymetrix_ID, sample_id, PC1, PC2, PC3, PC4) %>%
  # join to cell props
  left_join(interval_cibersort %>%
              mutate(Neutrophils_centered = scale(Neutrophils, center = TRUE, scale = FALSE),
                     Monocytes_centered = scale(Monocytes, center = TRUE, scale = FALSE)) %>%
              select(Mixture, Neutrophils_centered, Monocytes_centered),
            by = c("sample_id" = "Mixture")) %>%
  # re-name genotype PCs
  rename("geno_PC1" = "PC1", "geno_PC2" = "PC2", "geno_PC3" = "PC3", "geno_PC4" = "PC4") %>%
  select(-affymetrix_ID) %>%
  column_to_rownames("sample_id")

# show cell props were centered to mean
covariates %>% ggplot(aes(x = Neutrophils_centered)) + geom_histogram()
covariates %>% ggplot(aes(x = Monocytes_centered)) + geom_histogram()

interval_cpm_counts_using <- read.table(paste0(interval_path, 
                                               "POPS2_comparison/data_subset/Interval_counts_female_less45_log2cpm.txt"), 
                                        sep = "\t", header = TRUE)

counts.eQTL <- interval_cpm_counts_using %>% t()
# make sure they are in the same order
all(rownames(covariates) == rownames(counts.eQTL)) # TRUE

# write out covariate matrix
# The covariates file should be in csv or tab format, and have N rows and C columns, where N matches the number of samples in the expression file.
# WITH row and colnames for future reference of order
covariates %>% write.table(paste0(interval_path, "POPS2_comparison/PEER/PEER_covariates_colnames_update_10_27_25.csv"), quote = FALSE, sep=",")
# withOUT row and colnames (for inout to PEER)
covariates %>% write.table(paste0(interval_path, "POPS2_comparison/PEER/PEER_covariates_update_10_27_25.csv"), quote = FALSE, row.names = FALSE, sep=",",  col.names=FALSE)

# write out count matrix as csv
# The matrix is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes. 
# PEER can read comma-separated (.csv) or tab separated (.tab) files, and assumes single space separation if the extension is not .csv or .tab. 
counts.eQTL %>% write.table(file = paste0(interval_path, "POPS2_comparison/PEER/log2cpm_gene_counts_t_for_PEER_update_10_27_25.csv"), quote = FALSE, row.names = FALSE, sep=",",  col.names=FALSE)

################################################################################
#                     RUN IN BASH
################################################################################
# # run peer factors by running the following code as a script with bsub
# takes ~4.5 hours 
# # cd to the correct dir - output_dir
# cd interval_rna/POPS2_comparison/PEER
# # load the module 
# module load HGI/softpack/users/sh50/sh50_eQTL/1
# # run peer 
# peertool -f log2cpm_gene_counts_t_for_PEER_update_10_27_25.csv -c PEER_covariates_update_10_27_25.csv -n 50 --no_res_out --add_mean -i 10000 
################################################################################

# peer factors from https://github.com/PMBio/peer/wiki/Tutorial
# The output is written to directory peer_out by default, creating csv files for the residuals after accounting for the factors (residuals.csv, NxG matrix), 
#       the inferred factors (X.csv, NxK), the weights of each factor for every gene (W.csv, GxK), 
#       and the inverse variance of the weights (Alpha.csv, Kx1).
covariates <- read.csv(paste0(interval_path, "POPS2_comparison/PEER/PEER_covariates_colnames_update_10_27_25.csv"))
interval_peer_factors <- t(read.csv(paste0(interval_path, "POPS2_comparison/PEER/peer_out/X.csv"), header = FALSE))
colnames(interval_peer_factors) <- c("geno_1","geno_2","geno_3","geno_4", "Neutrophils_centered","Monocytes_centered", "1s", paste0("peer_", c(1:50)))
interval_peer_factors <- interval_peer_factors %>% as.data.frame() %>% select(contains("peer_"))
rownames(interval_peer_factors) <- rownames(covariates)
interval_peer_factors <- as.matrix(interval_peer_factors)

# write this out 
interval_peer_factors %>% write.csv(paste0(interval_path, "POPS2_comparison/PEER/Interval_peer_factors.csv"), quote = FALSE)


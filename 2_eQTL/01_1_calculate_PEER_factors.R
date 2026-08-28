# 01_1_calculate_PEER_factors.R 

################################################################################

# 1.1. Calculate PEER factors 

################################################################################

# Aim: Calculate PEER factors for use in eQTL modeling 

########################### Output paths ##########################
output_dir <- "genotyping/analysis/PEER/outputs/"

########################### Input paths ###########################
cpm_counts_inpath <- "rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt"
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
cibersort_inpath <- "rna-seq/analysis/cibersort/outputs/CIBERSORTx_Adjusted.txt"
pca_inpath <- "genotyping/analysis/eQTL/input_data/POPS2_genotyping_imputed_eQTL_pca.eigenvec"
geno_samples_inpath <- "genotyping/analysis/eQTL/input_data/genotyping_samples_filt_unrelated_w_rna_inds.txt"
geno_rna_key <- "genotyping/data/genotyping_ANON_ID_key.csv"

########################### Parameters ############################

########################### Load packages ###########################

########################### Load data ###########################
counts <- read.table(cpm_counts_inpath, header = TRUE) 
sample.info <- read.csv(sample_info_inpath) 
cibersort <- read.delim(cibersort_inpath)
plink.pca <- read.table(pca_inpath)
geno_samples_eQTL <- read.delim(geno_samples_inpath)
geno_anonid_key <-  read.csv(geno_rna_key) 

########################### Analysis ###########################

# First we prepare the data that will be used to run PEER

# what is required (from PEER package)
# 1. Expression matrix NxG (G is genes, N is n samples). Can be log normalized and should be variance-stabilized (voom)
# 2. Covariates matrix NxC (C is covariates) 
#     - Time-point
#     - Genotype PCs
#     - Imputed neutrophil proportions (scaled to mean = 0)
#     - Imputed monocyte proportions (scaled to mean = 0)

counts[1:5, 1:5]
dim(counts)
dim(sample.info)
# set PCA colnames
colnames(plink.pca) <- c("sampleid", "sampleid2", paste0("PC", 1:20))
# format geno anon ID key 
geno_anonid_key <-  geno_anonid_key %>% 
  mutate(eQTL = ifelse(Genotyping_sample_id %in% gsub(".CEL", "", geno_samples_eQTL$FID), TRUE, FALSE)) %>%
  filter(eQTL == TRUE)

# subset count matrix and sample info df to only include samples from individuals with genotyping data we will use
sample.info.eQTL <- sample.info %>% filter(ANON_ID %in% geno_anonid_key$ANON_ID)
counts.eQTL <- counts %>% select(sample.info.eQTL$RNA_sanger_sample_id) %>% t()
dim(sample.info.eQTL)
dim(counts.eQTL)

# how many genotyping PCs to include? 6 
plink.pca %>% 
  left_join(geno_anonid_key %>% left_join(sample.info %>% select(ANON_ID, an_eth_cat_str) %>% distinct(), by = "ANON_ID") %>% 
              mutate(sampleid = paste0(Genotyping_sample_id, ".CEL"),
                     an_eth_cat_str = as.factor(an_eth_cat_str))) %>%
  select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, an_eth_cat_str) %>%
  pivot_longer(-c(an_eth_cat_str, PC1), names_to = "PC", values_to = "value") %>%
  ggplot(aes(PC1, value, color = an_eth_cat_str)) + theme_bw() + geom_point() +
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 6, name = "Accent")) + facet_wrap(~PC, ncol = 4)

# make covariates combined matrix 
# need neutrophils, monocytes, genotype PCs, time points 
covariates <- sample.info.eQTL %>% 
  mutate(tp_val = 1) %>%
  select(RNA_sanger_sample_id, Sampletakenat, ANON_ID, tp_val) %>% 
  # separate sample time point 
  pivot_wider(names_from = Sampletakenat, values_from = tp_val) %>%
  mutate(`12 weeks` = ifelse(is.na(`12 weeks`), 0, `12 weeks`),
         `20 weeks` = ifelse(is.na(`20 weeks`), 0, `20 weeks`),
         `28 weeks` = ifelse(is.na(`28 weeks`), 0, `28 weeks`),
         `36 weeks` = ifelse(is.na(`36 weeks`), 0, `36 weeks`)) %>%
  left_join(cibersort %>% select(Mixture, Neutrophils, Monocytes), by = c("RNA_sanger_sample_id" = "Mixture")) %>%
  # center cell props to mean = 0
  mutate(Neutrophils_centered = scale(Neutrophils, center = TRUE, scale = FALSE),
         Monocytes_centered = scale(Monocytes, center = TRUE, scale = FALSE)) %>%
  # add genotyping PCs
  left_join(plink.pca %>% 
              mutate(Genotyping_sample_id = gsub(".CEL", "", sampleid)) %>% 
              select(Genotyping_sample_id, PC1, PC2, PC3, PC4, PC5, PC6) %>% 
              left_join(geno_anonid_key %>% select(Genotyping_sample_id, ANON_ID), 
                        by = "Genotyping_sample_id") %>% select(-Genotyping_sample_id),
            by = "ANON_ID") %>%
  select(-ANON_ID, -Neutrophils, -Monocytes, -`12 weeks`) %>%
  as.data.frame() %>%
  column_to_rownames("RNA_sanger_sample_id")


# show cell props were scaled to mean
covariates %>% ggplot(aes(x = Neutrophils_centered)) + geom_histogram()
covariates %>% ggplot(aes(x = Monocytes_centered)) + geom_histogram()

# don't include individual or batch because PEER can't handle random effects

# make sure they are in the same order
all(rownames(covariates) == rownames(counts.eQTL)) # TRUE

# write out covariate matrix
# The covariates file should be in csv or tab format, and have N rows and C columns, where N matches the number of samples in the expression file.
# WITH row and colnames
covariates %>% write.table(paste0(output_dir, "PEER_covariates_colnames.csv"), quote = FALSE, sep=",")
# withOUT row and colnames 
covariates %>% write.table(paste0(output_dir, "PEER_covariates.csv"), quote = FALSE, row.names = FALSE, sep=",",  col.names=FALSE)
# write out count matrix as csv
# The matrix is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes. 
# PEER can read comma-separated (.csv) or tab separated (.tab) files, and assumes single space separation if the extension is not .csv or .tab. 
counts.eQTL %>% write.table(file = paste0(output_dir, "logcpm_gene_counts_t_for_PEER.csv"), quote = FALSE, row.names = FALSE, sep=",",  col.names=FALSE)

################################################################################
#                     RUN PEER IN BASH
################################################################################
# # run peer factors by running the following code in bash 
# # cd to the correct dir
# cd genotyping/analysis/PEER/output
# # load the module with peertools v20120508
# module load HGI/softpack/users/sh50/sh50_eQTL/1
# # run peer 
# peertool -f logcpm_gene_counts_t_for_PEER.csv -c PEER_covariates.csv -n 50 --no_res_out --add_mean -i 10000 
################################################################################

# read in the peer factors
peer_factors <- t(read.csv(paste0(output_dir, "peer_out/X.csv"), header = FALSE))
dim(peer_factors)
peer_factors[1:5, 1:13]
colnames(peer_factors) <- c("20_weeks","28_weeks","36_weeks","Neutrophils_centered","Monocytes_centered","geno_1","geno_2","geno_3","geno_4","geno_5","geno_6", "1s", paste0("peer_", c(1:50)))

########################### Plot ###########################

# Plot these PCs (Supplementary Figure 29)
PC_plt1 <- plink.pca %>% 
  as.data.frame() %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 0.5) +
  labs(x = "PC1", y = "PC2")
PC_plt2 <- plink.pca %>% 
  as.data.frame() %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(size = 0.5) +
  labs(x = "PC3", y = "PC4")
PC_plt3 <- plink.pca %>% 
  as.data.frame() %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point(size = 0.5) +
  labs(x = "PC5", y = "PC6")
PC_plt4 <- plink.pca %>% 
  as.data.frame() %>%
  ggplot(aes(x = PC7, y = PC8)) +
  geom_point(size = 0.5) +
  labs(x = "PC7", y = "PC8")

PCA_plt <- gridExtra::grid.arrange(PC_plt1,# + labs(tag = "C"),
                                   PC_plt2,# + labs(tag = ""),
                                   PC_plt3,# + labs(tag = ""),
                                   PC_plt4,# + labs(tag = ""),
                                   ncol = 2, nrow = 2, 
                                   widths = c(1, 1), 
                                   heights = c(1, 1),
                                   layout_matrix = rbind(c(1, 2),
                                                         c(3, 4)),
                                   top = "Top 8 genotyping PCs")

# Neutrophil and monocte proportions
cibersort %>% mutate(Neutrophils = scale(Neutrophils, center = TRUE, scale = FALSE),
                     Monocytes = scale(Monocytes, center = TRUE, scale = FALSE)) %>%
  select(Neutrophils, Monocytes, Mixture) %>%
  pivot_longer(cols = -c(Mixture), names_to = "celltype", values_to = "value") %>%
  ggplot(aes(x = value)) + 
  geom_histogram(fill = "black", color = "black") + 
  facet_wrap(~celltype, scales = "free") + 
  labs(title = "Distributions of imputed cell proportions centered around zero",
       x = "Centered imputed cell proportion")
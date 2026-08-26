# 01.2_PCA.R

################################################################################

# 1.3. PCA

################################################################################

# Aim: generate PCA from count matrix

########################### Output paths ###########################
outpath <- "rna-seq/analysis/data_preprocessing/expression_PCA.rds"

########################### Input paths ###########################
voom_inpath <- "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"


########################### load packages ###########################
library(tidyverse)
library(stats)

# read in voom object
vobj <- readRDS(voom_inpath)

# transpose counts
cpm_counts_t <-t(vobj$E)

print("Doing PCA...")
Sys.time()
#PCA 
exp.pca <- cpm_counts_t %>% prcomp(center = TRUE,scale. = TRUE) # takes ~15 mins

print("PCA done")
Sys.time()
print("writing data object")
saveRDS(exp.pca, outpath)

###########################  plot ########################### 
# read PCA back in
exp.pca <- readRDS(outpath)
# read in sample info for plotting 
sample.info <- read.csv(sample_info_inpath) 

# extract PCs
pcs <- data.frame(exp.pca$x)

# Extract eigs
eigs <- exp.pca$sdev^2

# Calculate proportion of variation explained in PCs
Proportion = eigs/sum(eigs)

# Plot
pcs %>% rownames_to_column() %>%
  left_join(sample.info, by = c("rowname" = "RNA_sanger_sample_id")) %>%
  mutate(`Time-point` = gsub("_", " ", Sample_taken_at)) %>%
  select(PC1, PC2, `Time-point`) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Time-point`)) +
  geom_point(alpha = 1, size = 0.8) +
  scale_color_manual(values=c("#BCE4D8", "#83C4CB", "#439FB7", "#32769B")) +
  xlab(paste0("PC1 (", format(round(Proportion[1]*100, 1)), "%)")) +
  ylab(paste0("PC2 (", format(round(Proportion[2]*100, 1)), "%)")) +
  ggtitle("Gene expression PCA") +
  labs(color='Time-point') + 
  guides(color = guide_legend(override.aes = list(size = 3)))

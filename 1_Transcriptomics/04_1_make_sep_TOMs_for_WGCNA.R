# 04_1_make_sep_TOMs_for_WGCNA.R 

################################################################################

# 4.1. Make separate time matrices by time-point for WGCNA

################################################################################

# Aim: Prepare data for WGCNA consensus module analysis up to preparing TOM matrices
# Note: This WGCNA analysis was guided majorly by Nikhil Milind's WGCNA analysis. 
#       https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/04_Expression/scripts/wgcna_01_gene_expression.ipynb
#       Consensus module creation was guided by consensus module WGCNA tutorial. 
#       https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?e=4&preview=Consensus-NetworkConstruction-man.pdf&rlkey=3za224679kv8ulbpl5s9mwrkt&dl=0
#       Longitudinal modeling was guided by Jack Gisby's WGCNA analysis.
#       https://github.com/jackgisby/covid-longitudinal-multi-omics/blob/v1.0/notebooks/4_module_analysis.md

########################### Output paths ##########################
outfiles_path <- "rna-seq/analysis/WGCNA/outputs/sep_preprocessing/"

########################### Input paths ###########################
# sample info input 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
cpm_counts_inpath <- "rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt"

########################### Parameters ############################

########################### Load packages ###########################
print("Loading packages...")
library(WGCNA)
library(tidyverse)
library(sva)
library(spqn)
print("Packages loaded.")

########################### Load data ###########################
# normalized counts 
print("Reading in count data...")
counts <- read.table(cpm_counts_inpath, header = TRUE) 
dim(counts)
counts[1:5, 1:5]
print("Count data loaded.")

# sample info 
print("Reading in sample info...")
sample.info <- read.csv(sample_info_inpath)
dim(sample.info)
print("Sample info loaded.")

########################### Analysis ###########################

# Prep for WGCNA
print("Prep for WGCNA")
# This code comes from https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/04_Expression/scripts/wgcna_01_gene_expression.ipynb

# This function runs all steps to generate a TOM matrix on a data subset (from one time point)
get_TOM <- function(counts_sub, tp){
  # print the dimention of the count matrix subset
  print(dim(counts_sub))
  # print the time point 
  print(tp)
  print("get sva")
  # Run sva network
  # From Nikhil:  Gene expression data may contain counfounding technical 
  #               variation that can lead to spurious associations between genes. 
  #               This is based on an analysis conducted by Parsana et al. 
  #               This can be corrected using their method implemented in the sva package.
  # I use the same parameter (20) as Nikhil did.
  gene.exp <- sva::sva_network(counts_sub, 20)
  # transpose the matrix
  gene.exp <- t(gene.exp)
  print("check ok")
  # Run WGCNA checks
  # From Nikhil:  The following function identifies genes with excessive missing 
  #               values and outlier samples. It is generally more important for 
  #               microarray data.
  gsg <- WGCNA::goodSamplesGenes(gene.exp)
  # Output whether all genes are OK
  gsg$allOK
  # Save this pre-processed RNA-seq dataset as a file
  saveRDS(gene.exp, paste0(outfiles_path, "gene.exp_tp", tp, ".rds"))
  print("run bicor")
  # Generate correlation matrix from (pre-processed) gene count matrix 
  cor.mtx <- bicor(gene.exp, use="pairwise.complete.obs", pearsonFallback="individual")
  print("write out")
  # Save this correlation matrix 
  saveRDS(cor.mtx, paste0(outfiles_path, "bicor.mtx_tp", tp, ".rds"))
  print("done")
  # Read back in the correlation matrix 
  cor.mtx <- readRDS(paste0(outfiles_path, "bicor.mtx_tp", tp, ".rds"))
  # Get average expression of each gene across the raw count matrix 
  avg.exp <- apply(counts_sub, 1, mean)
  # Correct mean-correlation bias 
  # From Nikhil:  There is a mean-correlation bias when building co-expression 
  #               modules, as discussed by Wang et al. 
  #               This can be corrected using a method they developed.
  # Plot original data
  plot_signal_condition_exp(cor.mtx, avg.exp, signal=0.001)
  # Normalize correlation matrix (I use parameters based on optimizing 
  # similarity of dist across all samples)
  cor.mtx.norm <- normalize_correlation(cor.mtx, ave_exp=avg.exp, ngrp=10, size_grp=2000, ref_grp=10)
  rownames(cor.mtx.norm) <- rownames(cor.mtx)
  colnames(cor.mtx.norm) <- colnames(cor.mtx)
  # Plot updated distribution 
  plot_signal_condition_exp(cor.mtx.norm, avg.exp, signal=0.001)
  # Save the corrected correlation matrix 
  saveRDS(cor.mtx.norm, file=paste0(outfiles_path, "cor.mtx.norm_tp", tp, ".rds"))
  # Soft Threshold
  # From Nikhil:  The adjacency matrix generated using the correlation function 
  #               does not fit a scale-free topology. A scaling function is used 
  #               to generate the appropriate gene degree distribution. A 
  #               parameter of this scaling function, is estimated by trying 
  #               multiple values and checking how well the resulting network 
  #               fits the assumptions of a scale-free network. 
  #               Here, I try all integer values in 1-20.
  # Set powers to test
  powers = seq(1, 20, by=1)
  # Back up correlation matrix 
  cor.mtx.norm_for_signedhybrid = cor.mtx.norm
  # Update correlation matrix so any negative correlation values are 0 
  # to suit signed hybrid module type
  cor.mtx.norm_for_signedhybrid[cor.mtx.norm_for_signedhybrid < 0] = 0
  # Run pick soft thresholds
  soft.thresholds = suppressWarnings(
    pickSoftThreshold.fromSimilarity(cor.mtx.norm_for_signedhybrid, powerVector=powers, verbose=4)
  )
  # Get power estimate
  soft.threshold = soft.thresholds$powerEstimate
  # Plot results  
  options(repr.plot.width=12, repr.plot.height=6)
  p1 <- soft.thresholds$fitIndices %>%
    ggplot() +
    geom_text(aes(x=Power, y=SFT.R.sq, label=Power)) +
    xlab("Soft Threshold Power") + ylab(bquote(R^2*" of Fit to Scale-Free Network")) + geom_hline(yintercept = .85)
  p2 <- soft.thresholds$fitIndices %>%
    ggplot() +
    geom_text(aes(x=Power, y=mean.k., label=Power)) +
    xlab("Soft Threshold Power") + ylab("Mean Connectivity") + geom_hline(yintercept = 100)
  # Write out the plots to evaluate later (Supplementary Figure 24)
  soft_thresh <- gridExtra::grid.arrange(p1 + labs(tag = "A"), p2 + labs(tag = "B"), nrow=1)
  ggsave(paste0(outfiles_path, "soft.thresh.plt_tp", tp, ".pdf"), soft_thresh, width = 9, height = 4)
  # Print soft threshold value
  print(soft.threshold)
  # Across all time points, the threshold of 13 was chosen, so we will use that value 
  # Make adjacency matrix signed hybrid 
  # From Nikhil:  The adjacency matrix defines the similarity between all genes 
  #               in the sample data. I am using the biweight midcorrelation 
  #               (bicor) function to estimate similarity. It is a median-based 
  #               approach that is less sensitive to outliers.
  adjacency.matrix = adjacency.fromSimilarity(cor.mtx.norm, type="signed hybrid", power=13)
  # Write this out 
  saveRDS(adjacency.matrix, file=paste0(outfiles_path, "adjacency.mtx_st13_tp", tp, ".rds"))
  # Make tom matrix (type signed to match signed hybrid)
  # From Nikhil:  the topological overlap metric (TOM) matrix is generated from 
  #               the adjacency matrix. For details, check the original WGCNA paper.
  TOM.matrix = TOMsimilarity(adjacency.matrix, TOMType = "signed", verbose = 4)
  # Write this out 
  saveRDS(TOM.matrix, file=paste0(outfiles_path, "tom.mtx_st13_tp", tp, ".rds"))
}

# Run the script on all time-point subsets
get_TOM(counts_sub = (counts %>% select((sample.info %>% filter(Sample_taken_at == "12_weeks"))$RNA_sanger_sample_id)), tp = 1)
get_TOM(counts_sub = (counts %>% select((sample.info %>% filter(Sample_taken_at == "20_weeks"))$RNA_sanger_sample_id)), tp = 2)
get_TOM(counts_sub = (counts %>% select((sample.info %>% filter(Sample_taken_at == "28_weeks"))$RNA_sanger_sample_id)), tp = 3)
get_TOM(counts_sub = (counts %>% select((sample.info %>% filter(Sample_taken_at == "36_weeks"))$RNA_sanger_sample_id)), tp = 4)
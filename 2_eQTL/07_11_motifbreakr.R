# 07_11_motifbreakr.R

################################################################################

# 7.11. Run motifbreakR on eQTL interaction eSNPs

################################################################################

# Aim: Identify eSNPs lying in transcription factor binding motifs 

########################### Paths ##########################
path <- "genotyping/"
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Input paths ###########################

########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)

########################### Load data ###########################
# all time interaction eQTL effects for enrichment
time_interactions <- read.csv(paste0(path,"analysis/eQTL/output_data/interactions_flanders_2/time_interactions_flanders_2.csv"))
# id name conversion
gtf <- rtracklayer::import(gtf_inpath)
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)
# id name conversion
gtf <- rtracklayer::import(gtf_inpath)
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)
# bed file for motifbreakr
BED_file <-  readRDS(paste0(path, "/analysis/eQTL/input_data/interval_comp/full_pops_BED_file_for_motifbreakR.rds"))
# time interaction clusters
time_int_clust <- read.csv(paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/time_interaction_clusters_2.csv")) %>%
  mutate(cluster_label = cluster_label + 1) %>%
  separate(gsp, sep = "_", into = c("Gene", "SNP"), remove = FALSE, extra = "merge")

########################### Analysis ###########################

# Run motifbreakR on all cluster t4 eSNPs
# Generate list of cluster t4 eSNPs
cluster_t4_eSNPs <- time_interactions %>% mutate(gsp = paste0(Gene, "_", SNP)) %>%
  filter(gsp %in% (time_int_clust %>% filter(cluster_label == 4) %>% pull(gsp))) %>%
  pull(SNP)
length(cluster_t4_eSNPs) # [1] 309

# Make the BED file format but keep SNP IDs so I can filter this file based on the SNPs we eventually want to include in motifbreakR.
BED_file_formatted <- BED_file %>%
  # remove SNPs that weren't in the dbSNP ref so don't have a ref and alt
  dplyr::filter(!is.na(REF),
                # filter to only SNPs mashed and high LD snps for motifbreakR
                snpID %in% cluster_t4_eSNPs) %>%
  dplyr::mutate(
    chr = ifelse(chr == 23, "X", chr),
    chrom   = paste0("chr", chr),
    chromStart = as.integer(pos - 1),  # convert 1-based to 0-based
    chromEnd   = as.integer(pos),       # end is exclusive in BED
    ID_noRS = paste0(chrom, ":", pos, ":", REF, ":", ALT),# could use this if you want rsIDs: ifelse(grepl("rs", dbSNP_rsID), dbSNP_rsID, paste0(chrom, ":", pos, ":", REF, ":", ALT)),  
    score = 0,
    strand = "+"
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, ID_noRS, score, strand) %>%
  distinct()
head(BED_file_formatted)
dim(BED_file_formatted) # 296   6 # the 13 missing SNPs from the cluster t4 eSNP list weren't in the dbSNP file so didn't have REF and ALT alleles 

# write this out 
BED_file_formatted %>% write.table(paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/cluster_t4_signal_BED_file_for_motifbreakR_hg38.bed"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Run motifbreakR
# read in the snp file
snp.bed.file.path = paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/cluster_t4_signal_BED_file_for_motifbreakR_hg38.bed")
print(paste("bed file is ", snp.bed.file.path))
# read in bed file (just to check it is correct)
BED_file_formatted = read.table(snp.bed.file.path, header = FALSE)
print(paste("N SNPs being tested: ", dim(BED_file_formatted)[1]))
print("Head:")
print(head(BED_file_formatted))
print("Tail:")
print(tail(BED_file_formatted))
#import the BED file
print("Importing bed file")
Sys.time()
snps.mb.frombed <- snps.from.file(file = snp.bed.file.path,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")
print("bed file imported")
Sys.time()

# run motifbreakR 
data(motifbreakR_motif)
print(paste("n cores to run in parallel: ", BiocParallel::bpworkers(BiocParallel::bpparam())))
print("run motifbreakR")
Sys.time()
motifbreakr.results <- motifbreakR(snpList = snps.mb.frombed, # SNP list
                                   pwmList = motifbreakR_motif, # motifs you want to test - in our case homo sapiens 
                                   threshold = 0.85, # the default, indicating the sequence matches the PMW very strongly (out of 1) - so it is only reporting a motif hit if either the ref or alt allele make it a strong binding site. Best for avoiding false positives... 
                                   verbose = TRUE, 
                                   show.neutral = TRUE) # don't filter out motifs that don't change affinity bewteen alleles

print("finished running motifbreakR")
Sys.time() # Takes approx. 40 mins

# write this out as an rdata object
path_out <- paste0(paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/motifbreakR_results_cluster_t4_eSNPs.rds"))
print(paste("saving output as rds to: ", path_out))
saveRDS(motifbreakr.results, path_out)
print("finished saving")

# Read in motifbreakR results
motifbreakR_results <- readRDS(paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/motifbreakR_results_cluster_t4_eSNPs.rds")) 

motifbreakR_results <- as.data.frame(mcols(motifbreakr.results)) %>%
  tidyr::separate(SNP_id, 
           into = c("chr", "pos", "ref", "alt"),
           sep = ":", 
           remove = FALSE) 

# subset motifbreakR results to SMAD3 and STAT1 motifs in eQTL with eGenes downstream of the respective TFs 
# read in info about eQTL with eGenes downstream of TFs
SMAD3_examples <- read.csv(paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/time_int_cluster_4_SMAD3_downstream_eQTL.csv"))
STAT1_examples <- read.csv(paste0(path, "analysis/eQTL/output_data/interactions_flanders_2/time_int_cluster_4_STAT1_downstream_eQTL.csv"))

# Makes sup table 10
SMAD3_examples %>% mutate(chr = paste0("chr", chr),
                                            SNPpos = as.character(SNPpos)) %>% 
  left_join(motifbreakR_results, by = c("SNPpos" = "pos", "chr")) %>%
  filter(grepl("SMAD3", geneSymbol)) %>% 
  select(Gene, SNP, SNP_id, geneSymbol, dataSource, providerId, alleleEffectSize)

STAT1_examples %>% mutate(chr = paste0("chr", chr),
                          SNPpos = as.character(SNPpos)) %>% 
  left_join(motifbreakR_results, by = c("SNPpos" = "pos", "chr")) %>%
  filter(grepl("STAT1", geneSymbol)) %>% 
  select(Gene, SNP, SNP_id, geneSymbol, dataSource, providerId, alleleEffectSize)



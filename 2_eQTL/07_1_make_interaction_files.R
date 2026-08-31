# 07_1_make_interaction_files.R

################################################################################

# 7.1.  Make files for eQTL interaction analysis

################################################################################

# Aim: Format eQTL interaction input files 

########################### Paths ##########################
path <- "genotyping/analysis/eQTL/"

########################### Input paths ###########################

########################### Parameters ############################

########################### Load packages ###########################

########################### Load data ###########################

########################### Analysis ###########################

# Take conditionally independent eQTL signals to test for interactions 
res <-  readRDS(paste0(path, "output_data/flanders_conditional_eQTL/all_conditional_signals_pops.rds"))

snps.sig <- unique(res$snps)

# Subset the genotyping data to just those SNPs
geno.int <- do.call(cbind, lapply(c(1:22, "X"), function(i){
  print(i)
  load(paste0(path, "input_data/eqtl_files_", i, ".rda"))
  geno[, which(colnames(geno) %in% snps.sig)]
}))

# Load the rest of the required data for interaction mapping
load(paste0(path, "input_data/eqtl_files_22.rda"))

# Make a list of SNP-gene pairs to test
pairs.int <- res[, c("phenotype_id", "snps")]
pairs.int$Gene <- as.character(pairs.int$phenotype_id)
pairs.int$SNP <- as.character(pairs.int$snps)
pairs.int <- pairs.int[, c("Gene", "SNP")]
dim(pairs.int) # 23859     2 #  31734     2
rownames(pairs.int) <- 1:nrow(pairs.int)

# read in the outcome information 
load(paste0(path, "input_data/eqtl_outcome_data.rda"))
# ensure the outcome data is in the proper format 
pe_char <- pe
pe <- as.logical(pe)
gdm_char <- gdm
gdm <- as.logical(gdm)
sga_char <- sga
sga <- as.logical(sga)
ptb_char <- ptb
ptb <- as.logical(ptb)

# Make the R data file containing everything needed for interaction testing
int_file <- paste0(path, "input_data/eqtl_int_files_flanders.rda")
save(list=c("gene_exp", "geno.int", "covariates", "individual", "batch", "peer_factors", "timepoint", "gdm", "pe", "sga", "ptb", "pairs.int"),
     file = int_file)

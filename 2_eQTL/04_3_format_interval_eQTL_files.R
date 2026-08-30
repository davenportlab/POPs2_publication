# 04_3_format_interval_eQTL_files.R

################################################################################

# 4.3. Format interval files for eQTL mapping 

################################################################################

# Aim: Prepare interval files for eQTL mapping 

########################### Paths ##########################
lustre_path <- "interval_rna/" 
output_dir <- paste0(lustre_path, "POPS2_comparison/eQTL_input_data/")
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Load packages ###########################
library(tidyverse)

########################### Analysis ###########################

# 1. Read in normalized count matrix
counts.eQTL <- read.table(paste0(lustre_path, "POPS2_comparison/data_subset/Interval_counts_female_less45_log2cpm.txt"), header = TRUE)
counts.eQTL[1:5, 1:5]
dim(counts.eQTL)

# read in covariates from PEER analysis
peer_covariates <- read.csv(paste0(lustre_path, "POPS2_comparison/PEER/PEER_covariates_colnames_update_10_27_25.csv"))

# read in interval metadata
metadata <- read.csv(paste0(lustre_path, "POPS2_comparison/data_subset/Interval_metadata_female_less45.csv"))

covariates <- metadata %>% 
  select(affymetrix_ID, sample_id, sequencingBatch) %>%
  left_join(peer_covariates %>% rownames_to_column(var = "sample_id"))

# write out final covariates file 
write.csv(covariates, file = paste0(output_dir, "eQTL_covariates_update_10_27_25.csv"), row.names = FALSE, quote = FALSE)

# check count matrix and covariates have same samples in same order
all(colnames(counts.eQTL) == covariates$sample_id)
# check all genotyping samples match 
fam <- read.table(paste0(lustre_path, "POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt.fam"))
setequal(fam$V1, covariates$affymetrix_ID)

################################################################################

# Make a list of SNP-gene pairs to test (SNPs 1Mb either side of TSS for each gene)

################################################################################
# prepare gene position information file
gtf <- as.data.frame(rtracklayer::import(gtf_inpath))
gene.info.filt <- gtf[match(rownames(counts.eQTL), gtf$gene_id), ]
rownames(gene.info.filt) <- rownames(counts.eQTL)
all(gene.info.filt$gene_id == rownames(counts.eQTL))
gene.info.filt <- gene.info.filt[, c("seqnames", "start", "end", "strand", 
                                     "gene_id", "gene_name")]
gene.info.filt$seqnames <- as.character(gene.info.filt$seqnames)
gene.info.filt$seqnames[gene.info.filt$seqnames == "X"] <- 23
write.table(gene.info.filt, paste0(output_dir, "gene_info_update_10_27_25.txt"), sep = "\t", quote = FALSE)

# edit the non-rsID SNP IDs for the X chromosome to be in the same format as all of the others (alphabetical)
bim <- data.frame(data.table::fread(paste0(lustre_path, "POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt.bim"), 
                                    header=F, stringsAsFactors = F))
bim_update_x <- bim %>%
  mutate(V2 = ifelse(
    # X chr not rsID both alleles length 1
    V1 == 23 & !grepl("^rs", V2) & nchar(V5) == 1 & nchar(V6) == 1, 
    paste0("X_", V4, "_", pmin(V5, V6), "_", pmax(V5, V6)),
    V2), 
    # replace : with _ for the indels too
    V2 = gsub(":", "_", V2))

bim_update_x %>% write.table(paste0(lustre_path, "POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt.bim"),
                             quote = FALSE,
                             sep = "\t",
                             row.names = FALSE,
                             col.names = FALSE)

# filter the plink file so it is only biallelic SNPs and MAF 0.05
################################################################################
#                     RUN IN BASH
################################################################################
# # start interactive job
# fash 20
# # load module with plink 
# module load HGI/softpack/groups/trynka/popgen/1.0
# cd interval_rna/POPS2_comparison/data_subset/genotyping
# plink \
# --bfile INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt \
# --biallelic-only strict --snps-only just-acgt --make-bed \
# --out INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_snps
# # filter to MAF 0.05 for the eQTL analysis 
# plink \
# --bfile INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_snps \
# --geno 0.05 --hwe 0.00001 --maf 0.05 --make-bed \
# --out INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_maf0.05_snps
################################################################################

# read in SNP and gene positional information
snp.info <- data.frame(data.table::fread(paste0(lustre_path, "POPS2_comparison/data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_maf0.05_snps.bim"), 
                                         header=F, stringsAsFactors = F))

genes <- read.delim(paste0(output_dir, "gene_info_update_10_27_25.txt"), stringsAsFactors = F)

# flip start and end for genes on - strand
is_neg <- genes$strand == "-"
tmp <- genes$start[is_neg]
genes$start[is_neg] <- genes$end[is_neg]
genes$end[is_neg] <- tmp

# divide up genes and SNPs by chromosome
genes.list <- vector("list", 23)
chr.bim.list <- vector("list", 23)

for(i in 1:23){
  genes.list[[i]] <- subset(genes, seqnames == i)
  rownames(genes.list[[i]]) <- genes.list[[i]][, "gene_id"]
  
  chr.bim.list[[i]] <- snp.info[which(snp.info$V1 == i), ]
  chr.bim.list[[i]] <- chr.bim.list[[i]][, c(1, 2, 4)]
  rownames(chr.bim.list[[i]]) <- chr.bim.list[[i]][, 2]
}

# make a list of SNP-gene pairs to test
results <- vector("list", 23)
# for each chromosome
for(c in 1:23){
  results.c <- vector("list", nrow(genes.list[[c]]))
  # for each gene get SNPs in a 1Mb window
  for (i in 1:nrow(genes.list[[c]])) {
    snp.names <- rownames(chr.bim.list[[c]])[which(as.numeric(genes.list[[c]][i, 1]) == chr.bim.list[[c]][, 1] # match chromosome
                                                   & (chr.bim.list[[c]][, 3] > (genes.list[[c]][i, 2] - 1000000) # position greater than 1mb less than gene start
                                                      &  chr.bim.list[[c]][, 3] < (genes.list[[c]][i, 2] + 1000000)))] # position less than 1mp greater than gene start
    if(length(snp.names) > 1){
      results.c[[i]] <- data.frame(rownames(genes.list[[c]])[i], snp.names, 
                                   rep(as.numeric(genes.list[[c]][i, 1]), length(snp.names)))
    }
  }
  results[[c]] <- data.table::rbindlist(results.c)
  colnames(results[[c]]) <- c("Gene", "SNP", "Chr")
}

results.2 <- data.table::rbindlist(results)
results.2[1:5, 1:3] # dim 63867702        3
dim(results.2) # 20077 distinct genes 
write.table(results.2, paste0(output_dir, "gene_snp_pairs_cis_rnaseq_update_10_27_25.txt"), sep="\t", 
            quote=FALSE, row.names=FALSE)
# list of snps 
snps <- unique(as.character(results.2$SNP)) # length = 4279133
write.table(snps, paste0(output_dir, "snps_for_cis_eqtl_rnaseq_update_10_27_25.txt"), sep="\t", quote=F, 
            row.names = F, col.names = F)

# # make genotyping files with only the snps we will use 
################################################################################
#                     RUN IN BASH
################################################################################
# # start interactive job 
# fash 20
# # load module
# module load HGI/softpack/users/sh50/sh50_geno_qc/1
# # cd to correct place
# cd interval_rna/POPS2_comparison/data_subset/genotyping

# #separate them into different files 
# mkdir Imputed_Genotyping_Data_byChr
# for P in {1..22} ; do plink --bfile INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_maf0.05_snps --chr $P --make-bed --out Imputed_Genotyping_Data_byChr/Interval_f_u45_chr.${P}.imputed.eqtl; done
# plink --bfile INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_maf0.05_snps --chr X --make-bed --out Imputed_Genotyping_Data_byChr/Interval_f_u45_chr.X.imputed.eqtl

# cd ../../eQTL_input_data
# mkdir geno_files_update_10_27_25

# # extract the snps we will use and recode the files as raw format
# extract only biallelic snps for eQTL mapping 
# for CHR in {1..22} ; do plink --bfile /path/to/interval_rna/POPS2_comparison/data_subset/genotyping/Imputed_Genotyping_Data_byChr/Interval_f_u45_chr.${CHR}.imputed.eqtl --extract snps_for_cis_eqtl_rnaseq_update_10_27_25.txt --biallelic-only strict --snps-only just-acgt --recode A --out geno_files_update_10_27_25/Interval_f_u45_genotyping_for_rna-seq_eqtl_${CHR} --allow-extra-chr --chr ${CHR} ; done
# plink --bfile /path/to/interval_rna/POPS2_comparison/data_subset/genotyping/Imputed_Genotyping_Data_byChr/Interval_f_u45_chr.X.imputed.eqtl --extract snps_for_cis_eqtl_rnaseq_update_10_27_25.txt --biallelic-only strict --snps-only just-acgt --recode A --out geno_files_update_10_27_25/Interval_f_u45_genotyping_for_rna-seq_eqtl_X --allow-extra-chr --chr X
################################################################################

# assemble all files

# read in gene expression, covariates, sample key, peer factors and make sure they are all in the same order
covariates <- read.csv(paste0(output_dir, "eQTL_covariates_update_10_27_25.csv")) %>% 
  column_to_rownames(var = "sample_id") #%>%
#select(-affymetrix_ID)
geno_id <- covariates$affymetrix_ID
gene_exp <- read.table(paste0(lustre_path, "POPS2_comparison/data_subset/Interval_counts_female_less45_log2cpm.txt"), header = TRUE)
all(colnames(gene_exp) == rownames(covariates))
batch <- paste0("batch_", covariates$sequencingBatch)
peer_factors <- read.csv(paste0(lustre_path, "POPS2_comparison/PEER/Interval_peer_factors.csv"), header = TRUE, row.names = 1)
# check the order is OK 
all(rownames(covariates) == rownames(peer_covariates))
peer_factors <- as.matrix(peer_factors)
# remove genotyping ID and batch from covariates
covariates <- covariates %>% select(-affymetrix_ID, -sequencingBatch) %>% as.matrix()

# set up genotyping files
for(i in 1:22){
  print(i)
  geno <- data.frame(data.table::fread(paste(output_dir, "geno_files_update_10_27_25/Interval_f_u45_genotyping_for_rna-seq_eqtl_", i, ".raw", sep=""),
                                       sep=" ", drop = 2:6))
  m1 <- match(geno_id, geno$FID)
  # make one row for each RNA-seq sample based on matched genotyping sample - in the same order
  geno <- geno[m1, ]
  # add the paired RNA sample id as rowname
  rownames(geno) <- rownames(covariates)
  # remove X from the non-rsid snps 
  colnames(geno) <- gsub("X", "", colnames(geno))
  # remove allele
  colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
  # remove genotyping sample ID
  geno[, 1] <- NULL
  geno <- as.matrix(geno)
  print("all snps in pairs same as those in genotyping matrix:")
  print(setequal(colnames(geno), results.2 %>% filter(Chr == i) %>% pull(SNP)))
  eQTL_files <- paste(output_dir, "eqtl_files_", i, ".rda", sep="")
  save(list=c("gene_exp", "geno", "covariates", "batch", "peer_factors"),
       file = eQTL_files)
}

# needs to be slightly adjusted for X chromosome for the non-rsID SNP IDs
i = "X"
print(i)
geno <- data.frame(data.table::fread(paste(output_dir, "geno_files_update_10_27_25/Interval_f_u45_genotyping_for_rna-seq_eqtl_", i, ".raw", sep=""),
                                     sep=" ", drop = 2:6))
m1 <- match(geno_id, geno$FID)
# make one row for each RNA-seq sample based on matched genotyping sample - in the same order
geno <- geno[m1, ]
# add the paired RNA sample id as rowname
rownames(geno) <- rownames(covariates)
# remove allele
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
# remove genotyping sample ID
geno[, 1] <- NULL
geno <- as.matrix(geno)
print("all snps in pairs same as those in genotyping matrix:")
print(setequal(colnames(geno), results.2 %>% filter(Chr == 23) %>% pull(SNP)))
eQTL_files <- paste(output_dir, "eqtl_files_", i, ".rda", sep="")
save(list=c("gene_exp", "geno", "covariates", "batch", "peer_factors"),
     file = eQTL_files)
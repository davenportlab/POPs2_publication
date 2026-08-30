# 05_1_prep_single_tp_eQTL_files.R

################################################################################

# 5.1. Prepare input files for eQTL mapping separately at each time-point 

################################################################################

# Aim: subset metadata to the individuals who have only one sample per time-point 

########################### Paths ###########################
lustre_path <- "genotyping/analysis/"
dir <- paste0(lustre_path, "eQTL/")
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Load packages ###########################
library(dplyr)

########################### Analysis ###########################

# read in gene expression, covariates, sample key, peer factors and make sure they are all in the same order
covariates <- read.csv(paste0(dir, "input_data/eQTL_covariates.csv"))

# make a list of individulas with all four time-points
covariates_sub <- covariates %>%
  group_by(ANON_ID) %>%
  filter(all(c("12_weeks", "20_weeks", "28_weeks", "36_weeks") %in% Sample_taken_at)) %>%
  ungroup()

# make a list of genotyping samples to filter to 
geno_sample_ids_fourtp <- covariates_sub %>% pull(Genotyping_sample_id) %>% unique()
cbind(paste0(geno_sample_ids_fourtp, ".CEL"), paste0(geno_sample_ids_fourtp, ".CEL")) %>%
  write.table(paste0(dir, "input_data/single_tp_eQTL/genotyping_sample_ids_all_four_RNAseq.txt"), sep="\t", quote=F, 
              row.names = F, col.names = F)

# Subset genotyping data to individuals with RNA-seq samples at all four time-points
# And re-filter to MAF < 0.05
################################################################################
#           Run in bash 
################################################################################
# # start interactive job 
# fash 20
# # load module
# module load HGI/softpack/users/sh50/sh50_geno_qc/1
# # cd to correct place
# cd genotyping/analysis/eQTL/input_data/single_tp_eQTL
# filter to the individuals of interest 
# plink \
# --bfile genotyping/data/imputed_genotyping_dataset_for_eQTL/POPS2_genotyping_imputed_eQTL_rsid \
# --keep genotyping_sample_ids_all_four_RNAseq.txt \
# --make-bed \
# --out POPS2_genotyping_imputed_eQTL_rsid_fourRNA
# #now filter variants including MAF 0.05
# plink \
# --bfile POPS2_genotyping_imputed_eQTL_rsid_fourRNA \
# --geno 0.05 --hwe 0.00001 --maf 0.05 --make-bed \
# --out POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps
################################################################################

# 1. Read in eQTL count matrix
counts.eQTL <- readRDS(paste0(dir, "input_data/eQTL_counts.rds"))
dim(counts.eQTL)
counts.eQTL[1:5, 1:5]

# subset count matrix and sample info df to only include samples from the individuals we will include
counts.eQTL <- counts.eQTL %>% select(covariates_sub$RNA_sanger_sample_id)
dim(counts.eQTL)

# write this out 
saveRDS(counts.eQTL, file = paste0(dir, "input_data/single_tp_eQTL/eQTL_counts_fourRNA.rds"))

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
write.table(gene.info.filt, paste0(dir, "input_data/single_tp_eQTL/gene_info_fourRNA.txt"), sep = "\t")

# read in SNP and gene positional information
snp.info <- data.frame(data.table::fread(paste0(dir, "input_data/single_tp_eQTL/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps.bim"), 
                                         header=F, stringsAsFactors = F))
genes <- read.delim(paste0(dir, "input_data/single_tp_eQTL/gene_info_fourRNA.txt"), stringsAsFactors = F)

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
results.2[1:5, 1:3] # dim 58259322        3
write.table(results.2, paste0(dir, "input_data/single_tp_eQTL/gene_snp_pairs_cis_rnaseq_fourRNA.txt"), sep="\t", 
            quote=FALSE, row.names=FALSE)
# list of snps 
snps <- unique(as.character(results.2$SNP))
write.table(snps, paste0(dir, "input_data/single_tp_eQTL/snps_for_cis_eqtl_rnaseq_fourRNA.txt"), sep="\t", quote=F, 
            row.names = F, col.names = F)

# # format genotyping files with only the snps we will use 
################################################################################
#           Run in bash 
################################################################################

# # start interactive job 
# fash 20
# # load module
# module load HGI/softpack/users/sh50/sh50_geno_qc/1
# # cd to correct place
# cd genotyping/analysis/eQTL/input_data/single_tp_eQTL

# #separate them into different files 
# mkdir Imputed_Genotyping_Data_byChr
# for P in {1..22} ; do plink --bfile POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps --chr $P --make-bed --out Imputed_Genotyping_Data_byChr/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_chr.${P}.imputed.eqtl; done
# plink --bfile POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps --chr X --make-bed --out Imputed_Genotyping_Data_byChr/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_chr.X.imputed.eqtl

# mkdir geno_files

# # extract the snps we will use and recode the files as raw format
# extract only biallelic snps for eQTL mapping 
# for CHR in {1..22} ; do plink --bfile Imputed_Genotyping_Data_byChr/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_chr.${CHR}.imputed.eqtl --extract snps_for_cis_eqtl_rnaseq_fourRNA.txt --recode A --out geno_files/genotyping_for_rna-seq_eqtl_fourRNA_${CHR} --allow-extra-chr --chr ${CHR} ; done
# plink --bfile Imputed_Genotyping_Data_byChr/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_chr.X.imputed.eqtl --extract snps_for_cis_eqtl_rnaseq_fourRNA.txt --recode A --out geno_files/genotyping_for_rna-seq_eqtl_fourRNA_X --allow-extra-chr --chr X

################################################################################

# assemble all files
# read in gene expression, covariates, sample key, peer factors and make sure they are all in the same order
covariates <- read.csv(paste0(dir, "input_data/eQTL_covariates.csv"))

# make a list of individulas with all four time-points
covariates_sub <- covariates %>%
  group_by(ANON_ID) %>%
  filter(all(c("12_weeks", "20_weeks", "28_weeks", "36_weeks") %in% Sample_taken_at)) %>%
  ungroup() %>%
  column_to_rownames(var = "RNA_sanger_sample_id")

# 1. Read in eQTL count matrix
gene_exp <- readRDS(paste0(dir, "input_data/single_tp_eQTL/eQTL_counts_fourRNA.rds"))
dim(gene_exp)
all(colnames(gene_exp) == rownames(covariates_sub))
geno_id_CEL <- paste0(covariates_sub$Genotyping_sample_id, ".CEL")
individual <- covariates_sub$ANON_ID
batch <- covariates_sub$id_run_position
timepoint <- as.factor(covariates_sub$Sample_taken_at)
peer_factors <- t(read.csv(paste0(lustre_path, "PEER/outputs/peer_out/X.csv"), header = FALSE))
colnames(peer_factors) <- c("20_weeks","28_weeks","36_weeks","Neutrophils_centered","Monocytes_centered","geno_1","geno_2","geno_3","geno_4","geno_5","geno_6", "1s", paste0("peer_", c(1:50)))
peer_factors <- peer_factors %>% as.data.frame() %>% select(contains("peer_"))
rownames(peer_factors) <- covariates$RNA_sanger_sample_id
peer_factors_sub <- peer_factors[rownames(peer_factors) %in% rownames(covariates_sub),] 
peer_factors_sub <- as.matrix(peer_factors_sub)

covariates_sub <- covariates_sub %>% select(-Genotyping_sample_id, -ANON_ID, -id_run_position, -Sample_taken_at) %>% as.matrix()

# set up genotyping files
for(i in c(1:22, "X")){
  print(i)
  geno <- data.frame(data.table::fread(paste(dir, "input_data/single_tp_eQTL/geno_files/genotyping_for_rna-seq_eqtl_fourRNA_", i, ".raw", sep=""),
                                       sep=" ", drop = 2:6))
  m1 <- match(geno_id_CEL, geno$FID)
  # make one row for each RNA-seq sample based on matched genotyping sample - in the same order
  geno <- geno[m1, ]
  # add the paired RNA sample id as rowname
  rownames(geno) <- rownames(covariates_sub)
  # remove X from the non-rsid snps 
  colnames(geno) <- gsub("X", "", colnames(geno))
  # remove allele
  colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
  # remove genotyping sample ID
  geno[, 1] <- NULL
  geno <- as.matrix(geno)
  
  eQTL_files <- paste(dir, "input_data/single_tp_eQTL/eqtl_files_", i, ".rda", sep="")
  save(list=c("gene_exp", "geno", "covariates_sub", "individual", "batch", "
              ", "timepoint"),
       file = eQTL_files)
}
# 01_2_prepare_files_for_ciseqtl_mapping.R 

################################################################################

# 1.2. Prepare files for cis-eQTL mapping

################################################################################

# Aim: Format data for input into eQTL models

########################### Output paths ##########################
covariate_outpath <- "genotyping/analysis/eQTL/input_data/eQTL_covariates.csv"
eQTL_counts_outpath <- "genotyping/analysis/eQTL/input_data/eQTL_counts.rds"
output_dir <- "genotyping/analysis/eQTL/"

########################### Input paths ###########################
cpm_counts_inpath <- "rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt"
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
geno_rna_inpath <- "rna-seq/data/07_03_24_unswapped_clinical_rna_genotyping_ind_key_final3.csv"
geno_samples_inpath <- "genotyping/analysis/eQTL/input_data/genotyping_samples_filt_unrelated_w_rna_inds.txt"
peer_covariate_inpath <- "genotyping/analysis/PEER/outputs/PEER_covariates_colnames.csv"
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"
bim_inpath <- "genotyping/data/imputed_genotyping_dataset_for_eQTL/POPS2_genotyping_imputed_eQTL_rsid.bim"
peer_inpath <- "genotyping/analysis/PEER/outputs/peer_out/X.csv"
geno_files_inpath <- "genotyping/analysis/eQTL/input_data/geno_files/"

########################### Parameters ############################

########################### Load packages ###########################

########################### Load data ###########################
counts <- read.table(cpm_counts_inpath, header = TRUE) 
sample.info <- read.csv(sample_info_inpath) 
geno_rna_key <- read.csv(geno_rna_inpath) 
geno_samples_eQTL <- read.delim(geno_samples_inpath)
peer_covariates <- read.csv(peer_covariate_inpath)
gtf <- rtracklayer::import(gtf_inpath)

########################### Analysis ###########################

# Prepare sample ID key 
geno_rna_key <- geno_rna_key %>% # genotyping RNA-seq key 
  select(RNA_sanger_sample_id, ANON_ID, Sample_taken_at, Genotyping_sample_id) %>%
  left_join(sample.info %>% # clin tech var 
              select(RNA_sanger_sample_id, id_run_position), by = "RNA_sanger_sample_id")
dim(geno_rna_key) # [1] 3776    5
# subset key to match samples being used for eQTL 
eQTL_key <-  geno_rna_key %>% 
  filter(Genotyping_sample_id %in% gsub(".CEL", "", geno_samples_eQTL$FID))
dim(eQTL_key) # [1] 3394    5
head(eQTL_key)
# RNA_sanger_sample_id   ANON_ID Sample_taken_at Genotyping_sample_id id_run_position
# 1        pops212971188 OBGYN0830        12 weeks            S00366301         45759_1
# 2        pops212971189 OBGYN0850        12 weeks            S00367194         45759_1
# 3        pops212971190 OBGYN0714        20 weeks            S00356904         45759_1
# 4        pops212971191 OBGYN0737        20 weeks            S00358706         45759_1
# 5        pops212971192 OBGYN1037        20 weeks            S00379035         45759_1
# 6        pops212971196 OBGYN0831        12 weeks            S00366050         45759_1

# join covariates from PEER analysis to the sample key 
covariates <- eQTL_key %>% left_join(peer_covariates %>% select(-X20.weeks, -X28.weeks, -X36.weeks) %>% rownames_to_column(var = "RNA_sanger_sample_id"))
dim(covariates) # [1] 3394   13

# write this out 
write.csv(covariates %>% mutate(Sample_taken_at = gsub(" ", "_", Sample_taken_at)), 
          file = covariate_outpath, row.names = FALSE, quote = FALSE)

# subset count matrix and sample info df to only include samples from individuals with genotyping data we will use
counts[1:5, 1:5]
dim(counts) # [1] 18826  3776
counts.eQTL <- counts %>% select(eQTL_key$RNA_sanger_sample_id)
# remove HLA-DRB5 == ENSG00000198502
# remove HLA-DRB3 and HLA-DRB4
counts.eQTL <- counts.eQTL[!(row.names(counts.eQTL) %in% c("ENSG00000198502", "HLA-DRB3", "HLA-DRB4")),]
dim(counts.eQTL) # [1] 18823  3394
# write this out 
saveRDS(counts.eQTL, file = eQTL_counts_outpath)

################################################################################

# Make a list of SNP-gene pairs to test (SNPs 1Mb either side of TSS for each gene)

################################################################################

# prepare gene position information file
gtf <- as.data.frame(gtf)
gene.info.filt <- gtf[match(rownames(counts.eQTL), gtf$gene_id), ]
rownames(gene.info.filt) <- rownames(counts.eQTL)
all(gene.info.filt$gene_id == rownames(counts.eQTL))
gene.info.filt <- gene.info.filt[, c("seqnames", "start", "end", "strand", 
                                     "gene_id", "gene_name")]
gene.info.filt$seqnames <- as.character(gene.info.filt$seqnames)
gene.info.filt$seqnames[gene.info.filt$seqnames == "X"] <- 23
write.table(gene.info.filt, paste0(output_dir, "input_data/gene_info.txt"), sep = "\t")

# read in SNP and gene positional information
snp.info <- data.frame(data.table::fread(bim_inpath, 
                                         header=F, stringsAsFactors = F))
genes <- read.delim(paste0(output_dir, "input_data/gene_info.txt"), stringsAsFactors = F)

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
                                                      &  chr.bim.list[[c]][, 3] < (genes.list[[c]][i, 2] + 1000000)))] # position less than 1mb greater than gene start
    if(length(snp.names) > 1){
      results.c[[i]] <- data.frame(rownames(genes.list[[c]])[i], snp.names, 
                                   rep(as.numeric(genes.list[[c]][i, 1]), length(snp.names)))
    }
  }
  results[[c]] <- data.table::rbindlist(results.c)
  colnames(results[[c]]) <- c("Gene", "SNP", "Chr")
}

results.2 <- data.table::rbindlist(results)
results.2[1:5, 1:3] # dim 84746712        3
write.table(results.2, paste0(output_dir, "input_data/gene_snp_pairs_cis_rnaseq.txt"), sep="\t", 
            quote=FALSE, row.names=FALSE)
# list of snps 
snps <- unique(as.character(results.2$SNP))
write.table(snps, paste0(output_dir, "input_data/snps_for_cis_eqtl_rnaseq.txt"), sep="\t", quote=F, 
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
# cd genotyping/analysis/eQTL/input_data
# mkdir geno_files
# # extract the snps we will use and recode the files as raw format
# for CHR in {1..22} ; do plink --bfile Imputed_Genotyping_Data_byChr/chr.${CHR}.imputed.eqtl --extract snps_for_cis_eqtl_rnaseq.txt --recode A --out geno_files/genotyping_for_rna-seq_eqtl_${CHR} --allow-extra-chr --chr ${CHR} ; done
# plink --bfile Imputed_Genotyping_Data_byChr/chr.X.imputed.eqtl --extract snps_for_cis_eqtl_rnaseq.txt --recode A --out geno_files/genotyping_for_rna-seq_eqtl_X --allow-extra-chr --chr X

################################################################################
# assemble all files

# read in gene expression, covariates, sample key, peer factors and make sure they are all in the same order
covariates <- read.csv(covariate_outpath) %>% 
  column_to_rownames(var = "RNA_sanger_sample_id")
gene_exp <- readRDS(eQTL_counts_outpath)
all(colnames(gene_exp) == rownames(covariates))
geno_id_CEL <- paste0(covariates$Genotyping_sample_id, ".CEL")
individual <- covariates$ANON_ID
batch <- covariates$id_run_position
timepoint <- as.factor(covariates$Sample_taken_at)
# peer factors from https://github.com/PMBio/peer/wiki/Tutorial
# The output is written to directory peer_out by default, creating csv files for the residuals after accounting for the factors (residuals.csv, NxG matrix), 
#       the inferred factors (X.csv, NxK), the weights of each factor for every gene (W.csv, GxK), 
#       and the inverse variance of the weights (Alpha.csv, Kx1).
peer_factors <- t(read.csv(peer_inpath, header = FALSE))
colnames(peer_factors) <- c("20_weeks","28_weeks","36_weeks","Neutrophils_centered","Monocytes_centered","geno_1","geno_2","geno_3","geno_4","geno_5","geno_6", "1s", paste0("peer_", c(1:50)))
peer_factors <- peer_factors %>% as.data.frame() %>% select(contains("peer_"))
rownames(peer_factors) <- rownames(covariates)
peer_factors <- as.matrix(peer_factors)
# remove genotyping ID, ANON_ID and id_run_poisiton, sampletakenat from covariates
covariates <- covariates %>% select(-Genotyping_sample_id, -ANON_ID, -id_run_position, -Sample_taken_at) %>% as.matrix()

# set up genotyping files
for(i in c(1:22, "X")){
  print(i)
  geno <- data.frame(data.table::fread(paste(geno_files_inpath, "genotyping_for_rna-seq_eqtl_", i, ".raw", sep=""),
                                       sep=" ", drop = 2:6))
  m1 <- match(geno_id_CEL, geno$FID)
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
  
  eQTL_files <- paste(output_dir, "input_data/eqtl_files_", i, ".rda", sep="")
  save(list=c("gene_exp", "geno", "covariates", "individual", "batch", "peer_factors", "timepoint"),
       file = eQTL_files)
}

# make an outcomes rda to read in in outcome interactions
outcomes <- read.csv(sample_info_inpath) %>% 
  mutate(pn_gdm = as.character(as.factor(pn_diabetes_3cat == 2)),
         sga = as.character(as.factor(BW_Centile_Br1990 < 10)),
         pn_petACOG13 = as.character(as.factor(pn_petACOG13 == 1)), 
         pn_ptd_sp = as.character(as.factor(pn_ptd_sp == 1))) %>%
  select(pn_gdm, sga, pn_petACOG13, pn_ptd_sp, RNA_sanger_sample_id) %>% 
  column_to_rownames(var = "RNA_sanger_sample_id")

# order these the same as the eQTL input files 
outcomes <- outcomes[rownames(covariates),]
all(rownames(outcomes) == rownames(covariates))
gdm <- outcomes$pn_gdm
pe <- outcomes$pn_petACOG13
sga <- outcomes$sga
ptb <- outcomes$pn_ptd_sp

# write this out as an rdata object
outcome_files <- paste0(output_dir, "input_data/eqtl_outcome_data.rda")
save(list=c("gdm", "pe", "sga", "ptb"),
     file = outcome_files)

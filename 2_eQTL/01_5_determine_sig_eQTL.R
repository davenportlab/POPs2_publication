# 01_5_determine_sig_eQTL.R 

################################################################################

# 1.5. Hierarchical multiple testing control: identify significant cis-eQTL from 1st pass

################################################################################

# Aim: Identify significance thresholds and significant main effects eQTL

########################### Output paths ##########################

########################### Input paths ###########################
bim_filepath <- "genotyping/analysis/eQTL/input_data/POPS2_genotyping_imputed_eQTL_rsid_refilter.bim"
gene_info_filepath <- "genotyping/analysis/eQTL/input_data/gene_info.txt"
filepath <- "genotyping/analysis/eQTL/"

########################### Parameters ############################
options(stringsAsFactors = FALSE)

########################### Load packages ###########################
library(data.table)

########################### Analysis ###########################

# combine results from first pass mapping by chromosome
cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, full.names = T)
  print(length(temp_eqtl))
  eqtl <- lapply(temp_eqtl, readRDS)
  eqtl <- rbindlist(eqtl)
  eqtl <- data.frame("snps" = eqtl$SNP, "gene" = eqtl$Gene,
                     "statistic"=eqtl$eQTL_t, "pvalue" = eqtl$eQTL_pval,
                     "beta"=eqtl$eQTL_beta, 'se'=eqtl$eQTL_SE)
  eqtl$snps <- gsub(":", ".", eqtl$snps)
  print(dim(eqtl))
  write.table(eqtl, output_filename, sep="\t", row.names=F, quote=F)
}

# make files for eigenMT
# This bim file is created by filtering the full plink dataset 
# (genotyping/data/imputed_genotyping_dataset_for_eQTL/POPS2_genotyping_imputed_eQTL_rsid)
# to contain only SNPs used for eQTL mapping (genotyping/analysis/eQTL/input_data/snps_for_cis_eqtl_rnaseq.txt)
bim <- data.frame(fread(bim_filepath,
                        header=FALSE))

gene.info <- read.delim(gene_info_filepath)
genes <- gene.info[, 1:3]
colnames(genes) <- c("chrom", "TSS", "TES")
genes$TSS[which(gene.info$strand == "-")] <- gene.info$end[which(gene.info$strand == "-")]
genes$TES[which(gene.info$strand == "-")] <- gene.info$start[which(gene.info$strand == "-")]
write.table(genes, paste0(filepath, "input_data/gene_pos_for_eigenMT.txt"))

for(i in c(1:23)){
  chr = ifelse(i == 23, "X", i)
  print(paste("start", i))
  # eqtl results
  cat.rds(paste0(filepath, "output_data/main_effects/chr", chr), "*rds",
          paste0(filepath, "output_data/eigenMT/eqtl_results_for_eigenMT_chr", chr, ".txt"))
  
  res <- data.frame(fread(paste0(filepath, "output_data/eigenMT/eqtl_results_for_eigenMT_chr", chr, ".txt")))
  
  # genotyping data
  load(paste0(filepath, "input_data/eqtl_files_", chr, ".rda"))
  geno <- as.data.table(geno)
  geno <- geno[which(!duplicated(individual))]
  geno <- t(geno)
  colnames(geno) <- individual[which(!duplicated(individual))]
  write.table(geno, paste0(filepath, "output_data/eigenMT/geno_for_eigenMT_chr", chr, ".txt"),
              sep="\t", quote=F)
  
  # snp positions
  chr.bim <- subset(bim, V1 == i)
  chr.bim$V2 <- gsub(":", ".", chr.bim$V2)
  geno.pos <- chr.bim[match(res$snps, chr.bim$V2), c(2, 1, 4)]
  write.table(geno.pos, paste0(filepath, "output_data/eigenMT/geno_pos_for_eigenMT_chr", chr, ".txt"), sep="\t",
              row.names=F, quote=F)
  
  # gene positions
  genes.pos <- subset(genes, chrom == i)
  genes.pos <- genes.pos[match(unique(res$gene), rownames(genes.pos)), ]
  write.table(genes.pos, paste0(filepath, "output_data/eigenMT/gene_pos_for_eigenMT_chr", chr, ".txt"), 
              sep="\t", quote=F)
  
  # add to res file
  res$chr <- i
  res$SNPpos <- geno.pos[match(res$snps, geno.pos$V2), 3]
  res$TSS <- genes.pos$TSS[match(res$gene, rownames(genes.pos))]
  write.table(res, paste0(filepath, "output_data/eigenMT/eqtl_results_for_eigenMT_chr", chr, ".txt"),
              sep="\t", quote=F, row.names=F)
  print(paste("end", i))
}

# Run eigenMT to determine local significance - bash 
###############################################################################
# #!/bin/bash
# #BSUB -o logfiles/eigenMT.log
# #BSUB -e logfiles/eigenMT.err
# #BSUB -J eigenMT
# #BSUB -q "normal"
# #BSUB -R "select[mem>5000] rusage[mem=5000] span[hosts=1]"
# #BSUB -M5000
# 
# # cd to correct dir
# cd genotyping/analysis/eQTL/output_data/eigenMT
# # load module
# module load HGI/softpack/groups/team282/eigenMT/1
# 
# for CHR in {1..22}
# do
# echo "$CHR"
# python /software/team282/cs54/eigenMT/eigenMT_fix.py --CHROM $CHR --QTL eqtl_results_for_eigenMT_chr$CHR.txt --GEN geno_for_eigenMT_chr$CHR.txt --GENPOS geno_pos_for_eigenMT_chr$CHR.txt --OUT eigenMT_chr$CHR --window 200 --PHEPOS gene_pos_for_eigenMT_chr$CHR.txt --cis_dist 1000000
# done
# 
# # do chr X
# echo "X"
# python /software/team282/cs54/eigenMT/eigenMT_fix.py --CHROM 23 --QTL eqtl_results_for_eigenMT_chrX.txt --GEN geno_for_eigenMT_chrX.txt --GENPOS geno_pos_for_eigenMT_chrX.txt --OUT eigenMT_chrX --window 200 --PHEPOS gene_pos_for_eigenMT_chrX.txt --cis_dist 1000000
###############################################################################

# Hierarchical correction to identify eGenes

# read in lead results and correct for number of genes
eigen.res <- list.files(path = paste0(filepath, "output_data/eigenMT"), pattern = "^eigenMT_chr*", full.names = T)
eigen.res <- rbindlist(lapply(eigen.res, read.delim))
eigen.res$BF.FDR <- p.adjust(eigen.res$BF, method="fdr")
eigen.res$Sig <- eigen.res$BF.FDR < 0.05
write.table(eigen.res, paste0(filepath, "output_data/eigenMT/cis-eQTL_eigenMT_corrected.txt"), sep="\t")

################################################################################

# Identify all significant associations for eGenes

# calculate the global significance threshold:
# the locally corrected p value corresponding to q-value=0.05
ub <- min(eigen.res$BF[eigen.res$BF.FDR > 0.05])  # smallest p-value above FDR
lb <- max(eigen.res$BF[eigen.res$BF.FDR <= 0.05])  # largest p-value below FDR
fdr.threshold <- (lb+ub)/2

# determine the nominal pvalue threshold for each gene
thresholds <- rbindlist(lapply(1:nrow(eigen.res), function(g){
  gene <- as.character(eigen.res[g, 2])
  n.tests <- eigen.res$TESTS[g]
  threshold <- fdr.threshold/n.tests
  return(data.frame(gene, n.tests, threshold))
}))
write.table(thresholds, paste0(filepath, "output_data/eigenMT/nominal_pval_thresholds.txt"), sep="\t",
            row.names=F, quote=F)
eigen.res$threshold <- thresholds$threshold
write.table(eigen.res, paste0(filepath, "output_data/eigenMT/ciseqtl_eigenMT_corrected.txt"), sep="\t")

full.res <- list.files(path = paste0(filepath, "output_data/eigenMT"),
                       pattern = "^eqtl_results_for_eigenMT_chr",
                       full.names = T)

full.res <- rbindlist(lapply(full.res, read.delim, stringsAsFactors=F))
full.res$threshold <- eigen.res$threshold[match(full.res$gene, eigen.res$gene)]
saveRDS(full.res, paste0(filepath, "output_data/main_effects/ciseqtl_all.rds"))

all.sig.eqtl <- full.res[which(full.res$pvalue <= full.res$threshold), ]
saveRDS(all.sig.eqtl, paste0(filepath, "output_data/main_effects/cisqtl_all_significant.rds"))
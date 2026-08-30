# 04_6_determine_sig_interval.R 

################################################################################

# 4.6. Hierarchical multiple testing control: identify significant interval cis-eQTL from 1st pass 

################################################################################

# Aim: Identify significance thresholds and significant main effects eQTL for interval FRA cohort

########################### Paths ##########################
filepath <- "interval_rna/POPS2_comparison/"

########################### Parameters ############################
options(stringsAsFactors = FALSE)

########################### Load packages ###########################
library(data.table)
library(tidyverse)
library(UpSetR)

########################### Analysis ###########################

cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, full.names = T)
  print(paste("Number of files:", length(temp_eqtl)))
  eqtl <- lapply(temp_eqtl, readRDS)
  eqtl <- rbindlist(eqtl)
  print(c("eQTL results:", dim(eqtl)))
  write.table(eqtl, output_filename, sep="\t", row.names=F, quote=F)
}

# make files for eigenMT
bim <- data.frame(fread(paste0(filepath, "data_subset/genotyping/INTERVAL_RNAseq_imputed_b38_AllChromosomes_female_u45_filt_maf0.05_snps.bim"),
                        header=FALSE))

gene.info <- read.delim(paste0(filepath, "eQTL_input_data/gene_info_update_10_27_25.txt"))
genes <- gene.info[, 1:3]
colnames(genes) <- c("chrom", "TSS", "TES")
genes$TSS[which(gene.info$strand == "-")] <- gene.info$end[which(gene.info$strand == "-")]
genes$TES[which(gene.info$strand == "-")] <- gene.info$start[which(gene.info$strand == "-")]
write.table(genes, paste0(filepath, "eQTL_input_data/interval_gene_pos_for_eigenMT.txt"))

for(i in c(1:23)){
  chr = ifelse(i == 23, "X", i)
  print(paste("start", i))
  # eqtl results
  cat.rds(paste0(filepath, "eQTL_main_effects/chr", chr), "*rds",
          paste0(filepath, "eQTL_main_effects/eqtl_results_chr", chr, ".txt"))
  
  res <- data.frame(fread(paste0(filepath, "eQTL_main_effects/eqtl_results_chr", chr, ".txt")))
  
  # genotyping data
  load(paste0(filepath, "eQTL_input_data/eqtl_files_", chr, ".rda"))
  geno <- as.data.table(geno)
  geno <- t(geno)
  colnames(geno) <- colnames(gene_exp)
  write.table(geno, paste0(filepath, "eigenMT/interval_geno_for_eigenMT_chr", chr, ".txt"),
              sep="\t", quote=F)
  
  # snp positions
  chr.bim <- subset(bim, V1 == i)
  #chr.bim$V2 <- gsub(":", ".", chr.bim$V2)
  geno.pos <- chr.bim[match(res$SNP, chr.bim$V2), c(2, 1, 4)]
  write.table(geno.pos, paste0(filepath, "eigenMT/interval_geno_pos_for_eigenMT_chr", chr, ".txt"), sep="\t",
              row.names=F, quote=F)
  
  # gene positions
  genes.pos <- subset(genes, chrom == i)
  genes.pos <- genes.pos[match(unique(res$Gene), rownames(genes.pos)), ]
  write.table(genes.pos, paste0(filepath, "eigenMT/interval_gene_pos_for_eigenMT_chr", chr, ".txt"), 
              sep="\t", quote=F)
  
  # add to res file
  res <- data.frame("snps" = res$SNP, "gene" = res$Gene,
                    "statistic"=res$eQTL_t, "pvalue" = res$eQTL_pval,
                    "beta"=res$eQTL_beta, 'se'=res$eQTL_SE)
  res$chr <- i
  res$SNPpos <- geno.pos[match(res$snps, geno.pos$V2), 3]
  res$TSS <- genes.pos$TSS[match(res$gene, rownames(genes.pos))]
  write.table(res, paste0(filepath, "eigenMT/interval_eqtl_results_for_eigenMT_chr", chr, ".txt"),
              sep="\t", quote=F, row.names=F)
  print(paste("end", i))
}

# Run eigenMT to determine local significance - bash 
# make commands to run eigenMT per chr 
eigenMT_cmd <- paste0("python /path/to/eigenMT_fix.py --CHROM ", c(1:23) , " --QTL interval_eqtl_results_for_eigenMT_chr", c(1:22, "X"), ".txt --GEN interval_geno_for_eigenMT_chr", c(1:22, "X"), ".txt --GENPOS interval_geno_pos_for_eigenMT_chr", c(1:22, "X"), ".txt --OUT eigenMT_chr", c(1:22, "X"), " --window 200 --PHEPOS interval_gene_pos_for_eigenMT_chr", c(1:22, "X"), ".txt --cis_dist 1000000")
write.table(eigenMT_cmd, file = paste0("/path/to/interval_rna/POPS2_comparison/eigenMT/interval_eigenMT_wr_cmds.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

###############################################################################
#                     Run in bash 
###############################################################################
# submit wr jobs
# cd /path/to/interval_rna/POPS2_comparison/eigenMT
# module load HGI/common/wr
# wr add -f interval_eigenMT_wr_cmds.txt -i interval_eigenMT -r 0 --cwd_matters --queue "normal" --memory "5G" -o 2 --modules HGI/softpack/groups/team282/eigenMT/1 --rerun
###############################################################################

# Hierarchical correction to identify eGenes

# read in lead results and correct for number of genes
eigen.res <- list.files(path = paste0(filepath, "eigenMT"), pattern = "^eigenMT_chr*", full.names = T)
eigen.res <- rbindlist(lapply(eigen.res, read.delim))
eigen.res$BF.FDR <- p.adjust(eigen.res$BF, method="fdr")
eigen.res$Sig <- eigen.res$BF.FDR < 0.05
write.table(eigen.res, paste0(filepath, "eigenMT/interval_cis-eQTL_eigenMT_corrected.txt"), sep="\t")

# Identify all significant associations for eGenes

# calculate the global significance threshold:
# the locally corrected p value corresponding to q-value=0.05
ub <- min(eigen.res$BF[eigen.res$BF.FDR > 0.05])  # smallest p-value above FDR
lb <- max(eigen.res$BF[eigen.res$BF.FDR <= 0.05])  # largest p-value below FDR
fdr.threshold <- (lb+ub)/2 # 0.02666985

# determine the nominal pvalue threshold for each gene
thresholds <- rbindlist(lapply(1:nrow(eigen.res), function(g){
  gene <- as.character(eigen.res[g, 2])
  n.tests <- eigen.res$TESTS[g]
  threshold <- fdr.threshold/n.tests
  return(data.frame(gene, n.tests, threshold))
}))
write.table(thresholds, paste0(filepath, "eigenMT/interval_nominal_pval_thresholds.txt"), sep="\t",
            row.names=F, quote=F)
eigen.res$threshold <- thresholds$threshold
write.table(eigen.res, paste0(filepath, "eigenMT/interval_ciseqtl_eigenMT_corrected.txt"), sep="\t")

full.res <- list.files(path = paste0(filepath, "eigenMT"),
                       pattern = "^interval_eqtl_results_for_eigenMT_chr",
                       full.names = T)
full.res <- rbindlist(lapply(full.res, read.delim, stringsAsFactors=F))
full.res$threshold <- eigen.res$threshold[match(full.res$gene, eigen.res$gene)]
saveRDS(full.res, paste0(filepath, "eQTL_main_effects/interval_ciseqtl_all.rds"))

all.sig.eqtl <- full.res[which(full.res$pvalue <= full.res$threshold), ]
saveRDS(all.sig.eqtl, paste0(filepath, "eQTL_main_effects/interval_cisqtl_all_significant.rds"))

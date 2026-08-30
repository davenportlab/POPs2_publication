# 05_4_determine_sig_single_tp.R

################################################################################

# 5.4. Hierarchical multiple testing control: identify significant cis-eQTL from 1st pass

################################################################################

# Aim: Identify significance thresholds and significant main effects eQTL for POPs2 TP cohort

########################### Paths ##########################
filepath <- "genotyping/analysis/eQTL/"

########################### Parameters ############################
options(stringsAsFactors = FALSE)

########################### Load packages ###########################
library(data.table)
library(tidyverse)

########################### Analysis ###########################

cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, full.names = T)
  print(paste("Number of files:", length(temp_eqtl)))
  eqtl <- lapply(temp_eqtl, readRDS)
  eqtl <- rbindlist(eqtl)
  print(c("eQTL results:", dim(eqtl)))
  write.table(eqtl, output_filename, sep="\t", row.names=F, quote=F)
}

# read in the new bim file
bim <- data.frame(fread(paste0(filepath, "input_data/single_tp_eQTL/POPS2_genotyping_imputed_eQTL_rsid_fourRNA_filt_maf0.05_snps.bim"),
                        header=FALSE))

gene.info <- read.delim(paste0(filepath, "input_data/single_tp_eQTL/gene_info_fourRNA.txt"))
genes <- gene.info[, 1:3]
colnames(genes) <- c("chrom", "TSS", "TES")
genes$TSS[which(gene.info$strand == "-")] <- gene.info$end[which(gene.info$strand == "-")]
genes$TES[which(gene.info$strand == "-")] <- gene.info$start[which(gene.info$strand == "-")]
write.table(genes, paste0(filepath, "input_data/single_tp_eQTL/gene_pos_for_eigenMT.txt"))

# for each time point, compile results
for(wk in c("12_weeks", "20_weeks", "28_weeks", "36_weeks")){
  for(i in 1:23){
    chr = ifelse(i == 23, "X", i)
    print(paste("start", wk, i))
    print("combine eqtl results by chromosome")
    cat.rds(paste0(filepath, "output_data/single_tp_analysis/", wk, "/chr", chr), "*rds",
            paste0(filepath, "output_data/single_tp_analysis/", wk, "/", wk, "_eqtl_results_chr", chr, ".txt"))
    print("read those eQTL results back in")
    res <- data.frame(fread(paste0(filepath, "output_data/single_tp_analysis/", wk, "/", wk, "_eqtl_results_chr", chr, ".txt")))
    print("res in")
    # only need to re-write the genotyping data once since it is the same across tps
    if(wk == "12_weeks"){
      print("set up data for eigenMT")
      # genotyping data
      load(paste0(filepath, "input_data/single_tp_eQTL/eqtl_files_", chr, ".rda"))
      rm(gene_exp)
      rm(covariates_sub)
      geno <- as.data.table(geno)
      geno <- geno[which(!duplicated(individual))] # subset to one sample per individual
      geno <- t(geno)
      colnames(geno) <- individual[which(!duplicated(individual))]
      write.table(geno, paste0(filepath, "input_data/single_tp_eQTL/geno_for_eigenMT_chr", chr, ".txt"),
                  sep="\t", quote=F)
      
      # snp positions
      chr.bim <- subset(bim, V1 == i)
      chr.bim$V2 <- gsub(":", ".", chr.bim$V2)
      geno.pos <- chr.bim[match(res$SNP, chr.bim$V2), c(2, 1, 4)]
      write.table(geno.pos, paste0(filepath, "input_data/single_tp_eQTL/geno_pos_for_eigenMT_chr", chr, ".txt"), sep="\t",
                  row.names=F, quote=F)
      
      # gene positions
      genes.pos <- subset(genes, chrom == i)
      genes.pos <- genes.pos[match(unique(res$Gene), rownames(genes.pos)), ]
      write.table(genes.pos, paste0(filepath, "input_data/single_tp_eQTL/gene_pos_for_eigenMT_chr", chr, ".txt"),
                  sep="\t", quote=F)
    } else { # otherwise read in the info we need
      print("read in info for eigenMT")
      geno.pos <- read.table(paste0(filepath, "input_data/single_tp_eQTL/geno_pos_for_eigenMT_chr", chr, ".txt"),
                             header = TRUE)
      genes.pos <- read.table(paste0(filepath, "input_data/single_tp_eQTL/gene_pos_for_eigenMT_chr", chr, ".txt"),
                              header = TRUE)
    }
    print("add info to res files")
    # add info to res file
    res <- data.frame("snps" = res$SNP, "gene" = res$Gene,
                      "statistic"=res$eQTL_t, "pvalue" = res$eQTL_pval,
                      "beta"=res$eQTL_beta, 'se'=res$eQTL_SE)
    res$chr <- i
    res$SNPpos <- geno.pos[match(res$snps, geno.pos$V2), 3]
    res$TSS <- genes.pos$TSS[match(res$gene, rownames(genes.pos))]
    print("write results as table")
    write.table(res, paste0(filepath, "output_data/single_tp_analysis/", wk, "/eigenMT/eqtl_results_for_eigenMT_chr", chr, ".txt"),
                sep="\t", quote=F, row.names=F)
    print(paste("end", i))
    
  }
  print(paste("end", wk))
}

# Run eigenMT to determine local significance - bash
for(wk in c("12_weeks", "20_weeks", "28_weeks", "36_weeks")){
  eigenMT_cmd <- paste0("python eigenMT_fix.py --CHROM ", c(1:23) , " --QTL ", wk, "/eigenMT/eqtl_results_for_eigenMT_chr", c(1:22, "X"), ".txt --GEN ../../input_data/single_tp_eQTL/geno_for_eigenMT_chr", c(1:22, "X"), ".txt --GENPOS ../../input_data/single_tp_eQTL/geno_pos_for_eigenMT_chr", c(1:22, "X"), ".txt --OUT ", wk, "/eigenMT/", wk, "_eigenMT_chr", c(1:22, "X"), " --window 200 --PHEPOS ../../input_data/single_tp_eQTL/gene_pos_for_eigenMT_chr", c(1:22, "X"), ".txt --cis_dist 1000000")
  write.table(eigenMT_cmd, file = paste0("/path/to/dir/genotyping/analysis/eQTL/output_data/single_tp_analysis/wr_scripts/", wk, "_eigenMT_wr_cmds.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

###############################################################################
#               Run in bash 
###############################################################################
# submit wr jobs
# cd genotyping/analysis/eQTL/output_data/single_tp_analysis
# module load HGI/common/wr
# module load HGI/softpack/groups/team282/eigenMT/1
# wr add -f wr_scripts/12_weeks_eigenMT_wr_cmds.txt -i 12wk_eigenMT -r 0 --cwd_matters --queue "normal" --memory "5G" -o 2 --modules HGI/softpack/groups/team282/eigenMT/1 --rerun
# wr status -i 12wk_eigenMT
#
# wr add -f wr_scripts/20_weeks_eigenMT_wr_cmds.txt -i 20wk_eigenMT -r 0 --cwd_matters --queue "normal" --memory "5G" -o 2 --modules HGI/softpack/groups/team282/eigenMT/1 --rerun
# wr status -i 20wks_eigenMT
#
# wr add -f wr_scripts/28_weeks_eigenMT_wr_cmds.txt -i 28wk_eigenMT -r 0 --cwd_matters --queue "normal" --memory "5G" -o 2 --modules HGI/softpack/groups/team282/eigenMT/1 --rerun
# wr status -i 28wks_eigenMT
#
# wr add -f wr_scripts/36_weeks_eigenMT_wr_cmds.txt -i 36wk_eigenMT -r 0 --cwd_matters --queue "normal" --memory "5G" -o 2 --modules HGI/softpack/groups/team282/eigenMT/1 --rerun
# wr status -i 36wks_eigenMT
###############################################################################

# loop through timepoints to make final outputs
for(wk in c("12_weeks", "20_weeks", "28_weeks", "36_weeks")){
  # Hierarchical correction to identify eGenes
  print(paste("Starting", wk))
  print("read in lead results and correct for number of genes")
  eigen.res <- list.files(path = paste0(filepath, "output_data/single_tp_analysis/", wk, "/eigenMT"), pattern = paste0("^", wk, "_eigenMT_chr*"), full.names = T)
  eigen.res <- rbindlist(lapply(eigen.res, read.delim))
  eigen.res$BF.FDR <- p.adjust(eigen.res$BF, method="fdr")
  eigen.res$Sig <- eigen.res$BF.FDR < 0.05
  print("start write out lead results")
  write.table(eigen.res, paste0(filepath, "output_data/single_tp_analysis/", wk, "/eigenMT/", wk, "_cis-eQTL_eigenMT_corrected.txt"), sep="\t")
  print("done write out lead results")
  # Identify all significant associations for eGenes
  
  # calculate the global significance threshold:
  # the locally corrected p value corresponding to q-value=0.05
  ub <- min(eigen.res$BF[eigen.res$BF.FDR > 0.05])  # smallest p-value above FDR
  lb <- max(eigen.res$BF[eigen.res$BF.FDR <= 0.05])  # largest p-value below FDR
  fdr.threshold <- (lb+ub)/2
  print(paste("FDR threshold", fdr.threshold))
  # 12 weeks =  0.02241992
  
  print("determine the nominal pvalue threshold for each gene")
  thresholds <- rbindlist(lapply(1:nrow(eigen.res), function(g){
    gene <- as.character(eigen.res[g, 2])
    n.tests <- eigen.res$TESTS[g]
    threshold <- fdr.threshold/n.tests
    return(data.frame(gene, n.tests, threshold))
  }))
  print("write out thresholds")
  write.table(thresholds, paste0(filepath, "output_data/single_tp_analysis/", wk, "/eigenMT/", wk, "_nominal_pval_thresholds.txt"), sep="\t",
              row.names=F, quote=F)
  eigen.res$threshold <- thresholds$threshold
  print("write out corrected ciseQTL")
  write.table(eigen.res, paste0(filepath, "output_data/single_tp_analysis/", wk, "/eigenMT/", wk, "_ciseqtl_eigenMT_corrected.txt"), sep="\t")
  print("read in full results")
  full.res <- list.files(path = paste0(filepath, "output_data/single_tp_analysis/", wk, "/eigenMT"),
                         pattern = paste0("^eqtl_results_for_eigenMT_chr"),
                         full.names = T)
  full.res <- rbindlist(lapply(full.res, fread))
  print("add thresholds to full results")
  full.res$threshold <- eigen.res$threshold[match(full.res$gene, eigen.res$gene)]
  print("save full results")
  saveRDS(full.res, paste0(filepath, "output_data/single_tp_analysis/", wk, "/", wk, "_ciseqtl_all.rds"))
  print("save significant results")
  all.sig.eqtl <- full.res[which(full.res$pvalue <= full.res$threshold), ]
  saveRDS(all.sig.eqtl, paste0(filepath, "output_data/single_tp_analysis/", wk, "/", wk, "_cisqtl_all_significant.rds"))
  print("done")
  
}

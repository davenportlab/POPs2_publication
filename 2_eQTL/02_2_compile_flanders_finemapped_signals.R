# 02_2_compile_flanders_finemapped_signals.R 

################################################################################

# 2.2. Compile fine-mapped signals from Flanders

################################################################################

# Aim: Compile finemapped loci for each study of interest from Flanders

########################### Output paths ##########################
path <- "genotyping/analysis/eQTL/"

########################### Load packages ###########################
library(tidyverse)
library(data.table)

########################### Analysis ###########################

##################### POPs2 ##################### 

pops_finemap_files <- list.files(path=paste0(path, "output_data/flanders_5/flanders_output/results/finemap"), pattern = "^POPs2_bulk.*\\.rds$", full.names = T)

pops_finemapped.loci <- lapply(pops_finemap_files, function(f) {
  x <- readRDS(f)
  
  rbindlist(lapply(names(x), function(cs_name) {
    cs <- x[[cs_name]]
    
    fm <- cs$finemapping_lABFs
    eff <- cs$effect
    meta <- cs$metadata
    
    ## compute p-values exactly as in flanders https://github.com/Biostatistics-Unit-HT/Flanders/blob/255d8cb92acebde818ba0f8ae9b0cd56237cac4a/bin/s07_rds2anndata.R#L254
    p_values <- pchisq((fm$bC / fm$bC_se)^2, df = 1, lower.tail = FALSE)
    top_pvalue <- min(p_values, na.rm = TRUE)
    
    data.table(
      credible_set_name = cs_name,
      credible_set_snps = paste(fm[is_cs == TRUE, snp], collapse = ","),
      
      study_id = meta$study_id,
      phenotype_id = meta$phenotype_id,
      chr = meta$chr,
      start = meta$start,
      end = meta$end,
      
      ## top p-value per credible set
      top_pvalue = top_pvalue,
      
      path_rds = f,
      
      snp = eff$snp,
      a1 = eff$a1,
      a0 = eff$a0,
      freq = eff$freq,
      N = eff$N,
      bC = eff$beta,
      bC_se = eff$se
      
    )
  }))
})

pops_finemapped.loci <- as.data.frame(rbindlist(pops_finemapped.loci))
dim(pops_finemapped.loci)

# write this out as a rda file
saveRDS(pops_finemapped.loci, paste0(path, "output_data/flanders_finemap_2/finemapped_loci_pops.rds"))

##################### GWAS ##################### 

GWAS_full_finemap_files <- list.files(path=paste0(path, "output_data/flanders_5/flanders_output/results/finemap"), pattern = "GWAS", full.names = T)

gwas_finemapped.loci <- lapply(GWAS_full_finemap_files, function(f) {
  x <- readRDS(f)
  
  rbindlist(lapply(names(x), function(cs_name) {
    cs <- x[[cs_name]]
    
    fm <- cs$finemapping_lABFs
    eff <- cs$effect
    meta <- cs$metadata
    
    ## compute p-values exactly as in flanders https://github.com/Biostatistics-Unit-HT/Flanders/blob/255d8cb92acebde818ba0f8ae9b0cd56237cac4a/bin/s07_rds2anndata.R#L254
    p_values <- pchisq((fm$bC / fm$bC_se)^2, df = 1, lower.tail = FALSE)
    top_pvalue <- min(p_values, na.rm = TRUE)
    
    data.table(
      credible_set_name = cs_name,
      credible_set_snps = paste(fm[is_cs == TRUE, snp], collapse = ","),
      
      study_id = meta$study_id,
      phenotype_id = meta$phenotype_id,
      chr = meta$chr,
      start = meta$start,
      end = meta$end,
      
      ## top p-value per credible set
      top_pvalue = top_pvalue,
      
      path_rds = f,
      
      snp = eff$snp,
      a1 = eff$a1,
      a0 = eff$a0,
      freq = eff$freq,
      N = eff$N,
      bC = eff$beta,
      bC_se = eff$se
      
    )
  }))
})

gwas_finemapped.loci <- as.data.frame(rbindlist(gwas_finemapped.loci))
dim(gwas_finemapped.loci)

# write this out as a rda file
saveRDS(gwas_finemapped.loci, paste0(path, "output_data/flanders_finemap_2/finemapped_loci_gwas.rds"))

##################### Interval ##################### 

Interval_full_finemap_files <- list.files(path=paste0(path, "output_data/flanders_5/flanders_output/results/finemap"), pattern = "^Interval_full_eQTL.*\\.rds$", full.names = T)

Interval_finemapped.loci <- lapply(Interval_full_finemap_files, function(f) {
  x <- readRDS(f)
  
  rbindlist(lapply(names(x), function(cs_name) {
    cs <- x[[cs_name]]
    
    fm <- cs$finemapping_lABFs
    eff <- cs$effect
    meta <- cs$metadata
    
    ## compute p-values exactly as in flanders https://github.com/Biostatistics-Unit-HT/Flanders/blob/255d8cb92acebde818ba0f8ae9b0cd56237cac4a/bin/s07_rds2anndata.R#L254
    p_values <- pchisq((fm$bC / fm$bC_se)^2, df = 1, lower.tail = FALSE)
    top_pvalue <- min(p_values, na.rm = TRUE)
    
    data.table(
      credible_set_name = cs_name,
      credible_set_snps = paste(fm[is_cs == TRUE, snp], collapse = ","),
      
      study_id = meta$study_id,
      phenotype_id = meta$phenotype_id,
      chr = meta$chr,
      start = meta$start,
      end = meta$end,
      
      ## top p-value per credible set
      top_pvalue = top_pvalue,
      
      path_rds = f,
      
      snp = eff$snp,
      a1 = eff$a1,
      a0 = eff$a0,
      freq = eff$freq,
      N = eff$N,
      bC = eff$beta,
      bC_se = eff$se
      
    )
  }))
})

Interval_finemapped.loci <- as.data.frame(rbindlist(Interval_finemapped.loci))
dim(Interval_finemapped.loci) # [1] 46006    16
# from interval paper: Stepwise conditional analyses for each cis-eQTL revealed 56,959 independent signals

# write this out as a rda file
saveRDS(Interval_finemapped.loci, paste0(path, "output_data/flanders_finemap_2/finemapped_loci_interval_full.rds"))
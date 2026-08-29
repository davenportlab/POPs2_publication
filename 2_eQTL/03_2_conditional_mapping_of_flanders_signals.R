# 03_2_conditional_mapping_of_flanders_signals.R

################################################################################

# 3.2. Conditionally independent eQTL mapping 

################################################################################

# Aim: Identify conditionally independent finemapped signals from Flanders

########################### Output paths ##########################
dir <- "genotyping/analysis/eQTL/"

########################### Input paths ###########################

########################### Parameters ############################
# n peer factors
n.peer <- 30 

########################### Load packages ###########################
library(lme4)
library(data.table)
library(tidyverse)

########################### Load data ###########################
signals_for_conditional_mapping <- readRDS(paste0(dir, "input_data/conditional_eQTL/gene_snp_pairs_for_cond_mapping.rds")) %>%
  arrange(chr, phenotype_id, signal_rank)
# read in p-value thresholds
thresholds <- read.delim(paste0(dir, "output_data/eigenMT/nominal_pval_thresholds.txt"))
# get list of genes 
genes <- unique(signals_for_conditional_mapping$phenotype_id)

########################### Analysis ###########################

gene_chr_map <- signals_for_conditional_mapping %>%
  distinct(phenotype_id, chr) %>%
  #rename("phenotype_id" = "Gene", "chr" = "Chr")
  rename("Gene" = "phenotype_id", "Chr" = "chr")

capture_warnings_and_output <- function(expr) {
  # Initialize empty lists to store warnings and output
  warnings_list <- list()
  
  # Capture output
  output <- capture.output(
    result <- withCallingHandlers(
      expr,
      warning = function(w) {
        warnings_list <<- c(warnings_list, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )
  
  # Return captured result, warnings, and output
  list(
    result = result,
    output = output,
    warnings = warnings_list
  )
}

# make function to update model if it fails to converge
update_model_fun <- function(model_with_warnings_obj){
  # save original model 
  model.original <- model_with_warnings_obj
  # get parameters
  pars <- getME(model_with_warnings_obj$result, "theta")
  # restart the fit from the original value 10 times
  for(u in 1:50){
    # set seed 
    set.seed(101)
    # update model
    model_with_warnings_obj <- capture_warnings_and_output(update(model.original$result, start=pars))
    # check whether it has worked 
    if(length(model_with_warnings_obj$warnings) > 0){
      pars <- getME(model_with_warnings_obj$result, "theta")
    }
    else{
      break
    }
  }
  # if this doesn't work, re-start from a slightly perturbed value 10 times 
  if(length(model_with_warnings_obj$warnings) > 0){
    # make tolerance
    strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
    # loop through and give 10 chances 
    for(u in 1:50){
      # change the parameters a bit 
      pars_x <- runif(length(pars),pars/1.01,pars*1.01)
      set.seed(101)
      # re-run 
      model_with_warnings_obj <- capture_warnings_and_output(update(model.original$result, start=pars_x,
                                                                    control=strict_tol))
      if(length(model_with_warnings_obj$warnings) > 0){
        pars <- getME(model_with_warnings_obj$result, "theta")
      }
      else{
        break
      }
    }
  }
  # return whatever most updated model is at the end, even if it hasn't converged... 
  # add convergence metric if it is still failing
  if(length(model_with_warnings_obj$warnings) > 0){
    conv_check <- performance::check_convergence(model_with_warnings_obj$result)
    model_with_warnings_obj <- append(model_with_warnings_obj, attr(conv_check, "gradient"))
  }
  return(model_with_warnings_obj)
}

# for each chromosome
full_cond_results <- data.table::rbindlist(lapply(rev(unique(gene_chr_map$Chr)), function(chr){
  message("Processing chromosome ", chr)
  # set chr 
  chr_file <- ifelse(chr == 23, "X", chr)
  # read in data
  message("loading data")
  load(paste0(dir, "input_data/eqtl_files_", chr_file, ".rda"))
  message("data loaded")
  
  # set genes
  genes_on_chr <- gene_chr_map %>%
    filter(Chr == chr) %>%
    pull(Gene)
  
  message("N genes in chr:", length(genes_on_chr))
  
  # run per gene
  data.table::rbindlist(lapply(genes_on_chr, function(gene_id){
    message("  Gene: ", gene_id)
    # get significance threshold 
    threshold <- thresholds %>% filter(gene == gene_id) %>% pull(threshold)
    # get all sig snps for this gene
    sig.snps <- signals_for_conditional_mapping %>% filter(phenotype_id == gene_id)
    pairs.eqtl <- sig.snps %>% 
      select(phenotype_id, snps, chr) %>% 
      #rename("phenotype_id" = "Gene", "snps" = "SNP", "chr" = "Chr")
      rename("Gene" = "phenotype_id", "Chr" = "chr", "SNP" = "snps")
    # re-name in case not rsID
    pairs.eqtl$SNP <- gsub(":", ".", pairs.eqtl$SNP)
    message("N SNPs:", nrow(pairs.eqtl))
    
    # get number of signals in the gene to test 
    irange <- 1:nrow(pairs.eqtl)
    
    # for each signal in this gene 
    # test whether it is still significant while controlling for other signals 
    results <- data.table::rbindlist(lapply(irange, function(j){
      # set pairs to test and pairs to control 
      pair_testing <- pairs.eqtl[j,]
      other_pairs <- pairs.eqtl[!(pairs.eqtl$SNP == pair_testing$SNP), ]
      
      # get snps from other_pairs
      snps.control <- geno[, other_pairs$SNP, drop=FALSE]
      
      model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene_id, ]) ~
                                                       0 + 
                                                       snps.control +
                                                       timepoint +
                                                       covariates +
                                                       peer_factors[, 1:n.peer] +
                                                       (1|individual) + (1|batch), 
                                                     subset = complete.cases(geno[, as.character(pair_testing$SNP)]),
                                                     REML=FALSE))
      #if there was a warning in the null model 
      if(length(model.null$warnings) > 0){
        model.null <- update_model_fun(model.null)
      }
      model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene_id, ]) ~
                                                       0 + 
                                                       as.numeric(geno[, as.character(pair_testing$SNP)]) +
                                                       snps.control +
                                                       timepoint +
                                                       covariates + 
                                                       peer_factors[, 1:n.peer] +
                                                       (1|individual) + (1|batch),
                                                     subset = complete.cases(geno[, as.character(pair_testing$SNP)]),
                                                     REML=FALSE))
      #if there was a warning in the null model 
      if(length(model.test$warnings) > 0){
        model.test <- update_model_fun(model.test)
      }
      data.frame(matrix(data=c(gene_id, 
                               pair_testing$SNP,
                               chr,
                               summary(model.test$result)$coefficients[1, ],
                               anova(model.null$result, 
                                     model.test$result)$'Pr(>Chisq)'[2],
                               length(model.null$warnings) > 0, # model warnings 
                               length(model.test$warnings) > 0, # model warnings 
                               any(grepl("failed to converge", model.null$result@optinfo$conv$lme4$messages)), # convergence warning
                               any(grepl("failed to converge", model.test$result@optinfo$conv$lme4$messages)), # convergence warning
                               any(grepl("Hessian", model.null$result@optinfo$conv$lme4$messages)), # specific convergence warning - hessian 
                               any(grepl("Hessian", model.test$result@optinfo$conv$lme4$messages)), # specific convergence warning - hessian 
                               any(grepl("singular", model.null$result@optinfo$conv$lme4$messages)), # singular fit
                               any(grepl("singular", model.test$result@optinfo$conv$lme4$messages)), # singular fit
                               ifelse(any(grepl("failed to converge", model.null$result@optinfo$conv$lme4$messages)), as.numeric(model.null[[4]]), NA),
                               ifelse(any(grepl("failed to converge", model.test$result@optinfo$conv$lme4$messages)), as.numeric(model.test[[4]]), NA)), nrow=1, ncol=17))
      
    }))
    
    colnames(results) <- c("Gene", "SNP", "Chr", "eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval", "model_warnings_null", "model_warnings_test", "conv_warnings_null", "conv_warnings_test", "hess_warnings_null", "hess_warnings_test", "singular_fit_null", "singular_fit_test", "conv_val_null", "conv_val_test")
    results %>%
      mutate(eQTL_beta = as.numeric(eQTL_beta),
             eQTL_SE = as.numeric(eQTL_SE),
             eQTL_t = as.numeric(eQTL_t),
             eQTL_pval = as.numeric(eQTL_pval), 
             conditionally_independent = eQTL_pval < threshold)
  }))
  
}))

# save this 
saveRDS(full_cond_results, file = paste0(dir, "output_data/flanders_conditional_eQTL/conditional_eQTL_results.rds"))
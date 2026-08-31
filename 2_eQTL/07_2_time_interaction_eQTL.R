# 07_2_time_interaction_eQTL.R

################################################################################

# 7.2. Run eQTL interaction with time 

################################################################################

# Aim: Map eQTL interaction models with time-point 

########################### Output paths ##########################
dir <- "genotyping/analysis/eQTL/"

########################### Input paths ###########################

########################### Parameters ############################
# set number of PEER factors (same as original eQTL mapping)
n.peer <- 30

########################### Load packages ###########################
library(lme4)
library(data.table)
library(tidyverse)

########################### Analysis ###########################

# load data
load(paste0(dir, "input_data/eqtl_int_files_flanders.rda"))

# Number of SNP-gene pairs that will be tested
irange <- 1:nrow(pairs.int) # 23859 
pairs.int <- as.data.frame(pairs.int)

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


# Test each SNP-gene pair in turn for an interaction
results <- do.call(rbind, lapply(irange, function(i){
  
  # Check that there are enough patients in each group that are minor allele homs
  if ((length(unique(individual[which(geno.int[, pairs.int[i, 2]] == 2
                                      & timepoint == levels(timepoint)[1])])) > 1 &
       length(unique(individual[which(geno.int[, pairs.int[i, 2]] == 2
                                      & timepoint == levels(timepoint)[2])])) > 1) &
      (length(unique(individual[which(geno.int[, pairs.int[i, 2]] == 2
                                      & timepoint == levels(timepoint)[3])])) > 1 & 
       length(unique(individual[which(geno.int[, pairs.int[i, 2]] == 2
                                      & timepoint == levels(timepoint)[4])])) > 1)){
    
    # is there more than one signal for this gene?
    gene <- pairs.int[i, 1]
    signals <- pairs.int[pairs.int$Gene == gene, ]
    
    # If only 1 signal
    if(nrow(signals) == 1){
      # null model includes genotype but not interaction 
      model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                                                       0 + # fake intercept - real intercept = timepoint
                                                       geno.int[, as.character(pairs.int[i, 2])] + # genotype
                                                       timepoint + # timepoint
                                                       covariates + # covariates
                                                       peer_factors[, 1:n.peer] + # peer factors
                                                       (1|individual) + (1|batch), # individual and batch 
                                                     subset = complete.cases(geno.int[, as.character(pairs.int[i, 2])]),
                                                     REML=FALSE))
      model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                                                       0 + # fake intercept - real intercept = timepoint
                                                       geno.int[, as.character(pairs.int[i, 2])] + # genotype
                                                       timepoint + # timepoint
                                                       covariates + # covariates
                                                       peer_factors[, 1:n.peer] + # peer factors
                                                       geno.int[, as.character(pairs.int[i, 2])]*timepoint + # interaction term 
                                                       (1|individual) + (1|batch), # individual and batch 
                                                     subset = complete.cases(geno.int[, as.character(pairs.int[i, 2])]),
                                                     REML=FALSE))
    } else{ 
      # If there are multiple signals for the gene include other SNPs in models
      other.snps <- signals$SNP[signals$SNP != pairs.int[i, 2]]
      model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                                                       0 + # fake intercept - real intercept = timepoint
                                                       geno.int[, as.character(pairs.int[i, 2])] + # genotype
                                                       timepoint + # timepoint
                                                       covariates + # covariates
                                                       peer_factors[, 1:n.peer] + # peer factors
                                                       (1|individual) + (1|batch) +  # individual and batch 
                                                       geno.int[, as.character(other.snps)], # other snps,
                                                     subset = complete.cases(geno.int[, as.character(pairs.int[i, 2])]),
                                                     REML=FALSE))
      model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~ # gene expression
                                                       0 + # fake intercept - real intercept = timepoint
                                                       geno.int[, as.character(pairs.int[i, 2])] + # genotype
                                                       timepoint + # timepoint
                                                       covariates + # covariates
                                                       peer_factors[, 1:n.peer] + # peer factors
                                                       geno.int[, as.character(pairs.int[i, 2])]*timepoint + # interaction term 
                                                       (1|individual) + (1|batch) +  # individual and batch 
                                                       geno.int[, as.character(other.snps)], # other snps
                                                     subset = complete.cases(geno.int[, as.character(pairs.int[i, 2])]),
                                                     REML=FALSE))
      
    }
    # look at the n minor alleles in each group that actually ended up in each model 
    indx_incl <- as.integer(rownames(model.frame(model.test$result)))
    geno_incl <- geno.int[indx_incl,as.character(pairs.int[i, 2])]
    tp_incl <- timepoint[indx_incl]
    ind_incl <- individual[indx_incl]
    tp1_nminoralleles <- length(unique(ind_incl[which(geno_incl == 2
                                                      & tp_incl == levels(timepoint)[1])]))
    tp2_nminoralleles <- length(unique(ind_incl[which(geno_incl == 2
                                                      & tp_incl == levels(timepoint)[2])]))
    tp3_nminoralleles <- length(unique(ind_incl[which(geno_incl == 2
                                                      & tp_incl == levels(timepoint)[3])]))
    tp4_nminoralleles <- length(unique(ind_incl[which(geno_incl == 2
                                                      & tp_incl == levels(timepoint)[4])]))
    
    c(summary(model.test$result)$coefficients[1,], # get results info for the genotype
      summary(model.test$result)$coefficients["timepoint12_weeks",], # get results for 12 week time-point (the different intercepts)
      summary(model.test$result)$coefficients["timepoint20_weeks",], # get results for 20 week time-point (the different intercepts)
      summary(model.test$result)$coefficients["timepoint28_weeks",], # get results for 28 week time-point (the different intercepts)
      summary(model.test$result)$coefficients["timepoint36_weeks",], # get results for 36 week time-point (the different intercepts)
      summary(model.test$result)$coefficients[dim(summary(model.test$result)$coefficients)[1]-2,], # get interaction results for 20 weeks
      summary(model.test$result)$coefficients[dim(summary(model.test$result)$coefficients)[1]-1,], # get interaction results for 28 weeks
      summary(model.test$result)$coefficients[dim(summary(model.test$result)$coefficients)[1],], # get interaction results for 36 weeks
      anova(model.null$result, model.test$result)$'Pr(>Chisq)'[2], # pval of interaction model 
      length(model.null$warnings) > 0, # model warnings 
      length(model.test$warnings) > 0, # model warnings 
      any(grepl("failed to converge", model.null$result@optinfo$conv$lme4$messages)), # convergence warning
      any(grepl("failed to converge", model.test$result@optinfo$conv$lme4$messages)), # convergence warning
      any(grepl("Hessian", model.null$result@optinfo$conv$lme4$messages)), # specific convergence warning - hessian 
      any(grepl("Hessian", model.test$result@optinfo$conv$lme4$messages)), # specific convergence warning - hessian 
      any(grepl("singular", model.null$result@optinfo$conv$lme4$messages)), # singular fit
      any(grepl("singular", model.test$result@optinfo$conv$lme4$messages)), # singular fit
      tp1_nminoralleles, # tp 1 n minor alleles in the final fit model 
      tp2_nminoralleles, # tp 2 n minor alleles in the final fit model 
      tp3_nminoralleles, # tp 3 n minor alleles in the final fit model 
      tp4_nminoralleles, # tp 4 n minor alleles in the final fit model 
      "tested")
    
  } else {
    # if there are not enouch patients in each group that are minor allele homs -> return NA
    c(rep(NA, 37), "untested_minor_allele_homs")
  }
  
  
}))

colnames(results) <- c("eQTL_beta_baseline",
                       "eQTL_SE_baseline",
                       "eQTL_t_baseline",
                       
                       "expr_mean_wk12",
                       "expr_se_wk12",
                       "expr_t_wk12",
                       
                       "expr_mean_wk20",
                       "expr_se_wk20",
                       "expr_t_wk20",
                       
                       "expr_mean_wk28",
                       "expr_se_wk28",
                       "expr_t_wk28",
                       
                       "expr_mean_wk36",
                       "expr_se_wk36",
                       "expr_t_wk36",
                       
                       "eqtl_diff_wk20_interaction",
                       "eqtl_diff_SE_wk20_interaction",
                       "eqtl_diff_t_wk20_interaction",
                       
                       "eqtl_diff_wk28_interaction",
                       "eqtl_diff_SE_wk28_interaction",
                       "eqtl_diff_t_wk28_interaction",
                       
                       "eqtl_diff_wk36_interaction",
                       "eqtl_diff_SE_wk36_interaction",
                       "eqtl_diff_t_wk36_interaction",
                       
                       "Interaction_pval",
                       "model_warnings_null", 
                       "model_warnings_test", 
                       "conv_warnings_null", 
                       "conv_warnings_test", 
                       "hess_warnings_null", 
                       "hess_warnings_test", 
                       "singular_fit_null", 
                       "singular_fit_test", 
                       
                       "tp1_nminoralleles", 
                       "tp2_nminoralleles", 
                       "tp3_nminoralleles", 
                       "tp4_nminoralleles",
                       
                       "reason_not_tested")

results <- data.frame(Gene = pairs.int[irange, 1],
                      SNP = pairs.int[irange, 2], 
                      results)

saveRDS(results, paste0(dir, "output_data/interactions_flanders_2/tp_int_results_flanders_2.rds"))
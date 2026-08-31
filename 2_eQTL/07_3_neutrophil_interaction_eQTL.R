# 07_3_neutrophil_interaction_eQTL.R

################################################################################

# 7.3. Run eQTL interaction with neutrophils 

################################################################################

# Aim: Identify eQTL with effects varying by neutrophil proportion 

########################### Paths ##########################
dir <- "genotyping/analysis/eQTL/"

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
  
  # Check that there are enough patients in each group (above and below 0) that are minor allele homs
  # In the case of cell proportions, an individual could appear in both groups 
  if (length(unique(individual[which(geno.int[, pairs.int[i, 2]] == 2
                                     & covariates[,1] > 0)])) > 1 &
      length(unique(individual[which(geno.int[, pairs.int[i, 2]] == 2
                                     & covariates[,1] < 0)])) > 1){
    
    # check that the main effect still remains with the inclusion of outcome (outcomes only - this is cell prop)
    
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
                                                       geno.int[, as.character(pairs.int[i, 2])]*covariates[,1] + # interaction term 
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
                                                       geno.int[, as.character(pairs.int[i, 2])]*covariates[,1] + # interaction term  
                                                       (1|individual) + (1|batch) +  # individual and batch 
                                                       geno.int[, as.character(other.snps)], # other snps
                                                     subset = complete.cases(geno.int[, as.character(pairs.int[i, 2])]),
                                                     REML=FALSE))
      
    }
    # look at the n minor alleles in each group that actually ended up in each model 
    indx_incl <- as.integer(rownames(model.frame(model.test$result)))
    geno_incl <- geno.int[indx_incl,as.character(pairs.int[i, 2])]
    neut_incl <- covariates[indx_incl,1]
    ind_incl <- individual[indx_incl]
    bottom_half_nminoralleles <- length(unique(ind_incl[which(geno_incl == 2
                                                              & neut_incl < 0)]))
    top_half_nminoralleles <- length(unique(ind_incl[which(geno_incl == 2
                                                           & neut_incl > 0)]))
    
    c(summary(model.test$result)$coefficients[1,], # get results info for the genotype
      summary(model.test$result)$coefficients["covariatesNeutrophils_centered",], # get results for neutrophils
      summary(model.test$result)$coefficients["geno.int[, as.character(pairs.int[i, 2])]:covariates[, 1]",], # get interaction results for neutrophils
      #summary(model.test$result)$coefficients[dim(summary(model.test$result)$coefficients)[1],], # get interaction results for neutrophils
      anova(model.null$result, model.test$result)$'Pr(>Chisq)'[2], # pval of interaction model 
      length(model.null$warnings) > 0, # model warnings 
      length(model.test$warnings) > 0, # model warnings 
      any(grepl("failed to converge", model.null$result@optinfo$conv$lme4$messages)), # convergence warning
      any(grepl("failed to converge", model.test$result@optinfo$conv$lme4$messages)), # convergence warning
      any(grepl("Hessian", model.null$result@optinfo$conv$lme4$messages)), # specific convergence warning - hessian 
      any(grepl("Hessian", model.test$result@optinfo$conv$lme4$messages)), # specific convergence warning - hessian 
      any(grepl("singular", model.null$result@optinfo$conv$lme4$messages)), # singular fit
      any(grepl("singular", model.test$result@optinfo$conv$lme4$messages)), # singular fit
      bottom_half_nminoralleles, # bottom half neutrophil distribution n minor alleles in the final fit model 
      top_half_nminoralleles, # top half neutrophil distribution n minor alleles in the final fit model 
      "tested")
    
  } else {
    # if there are not enouch patients in each group that are minor allele homs -> return NA
    c(rep(NA, 20), "untested_minor_allele_homs")
  }
  
  
}))

colnames(results) <- c("eqtl_beta_mean_neut", "eqtl_SE_mean_neut", "eqtl_t_mean_neut",
                       "neut_effect_expr", "neut_effect_expr_SE", "neut_eff_expr_t",
                       "eqtl_neut_interaction", "eqtl_neut_interaction_se", "eqtl_neut_interaction_t",
                       "Interaction_pval",
                       "model_warnings_null", 
                       "model_warnings_test", 
                       "conv_warnings_null", 
                       "conv_warnings_test", 
                       "hess_warnings_null", 
                       "hess_warnings_test", 
                       "singular_fit_null", 
                       "singular_fit_test", 
                       "bottom_half_nminoralleles",
                       "top_half_nminoralleles",
                       "reason_not_tested")

results <- data.frame(Gene = pairs.int[irange, 1],
                      SNP = pairs.int[irange, 2], 
                      results)

saveRDS(results, paste0(dir, "output_data/interactions_flanders_2/neutrophil_int_results_flanders_2.rds"))

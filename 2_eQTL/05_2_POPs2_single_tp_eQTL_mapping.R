# 05_2_POPs2_single_tp_eQTL_mapping.R

################################################################################

# 5.2. cis-eQTL mapping at a single tp: for one gene identified by job index in POPs2 TP

################################################################################

# Aim: Map eQTL in POPs2 at each time-point separately 
# Note: must be run across each time-point by varuing the "timepoint_in" input parameter

########################### Parameters ############################
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
#args <- c("1", "10", "12_weeks")
# first argument from the job array is the gene for which to map eQTL
jobindex <- as.numeric(args[1])
# second argument is the number of PEER factors to include in the model
n.peer <- as.numeric(args[2])
# third argument is the timepoint
timepoint_in <- as.character(args[3])
print("args loaded")
print(jobindex)
print(n.peer)
print(timepoint_in)

########################### Output paths ##########################
out_dir <- paste0("genotyping/analysis/eQTL/output_data/single_tp_analysis/", timepoint_in, "/chr") 

########################### Input paths ###########################
in_dir <- "genotyping/analysis/eQTL/input_data/single_tp_eQTL/"

########################### Load packages ###########################
library(Matrix)
library(lme4)
library(data.table)
library(performance)
library(dplyr)
print("packages loaded")

########################### Analysis ###########################

# read in the list of SNP-gene pairs
pairs <- data.table::fread(paste0(in_dir, "gene_snp_pairs_cis_rnaseq_fourRNA.txt"), stringsAsFactors = F) 

# use the job index to identify the gene for which to map eqtl
genes <- unique(pairs$Gene)
gene <- genes[jobindex]

# extract the list of SNPs to test for this gene
pairs.eqtl <- pairs[which(pairs$Gene == gene), ]
pairs.eqtl$SNP <- gsub(":", ".", pairs.eqtl$SNP)
irange <- 1:nrow(pairs.eqtl)
chr <- pairs.eqtl[1, 3]
if(chr == 23){
  chr <- "X"
}

print("snps extracted")
# read in the data
load(paste(in_dir, "eqtl_files_", chr, ".rda", sep=""))
print("data loaded")

# make a list of sample indices from the correct indiduals at the correct tp to keep
samples_idx <- which(timepoint == timepoint_in)

# subset consistently across all objects used downstream
gene_exp_sub   <- gene_exp[, samples_idx, drop = FALSE]
geno_sub       <- geno[samples_idx, , drop = FALSE]
covariates_sub_2 <- covariates_sub[samples_idx, , drop = FALSE]
peer_factors_sub_2 <- as.matrix(peer_factors_sub[samples_idx, , drop = FALSE])
batch_sub      <- batch[samples_idx]

# make function to capture warnings and output 
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

results <- data.table::rbindlist(lapply(irange, function(i){
  # for each SNP, fit a LMM for gene expression including genotyping and compare
  # to a null model without genotype
  if(n.peer == 0){
    model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp_sub[gene, ]) ~
                                                     covariates_sub_2 +
                                                     (1|batch_sub), 
                                                   subset = complete.cases(geno_sub[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
    #if there was a warning in the null model 
    if(length(model.null$warnings) > 0){
      model.null <- update_model_fun(model.null)
    }
    
    
    model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp_sub[gene, ]) ~
                                                     as.numeric(geno_sub[, as.character(pairs.eqtl[i, 2])]) +
                                                     covariates_sub_2 + 
                                                     (1|batch_sub),
                                                   subset = complete.cases(geno_sub[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
    #if there was a warning in the null model 
    if(length(model.test$warnings) > 0){
      model.test <- update_model_fun(model.test)
    }
  } else {
    model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp_sub[gene, ]) ~
                                                     covariates_sub_2 +
                                                     peer_factors_sub_2[, 1:n.peer] +
                                                     (1|batch_sub), 
                                                   subset = complete.cases(geno_sub[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
    #if there was a warning in the null model 
    if(length(model.null$warnings) > 0){
      model.null <- update_model_fun(model.null)
    }
    
    
    model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp_sub[gene, ]) ~
                                                     as.numeric(geno_sub[, as.character(pairs.eqtl[i, 2])]) +
                                                     covariates_sub_2 + 
                                                     peer_factors_sub_2[, 1:n.peer] +
                                                     (1|batch_sub),
                                                   subset = complete.cases(geno_sub[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
    #if there was a warning in the null model 
    if(length(model.test$warnings) > 0){
      model.test <- update_model_fun(model.test)
    }
  }
  
  # return the information 
  data.frame(matrix(data=c(summary(model.test$result)$coefficients[2, ],
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
                           ifelse(any(grepl("failed to converge", model.test$result@optinfo$conv$lme4$messages)), as.numeric(model.test[[4]]), NA)),
                    nrow = 1, ncol=14))
  
}))

# Make results table for this gene
colnames(results) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval", "model_warnings_null", "model_warnings_test", "conv_warnings_null", "conv_warnings_test", "hess_warnings_null", "hess_warnings_test", "singular_fit_null", "singular_fit_test", "conv_val_null", "conv_val_test")

results <- data.frame(Gene = gene, SNP = pairs.eqtl[irange, 2], results)

output_dir <- paste0(out_dir, chr)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
saveRDS(results, paste0(output_dir, "/", gene, ".rds"))
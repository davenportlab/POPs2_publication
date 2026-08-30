# 04_4_interval_eQTL_mapping.R

################################################################################

# 4.4. cis-eQTL mapping in interval: for one gene identified by job index

################################################################################

# Aim: Map eQTL in interval 

########################### Parameters ############################
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
# first argument from the job array is the gene for which to map eQTL
jobindex <- as.numeric(args[1])
# second argument is the number of PEER factors to include in the model
n.peer <- as.numeric(args[2])
# subdir for peer test 
subdir <- as.character(args[3])
print("args loaded")
print(jobindex)
print(n.peer)
print(subdir)

########################### Output paths ##########################
out_dir <- paste0("interval_rna/POPS2_comparison/", subdir, "/chr") 

########################### Input paths ###########################
in_dir <- "interval_rna/POPS2_comparison/eQTL_input_data/"

########################### Load packages ###########################
library(Matrix)
library(lme4)
library(data.table)
library(performance)
library(dplyr)
print("packages loaded")

########################### Load data and Analysis ###########################

# read in the list of SNP-gene pairs
pairs <- data.table::fread(paste0(in_dir, "gene_snp_pairs_cis_rnaseq_update_10_27_25.txt"), stringsAsFactors = F) 
print("pairs loaded")

# use the job index to identify the gene for which to map eqtl
genes <- unique(pairs$Gene)
gene <- genes[jobindex]

# extract the list of SNPs to test for this gene
pairs.eqtl <- pairs[which(pairs$Gene == gene), ]
pairs.eqtl$SNP <- gsub(":", ".", pairs.eqtl$SNP) # doesn't do anything in the interval format
chr <- pairs.eqtl[1, 3]
chr_num <- chr
if(chr == 23){
  chr <- "X"
}
print("snps extracted")

# read in the data
load(paste(in_dir, "eqtl_files_", chr, ".rda", sep=""))
print("data loaded")

irange <- 1:nrow(pairs.eqtl)

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
  # for each SNP, fit a LMM for gene expression including genotyping and compare to null model 
  if(n.peer == 0){
    model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~
                                                     covariates +
                                                     (1|batch), 
                                                   subset = complete.cases(geno[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
  } else{
    model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~
                                                     covariates +
                                                     peer_factors[, 1:n.peer] +
                                                     (1|batch), 
                                                   subset = complete.cases(geno[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
  }
  #if there was a warning in the null model 
  if(length(model.null$warnings) > 0){
    model.null <- update_model_fun(model.null)
  }
  if(n.peer == 0){
    model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~
                                                     as.numeric(geno[, as.character(pairs.eqtl[i, 2])]) +
                                                     covariates + 
                                                     (1|batch),
                                                   subset = complete.cases(geno[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
  } else {
    model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~
                                                     as.numeric(geno[, as.character(pairs.eqtl[i, 2])]) +
                                                     covariates + 
                                                     peer_factors[, 1:n.peer] +
                                                     (1|batch),
                                                   subset = complete.cases(geno[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
  }
  
  
  #if there was a warning in the null model 
  if(length(model.test$warnings) > 0){
    model.test <- update_model_fun(model.test)
  }
  
  # return the information 
  data.frame(matrix(data=c(summary(model.test$result)$coefficients[2, ], # genotype effect
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
print("models run")
# Make results table for this gene
colnames(results) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval", "model_warnings_null", "model_warnings_test", "conv_warnings_null", "conv_warnings_test", "hess_warnings_null", "hess_warnings_test", "singular_fit_null", "singular_fit_test", "conv_val_null", "conv_val_test")

results <- data.frame(Gene = gene, SNP = pairs.eqtl[irange, 2], results)
print("df prepared")
output_dir <- paste0(out_dir, chr)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
saveRDS(results, paste0(output_dir, "/", gene, ".rds"))
# 01_3_ciseqtl_mapping.R

################################################################################

# 1.3. cis-eQTL mapping: for one gene identified by job index

################################################################################

# Aim: Map cis-eQTL
# Note: This script is used first for varying the n peer factors (0-50 at increments of 5) on Chromosome 1 only. 
#       Then, the script is used for mapping main effects eQTL (where n peer factors = 30)

########################### Parameters ############################
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
# first argument from the job array is the gene for which to map eQTL
jobindex <- as.numeric(args[1])
# second argument is the number of PEER factors to include in the model
n.peer <- as.numeric(args[2])
# third argument is the subdir for this run of eQTL mapping 
subdir <- as.character(args[3])

########################### Load packages ###########################
suppressWarnings(library(Matrix))
suppressWarnings(library(lme4))
suppressWarnings(library(data.table))
suppressWarnings(library(performance))

########################### Output paths ##########################
out_dir <- paste0("genotyping/analysis/eQTL/output_data/main_effects/", subdir, "/chr") 

########################### Input paths ###########################
in_dir <- "genotyping/analysis/eQTL/input_data/"

########################### Load data and Analysis ###########################


# read in the list of SNP-gene pairs
pairs <- data.table::fread(paste0(in_dir, "gene_snp_pairs_cis_rnaseq.txt"), stringsAsFactors = F) 

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

# check if results file exists (in case of re-running failed jobs)
done <- file.exists(paste0(out_dir, chr, "/", gene, ".rds")) 
if(done == FALSE){ 
  # read in the data
  load(paste(in_dir, "eqtl_files_", chr, ".rda", sep=""))
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
  #print("data loaded")
  results <- data.table::rbindlist(lapply(irange, function(i){
    # for each SNP, fit a LMM for gene expression including genotyping and compare
    # to a null model (both including SRS, diagnosis, cell proportions, genotyping PCs and PEER factors)
    
    model.null <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~
                                                     0 + 
                                                     timepoint +
                                                     covariates +
                                                     peer_factors[, 1:n.peer] +
                                                     (1|individual) + (1|batch), 
                                                   subset = complete.cases(geno[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
    #if there was a warning in the null model 
    if(length(model.null$warnings) > 0){
      model.null <- update_model_fun(model.null)
    }
    
    
    model.test <- capture_warnings_and_output(lmer(as.numeric(gene_exp[gene, ]) ~
                                                     0 + 
                                                     as.numeric(geno[, as.character(pairs.eqtl[i, 2])]) +
                                                     timepoint +
                                                     covariates + 
                                                     peer_factors[, 1:n.peer] +
                                                     (1|individual) + (1|batch),
                                                   subset = complete.cases(geno[, as.character(pairs.eqtl[i, 2])]),
                                                   REML=FALSE))
    #if there was a warning in the null model 
    if(length(model.test$warnings) > 0){
      model.test <- update_model_fun(model.test)
    }
    
    # return the information 
    data.frame(matrix(data=c(summary(model.test$result)$coefficients[1, ],
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
}
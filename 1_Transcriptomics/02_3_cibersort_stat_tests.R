# 02_3_cibersort_stat_tests.R 

################################################################################

# 2.3. Statistical tests of imputed cell proportions

################################################################################

# Aim: Test for differences in cibersort output across time-points

########################### Output paths ##########################
# bonferoni adjusted p-values for cell proportion changes across time points 
cell_prop_change_bonf_pvals_outpath = "rna-seq/analysis/cibersort/outputs/cell_prop_change_timepoints_bonf_pvals.csv"
# bonferoni adjusted p-values and cell proportion mean differences for cell proportion changes across time points
cell_prop_change_all_results_outpath = "rna-seq/analysis/cibersort/outputs/cell_prop_change_timepoints_all_results.csv"
# bonferoni adjusted p-values for cell proportion changes across seasons
season_change_cell_prop_pvals_outpath = "rna-seq/analysis/cibersort/outputs/cell_prop_change_seasons_bonf_pvals.csv"

########################### Input paths ###########################
# cibersort output 
cibersort_inpath <- "rna-seq/analysis/cibersort/outputs/CIBERSORTx_Adjusted.txt"
# sample info path 
sample_info_inpath <- "rna-seq/data/sample_covariates_clin_tech.csv"

########################### Load packages ###########################
library(tidyverse)

########################### Load data ###########################
# cibersort output imputed cell proportions
cell_prop <- read.delim(cibersort_inpath)
# sample info 
sample.info <- read.csv(sample_info_inpath) 

########################### Analysis ###########################
# Join Metadata patients and time points with the cell proportions
cell_prop <- cell_prop %>% left_join(sample.info,  by=c("Mixture" ="RNA_sanger_sample_id"))

# Test for differences in cell type proportion between time points 
# Stats refresher for paired t test - https://cambiotraining.github.io/corestats/materials/cs1_practical_two-samples-paired.html

# We fail the assumption of normality for the two-sample Student’s t-test 
# Instead, we use Wilcoxon rank sum test 

# Due to missingness, for each pair of sample time points, 
#   subset to only individuals with a sample at BOTH time points
#   test whether differences in that sample subset != 0

# Make lists of individuals for each time point comparison
# Filter ANON_ID = OBGYN0496 and OBGYN0309 because they have two twelve week time points

ind_lists <- list(
  # list individuals with both 12 and 20 week sample
  individuals_12v20 <- (cell_prop %>% filter(Sample_taken_at %in% c("12_weeks", "20_weeks")) %>% group_by(ANON_ID) %>% summarize(n = n()) %>% filter(n == 2, !ANON_ID %in% c("OBGYN0496", "OBGYN0309")))$ANON_ID, 
  # list individuals with both 12 and 28 week sample
  individuals_12v28 <- (cell_prop %>% filter(Sample_taken_at %in% c("12_weeks", "28_weeks")) %>% group_by(ANON_ID) %>% summarize(n = n()) %>% filter(n == 2, !ANON_ID %in% c("OBGYN0496", "OBGYN0309")))$ANON_ID, 
  # list individuals with both 12 and 36 week sample
  individuals_12v36 <- (cell_prop %>% filter(Sample_taken_at %in% c("12_weeks", "36_weeks")) %>% group_by(ANON_ID) %>% summarize(n = n()) %>% filter(n == 2, !ANON_ID %in% c("OBGYN0496", "OBGYN0309")))$ANON_ID, 
  # list individuals with both 20 and 28 week sample
  individuals_20v28 <- (cell_prop %>% filter(Sample_taken_at %in% c("20_weeks", "28_weeks")) %>% group_by(ANON_ID) %>% summarize(n = n()) %>% filter(n == 2))$ANON_ID, 
  # list individuals with both 20 and 36 week sample
  individuals_20v36 <- (cell_prop %>% filter(Sample_taken_at %in% c("20_weeks", "36_weeks")) %>% group_by(ANON_ID) %>% summarize(n = n()) %>% filter(n == 2))$ANON_ID, 
  # list individuals with both 28 and 36 week sample
  individuals_28v36 <- (cell_prop %>% filter(Sample_taken_at %in% c("28_weeks", "36_weeks")) %>% group_by(ANON_ID) %>% summarize(n = n()) %>% filter(n == 2))$ANON_ID
)

# list time point comparisons 
tp_comps <- list(`12v20` = c("12_weeks", "20_weeks"),
                 `12v28` = c("12_weeks", "28_weeks"),
                 `12v36` = c("12_weeks", "36_weeks"),
                 `20v28` = c("20_weeks", "28_weeks"),
                 `20v36` = c("20_weeks", "36_weeks"), 
                 `28v36` = c("28_weeks", "36_weeks"))
# list cell types 
cell_types <- colnames(cell_prop)[2:23]

# prep columns to save results 
cell_prop_changes_pval <- c("cell type", "12v20", "12v28", "12v36", "20v28", "20v36", "28v36")
cell_prop_mean_diff <- cell_prop_changes_pval

# for each cell type
for(i in 1:length(cell_types)){
  # test for differences between all possible time point combinations
  # cell type list p-vals
  cell_type_pvals <- c(cell_types[i])
  # cell type list mean_diffs
  cell_type_mean_diffs <- c(cell_types[i])
  
  # for each time point pair  
  for(j in 1:length(tp_comps)){
    # list which time point comparison is being tested
    tps_test = unname(tp_comps[j])[[1]]
    # subset data to samples from individuals with both time points 
    subset_for_test <- cell_prop %>% 
      filter(Sample_taken_at %in% tps_test, ANON_ID %in% ind_lists[[j]]) %>% dplyr::select(cell_types[i], ANON_ID, Sample_taken_at)
    # put this in the right format for the wilcox test 
    subset_for_test <- subset_for_test %>% pivot_wider(id_cols = ANON_ID, names_from = Sample_taken_at, values_from = cell_types[i]) %>% 
      # change for each iteration!
      mutate(prop_change = cur_data()[[2]] - cur_data()[[3]])
    
    # run wilcox test if the variance is not 0 for that cell type 
    if(var(subset_for_test$prop_change) != 0 ){
      p_val <- (wilcox.test(subset_for_test$prop_change, mu = 0, alternative = "two.sided"))$p.value
      # add this to the list
      cell_type_pvals <- c(cell_type_pvals, p_val)
      # compute cell type mean diff
      cell_type_mean_diffs <- c(cell_type_mean_diffs, mean(subset_for_test$prop_change))
      # if the variance is 0, the test doesn't work, so just add NA
    }else{
      cell_type_pvals <- c(cell_type_pvals, NA)
      # compute cell type mean diff
      cell_type_mean_diffs <- c(cell_type_mean_diffs, mean(subset_for_test$prop_change))
    }
  }
  # put together all values from this cell type across time point comps 
  cell_prop_changes_pval <- rbind(cell_prop_changes_pval, cell_type_pvals)
  cell_prop_mean_diff <- rbind(cell_prop_mean_diff, cell_type_mean_diffs)
}

# format the results of the wilcox test 
cell_prop_changes_pval <- as.data.frame(cell_prop_changes_pval)
colnames(cell_prop_changes_pval) <- cell_prop_changes_pval[1,]
cell_prop_changes_pval <- cell_prop_changes_pval[-1,]
rownames(cell_prop_changes_pval) <- cell_prop_changes_pval$`cell type`
cell_prop_changes_pval <- cell_prop_changes_pval[,-1]
cell_prop_changes_pval <- cell_prop_changes_pval %>% mutate_all(., as.numeric)

# Correct for multiple testing burden using p.adjust
# do cell prop multiple testing correction
cell_prop_changes_adj_pval <- cell_prop_changes_pval %>% as.matrix() %>% as.vector() %>% p.adjust(method = "bonferroni") %>% matrix(nrow = 22) %>% as.data.frame()
# format 
colnames(cell_prop_changes_adj_pval) <- colnames(cell_prop_changes_pval)
rownames(cell_prop_changes_adj_pval) <- rownames(cell_prop_changes_pval)

# write this out 
write.csv(cell_prop_changes_adj_pval, file = cell_prop_change_bonf_pvals_outpath, quote = FALSE)

# format mean differences 
cell_prop_mean_diff <- as.data.frame(cell_prop_mean_diff)
colnames(cell_prop_mean_diff) <- cell_prop_mean_diff[1,]
cell_prop_mean_diff <- cell_prop_mean_diff[-1,]
rownames(cell_prop_mean_diff) <- cell_prop_mean_diff$`cell type`
cell_prop_mean_diff <- cell_prop_mean_diff[,-1]
cell_prop_mean_diff <- cell_prop_mean_diff %>% mutate_all(., as.numeric)

# Join cell prop differences and p-vals for large output df 
# start with p values 
cell_prop_changes_info <- cell_prop_changes_adj_pval %>% 
  # make rowname a column for cell type 
  rownames_to_column(var = "cell_type") %>% 
  # pivot longer
  pivot_longer(-cell_type, names_to = "time_point_comp", values_to = "bonf_adj_p_val") %>%
  # join this to mean diff in cell props per comp
  left_join((cell_prop_mean_diff %>% 
               # make rowname a column for cell type
               rownames_to_column(var = "cell_type") %>% 
               # pivot longer 
               pivot_longer(-cell_type, names_to = "time_point_comp", values_to = "cell_prop_mean_diff")), 
            # join by cel type and time point comp 
            by = c("cell_type", "time_point_comp")) %>%
  # add column about significance at 0.05
  mutate(sig_at_0.05 = bonf_adj_p_val < 0.05)

# write out this df
cell_prop_changes_info %>%
  write.csv(file = cell_prop_change_all_results_outpath, quote = FALSE, row.names = FALSE)

# Check for differences in distribution of cell prop at each time point by season 
# Use Kruskal-Wallis test
season_changes_pval <- c("cell type", "12_weeks", "20_weeks", "28_weeks", "36_weeks")
# list time points
time_points <- unique(cell_prop$Sample_taken_at)

# for each cell type
for(i in 1:length(cell_types)){
  # make list for p-vals
  cell_type_pvals <- c(cell_types[i])
  # test for differences between all possible time point combinations
  for(j in 1:length(time_points)){
    # make data subset
    dat_sub <- cell_prop %>% filter(Sample_taken_at == time_points[j]) %>% select(cell_types[i], Season)
    # perform kruskal wallace test
    p_val <- (kruskal.test(get(cell_types[i]) ~ Season,
                           data = dat_sub)[[3]])
    cell_type_pvals <- c(cell_type_pvals, p_val)
  }  
  season_changes_pval <- rbind(season_changes_pval, cell_type_pvals)
}

# format seasonal cell prop changes significance
season_changes_pval_df <- as.data.frame(season_changes_pval)
colnames(season_changes_pval_df) <- season_changes_pval_df[1,]
season_changes_pval_df <- season_changes_pval_df[-1,]
rownames(season_changes_pval_df) <- season_changes_pval_df[,1]
season_changes_pval_df <- season_changes_pval_df[,-1]
season_changes_pval_df <- season_changes_pval_df %>% mutate_all(., as.numeric)

# Correct for multiple testing burden - adjust using p.adjust
season_changes_adj_pval <- season_changes_pval_df %>% as.matrix() %>% as.vector() %>% p.adjust(method = "bonferroni") %>% matrix(nrow = 22) %>% as.data.frame()
# format 
colnames(season_changes_adj_pval) <- colnames(season_changes_pval_df)
rownames(season_changes_adj_pval) <- rownames(season_changes_pval_df)
# write this out
write.csv(season_changes_adj_pval, file = season_change_cell_prop_pvals_outpath, quote = FALSE)
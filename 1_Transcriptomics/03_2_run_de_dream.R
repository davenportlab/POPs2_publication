# 03_2_run_de_dream.R 

################################################################################

# 3.2. Run differential expression analysis 

################################################################################

# Aim: run dream to identify differentially expressed genes 

########################### Output paths ##########################
# fit
out_obj_path <- "rna-seq/analysis/dream/outputs/vobjDream_NOform"
# top table
toptable_outpath <- "rna-seq/analysis/dream/outputs/pops_toptable.rds"

########################### Input paths ###########################
# voom object path (made in RNA-seq_preprocessing.R)
voom_obj_inpath <- "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
# sample info input 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Parameters ############################
# specify dream formula
form <- "~0+Sample_taken_at+(1|ANON_ID)+(1|id_run_position)"
# specify variance partition formula
variance_partition_formula <- "~(1|Sample_taken_at)+(1|ANON_ID)+(1|id_run_position)"

########################### Load packages ###########################
print("loading packages")
library("edgeR")
library("BiocParallel")
library("variancePartition")
suppressPackageStartupMessages(library(tidyverse))
library(stats)
library(lme4)
print("packages loaded")

########################### Load data ###########################
print("Loading data...")
# load voom object
vobjDream <- readRDS(voom_obj_inpath)
# read in metadata that shows which samples we will use 
sample.info <- read.csv(sample_info_inpath)
# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
print("Data loaded.")
Sys.time()

########################### Analysis ###########################

set.seed(1)

print(c("starting at"))
Sys.time()

print(paste("Voom object:", voom_obj_inpath))
print(paste("Formula for dream:", form))
print(paste("Variance partition formula:", variance_partition_formula))
print(paste("Output object:", out_obj_path))

# 1. load data 
sample.info <- sample.info %>%
  # order by count matrix 
  arrange(factor(RNA_sanger_sample_id, levels = colnames(vobjDream$E)))

# add covariate 
rownames(sample.info) <- sample.info$RNA_sanger_sample_id
# check we have 3779 samples and 60678 genes
dim(sample.info)
dim(vobjDream$E)
# check things are in correct order
all(rownames(sample.info) == colnames(vobjDream$E))


# 2. make contrasts 
L <- makeContrastsDream(form, sample.info,
                        contrasts = c(
                          compare12_20 = "Sample_taken_at20_weeks - Sample_taken_at12_weeks",
                          compare12_28 = "Sample_taken_at28_weeks - Sample_taken_at12_weeks",
                          compare12_36 = "Sample_taken_at36_weeks - Sample_taken_at12_weeks",
                          compare20_28 = "Sample_taken_at28_weeks - Sample_taken_at20_weeks",
                          compare20_36 = "Sample_taken_at36_weeks - Sample_taken_at20_weeks", 
                          compare28_36 = "Sample_taken_at36_weeks - Sample_taken_at28_weeks"
                        )
)

# 3. fit dream model with contrasts
print("start fit dream model with contrasts")
Sys.time()
fit <- dream(vobjDream, form, sample.info, L)
fit <- eBayes(fit)
print("end fit dream model with contrasts")
Sys.time()
# get names of available coefficients and contrasts for testing
colnames(fit)
print("Finished dream model")
Sys.time()

#write data 
print("Writing out model fit data")
saveRDS(fit, file = paste0(out_obj_path, "_dream_fit.rds"))

# 4. Model variance partition of all covariates together 
print("start variance partition")
vp <- fitExtractVarPartModel(vobjDream, variance_partition_formula, sample.info)
# If a factor doesn't contribute a lot to the variance or differences between groups, it might be better to just ignore it.
#dev.off()
#plotVarPart(sortCols(vp))
print("Finished variance partition")
Sys.time()

#write data 
print("Writing out variance partition data")
saveRDS(vp, file = paste0(out_obj_path, "_vp.rds"))

print("finished all")
Sys.time()


# Create toptable with all DE results (Extract top tables and make a long version including all contrasts)
# read in the results
fit <- readRDS(paste0(out_obj_path, "_dream_fit.rds"))
# Create gene ID name key 
# edit gtf since we added two genes that are not in the reference 
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

# function to make top table
make_top_table_long_form_df <- function(fit, n_genes){
  # extract top tables
  top.tables <- list(
    top.table.12.20 = topTable(fit, coef = "compare12_20", number = n_genes),
    top.table.12.28 = topTable(fit, coef = "compare12_28", number = n_genes),
    top.table.12.36 = topTable(fit, coef = "compare12_36", number = n_genes),
    top.table.20.28 = topTable(fit, coef = "compare20_28", number = n_genes),
    top.table.20.36 = topTable(fit, coef = "compare20_36", number = n_genes),
    top.table.28.36 = topTable(fit, coef = "compare28_36", number = n_genes))
  # make reference 
  top.table.comparisons <- c("12 weeks vs. 20 weeks", 
                             "12 weeks vs. 28 weeks", 
                             "12 weeks vs. 36 weeks", 
                             "20 weeks vs. 28 weeks", 
                             "20 weeks vs. 36 weeks", 
                             "28 weeks vs. 36 weeks")
  top.table.comparisons2 <- gsub(" weeks", "", top.table.comparisons)
  # make long form df
  top_table_df <- rbind(top.tables$top.table.12.20 %>% mutate(comp = "12v20") %>% rownames_to_column(var = "gtf.gene_id"),
                        top.tables$top.table.12.28 %>% mutate(comp = "12v28") %>% rownames_to_column(var = "gtf.gene_id"),
                        top.tables$top.table.12.36 %>% mutate(comp = "12v36") %>% rownames_to_column(var = "gtf.gene_id"),
                        top.tables$top.table.20.28 %>% mutate(comp = "20v28") %>% rownames_to_column(var = "gtf.gene_id"),
                        top.tables$top.table.20.36 %>% mutate(comp = "20v36") %>% rownames_to_column(var = "gtf.gene_id"),
                        top.tables$top.table.28.36 %>% mutate(comp = "28v36") %>% rownames_to_column(var = "gtf.gene_id")) %>%
    left_join(id_name, by = "gtf.gene_id")
  
  return(top_table_df)
  
}
# make top table 
toptable <- make_top_table_long_form_df(fit = fit, n_genes = 18826)
dim(toptable)
head(toptable)

# write out the toptable
saveRDS(toptable, toptable_outpath)

########################### Plot ##########################

# Variance partition plot 
vp <- readRDS(paste0(out_obj_path, "_vp.rds"))
plotVarPart(sortCols(vp), , col = c(rep("white", 3), "grey85")) + 
  scale_x_discrete(labels=c("ANON_ID" = "Individual", "id_run_position" = "Sequencing\nbatch", "Sample_taken_at" = "Time-point")) + 
  xlab("Covariate") + 
  #theme_minimal() +
  labs(title="Variance partition plot",
       subtitle = paste("Covariates used in dream differential gene expression model")) + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks = element_blank()) 

# Volcano plot
# define deg cutoffs
p_val = 0.01
log_FC = log2(1.5)

# make annotation to add n DEGs
n_degs <- toptable %>% 
  filter(abs(logFC) > log_FC, 
         adj.P.Val < p_val) %>% 
  group_by(comp) %>%
  summarize(n = n())

# Volcano plot (Main figure 1D)
toptable %>% 
  left_join(n_degs, by = "comp") %>%
  mutate(reg = ifelse(((logFC > log_FC) & (adj.P.Val < p_val)), "UP", 
                      ifelse(((logFC < -log_FC) & (adj.P.Val < p_val)), "DOWN", "no_change")), 
         adj.P.Val = ifelse(adj.P.Val < .Machine$double.xmin, .Machine$double.xmin, adj.P.Val),
         comp = paste0(comp, ", n DEGs = ", n)
  ) %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val), text = paste("Symbol:", gtf.gene_name), 
             color = reg, label = gtf.gene_name)) +
  scale_color_manual(name = "reg",
                     values =c("no_change"= "#999999", 
                               "UP" = "black", 
                               "DOWN" = "black")) +
  #ggrastr::geom_point_rast(shape = 19, show.legend = FALSE, size = 0.5) +
  geom_point(shape = 19, show.legend = FALSE, size = 0.5) +
  geom_hline(yintercept = -log10(p_val), linetype="longdash", colour="grey", size=0.5) +
  geom_vline(xintercept = log_FC, linetype="longdash", colour="#999999", size=0.5) +
  geom_vline(xintercept = -log_FC, linetype="longdash", colour="#999999", size=0.5) + 
  theme(legend.position="none", panel.grid = element_blank(), strip.background = element_blank()) + 
  coord_cartesian(clip = "off") + 
  xlab("Log 2 fold change") + 
  ylab("-log10(adjusted p-value)") + 
  facet_wrap(~comp)

# DEG quantification barplot (Supplementary Figure 6)
n_degs %>% ggplot(aes(x = comp, y = n)) + 
  geom_bar(stat = "identity", fill = "black") + 
  geom_text(aes(label = n), vjust = -0.5, color = "black") +
  xlab("Time-point comparisons") + 
  ylab("Number of DEGs") +
  ggtitle("Number of DEGs between time-points") + 
  ylim(c(0, 1800))

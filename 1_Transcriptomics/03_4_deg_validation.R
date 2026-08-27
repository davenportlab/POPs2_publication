# 03_4_deg_validation.R 

################################################################################

# 3.4. Validate differentially expressed genes 

################################################################################

# Aim: Re-run DE analysis on external datataset and compare with POPs2 

########################### Output paths ##########################
toptable_outpath <- "rna-seq/analysis/dream/outputs/Gomez-Lopez_2019_DEGs_filt_dream.csv"
biomart_outpath <- "rna-seq/analysis/dream/outputs/Gomez-Lopez_2019_DEGs_filt_dream_biomart.csv"

########################### Input paths ###########################
voom_obj_inpath <- "rna-seq/analysis/data_preprocessing/vobjDream_NOform_2.rds"
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"
# pops toptable
toptable_inpath <- "rna-seq/analysis/dream/outputs/pops_toptable.rds"

########################### Parameters ############################
# define deg cutoffs
p_val = 0.01
log_FC = log2(1.5)

########################### Load packages ###########################
library(GEOquery)
library(limma)
library(umap)
library(tidyverse)
library(hta20transcriptcluster.db)
library(AnnotationDbi)
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(eulerr)

########################### Load data ###########################
# load pops data
vobjDream <- readRDS(voom_obj_inpath)
# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
# pops toptable
toptable <- readRDS(toptable_inpath)


########################### Analysis ###########################

# DEG Validation 
# Using [Gomez-Lopez 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6928201/#SM1) data for validation, 
# selecting samples that were taken within 2 weeks of our time points. 
#                            
#   - 12 (10-14) n = 34
#   - 20 (18-22) n = 43
#   - 28 (26-30) n = 56
#   - 36 (34-38) n = 41
#                          
# The following code is from NCBI, GEO accession GSE121974 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121974)

# 1. load series and platform data from GEO
gset <- getGEO("GSE121974", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("01X2X0X2301X23X01XX3XXX3XX1322302XX213XX2123XX1X00",
               "01X112XX32X20X322X13X23X322XX20332XX2X313112X021XX",
               "2321303X01X033213122X2031X010X2X213X1XX2X12X0X023X",
               "X0XX03X2211XX32XX2X11XX002X220112X3XX0X12X23XX00X2",
               "23031231XX0X12223X1210X12XXXX1X301X3X2X30X303XX112",
               "X0X321X2312")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("12_weeks","20_weeks","28_weeks","36_weeks"))
levels(gs) <- groups
gset$group <- gsub("X", "", gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

# Now, filter genes
expr <- exprs(gset)
dim(expr)
#[1] 70523   174
expr[1:5, 1:5]
# GSM3452134 GSM3452135 GSM3452137 GSM3452139 GSM3452141
# 2824546_st   9.047975   9.368754  10.711046   9.008713  10.413861
# 2824549_st   9.108015   9.067654   9.618197   8.912863   9.712411
# 2824551_st   8.532852   8.946585  10.476199   8.762939   9.912970
# 2824554_st   8.566534   7.944315  10.335718   8.876704  10.162149
# 2827992_st   8.206761   8.494563   9.766822   9.038061   9.525643
feat <- fData(gset)

# check same order across both 
all(rownames(expr) == feat$ID)

# remove control probes (marked by not having a chromosome)
keep <- feat$seqname != "---"
expr <- expr[keep, ]
dim(expr)
# [1] 67528   174
feat <- feat[keep, ]

# map probe names to genes 
annot <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = rownames(expr),
  columns = c("SYMBOL", "ENSEMBL", "ENTREZID", "GENETYPE"),
  keytype = "PROBEID"
)
dim(annot)
# [1] 82041     5
# remove probes mapping to multiple genes
annot <- annot %>%
  filter(!is.na(ENSEMBL)) %>%
  group_by(PROBEID) %>%
  filter(n() == 1) %>%
  ungroup()
dim(annot)
# [1] 24324     5
# Filter to probes mapping to genes in POPs2 count matrix
annot <- annot %>%
  filter(ENSEMBL %in% rownames(vobjDream$E))
dim(annot)
# [1] 15176     5
# filter count matrix to these probes that are annotated
expr <- expr[annot$PROBEID, ]
dim(expr)
# [1] 15176   174
annot <- annot[match(rownames(expr), annot$PROBEID), ]
# Collapse probes to genes
expr_gene <- expr %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>%
  left_join(annot %>% select(PROBEID, ENSEMBL), by = "PROBEID") %>%
  select(-PROBEID) %>%
  group_by(ENSEMBL) %>%
  summarise(across(everything(), mean), .groups = "drop")
dim(expr_gene)
# [1] 12815   175
expr_gene[1:5, 1:5]
# # A tibble: 5 × 5
# ENSEMBL         GSM3452134 GSM3452135 GSM3452137 GSM3452139
# <chr>                <dbl>      <dbl>      <dbl>      <dbl>
# 1 ENSG00000000003       2.85       2.78       2.73       2.79
# 2 ENSG00000000419       6.99       7.22       7.47       7.27
# 3 ENSG00000000457       6.32       6.51       6.50       6.11
# 4 ENSG00000000460       4.33       4.25       4.35       4.36
# 5 ENSG00000000938       9.29       9.39       9.34       9.19

# Go back to expr matrix
expr_gene <- as.data.frame(expr_gene)
rownames(expr_gene) <- expr_gene$ENSEMBL
expr_gene$ENSEMBL <- NULL

expr_gene <- as.matrix(expr_gene)
dim(expr_gene)
# [1] 12815   174
expr_gene[1:5, 1:5]
# GSM3452134 GSM3452135 GSM3452137 GSM3452139 GSM3452141
# ENSG00000000003   2.846883   2.782586   2.726770   2.786026   2.691569
# ENSG00000000419   6.988888   7.219235   7.473761   7.269492   7.211432
# ENSG00000000457   6.318894   6.514190   6.501372   6.110716   6.062079
# ENSG00000000460   4.332121   4.250908   4.351790   4.359466   4.290800
# ENSG00000000938   9.290916   9.388593   9.335676   9.187807   9.343550

####  Run dream ####

# make df of metadata - rowname as sample ID, column for each variable
sample.info <- as.data.frame(cbind(Sample_taken_at = gset$group,
                                   rna_sample_id = gset$geo_accession,
                                   ANON_ID = gsub("individual: ", "", gset$characteristics_ch1),
                                   gestational_age_wk = gsub("gestational age (wk):", "", gset$characteristics_ch1.1, fixed = TRUE)))
rownames(sample.info) <- sample.info$rna_sample_id
dim(sample.info)
# [1] 174   4
head(sample.info)
# Sample_taken_at rna_sample_id ANON_ID gestational_age_wk
# GSM3452134        12_weeks    GSM3452134  MI_829               13.3
# GSM3452135        20_weeks    GSM3452135  MI_829               19.7
# GSM3452137        28_weeks    GSM3452137  MI_829               29.4
# GSM3452139        12_weeks    GSM3452139  MI_293               12.4
# GSM3452141        28_weeks    GSM3452141  MI_293               29.3
# GSM3452142        36_weeks    GSM3452142  MI_293               35.6
# make formula
form <- ~ 0 + Sample_taken_at + (1 | ANON_ID)

# make contrasts
# make contrasts
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
# make dge list object
dge <- DGEList(expr_gene)
dge <- calcNormFactors(dge) # Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.

# run dream with contrasts
fit <- dream(dge, form, sample.info, L) # 11:22 start
fit <- eBayes(fit)

# load naming data
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

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
  top_table_df <- rbind(rbind(rbind(top.tables$top.table.12.20 %>% 
                                      mutate(comp = "12v20") %>% 
                                      rownames_to_column(var = "gtf.gene_id") , 
                                    top.tables$top.table.12.28 %>% 
                                      mutate(comp = "12v28")  %>% 
                                      rownames_to_column(var = "gtf.gene_id")), 
                              rbind(top.tables$top.table.12.36 %>% 
                                      mutate(comp = "12v36") %>% 
                                      rownames_to_column(var = "gtf.gene_id"), 
                                    top.tables$top.table.20.28 %>% 
                                      mutate(comp = "20v28") %>% 
                                      rownames_to_column(var = "gtf.gene_id"))), 
                        rbind(top.tables$top.table.20.36 %>% 
                                mutate(comp = "20v36") %>% 
                                rownames_to_column(var = "gtf.gene_id"), 
                              top.tables$top.table.28.36 %>% 
                                mutate(comp = "28v36") %>% 
                                rownames_to_column(var = "gtf.gene_id"))) %>% 
    left_join(id_name, by = "gtf.gene_id")
  
  return(top_table_df)
  
}

top_table_df <- make_top_table_long_form_df(fit = fit, n_genes = dim(dge$counts)[1])

# write it out before adding gene names
write.csv(top_table_df, file = toptable_outpath, 
          quote = FALSE, row.names = FALSE)

########################### Plot ###########################
# Compare to POPs2
# read in val DE results
validation_degs_subset <- read.csv(toptable_outpath)


# Subset to genes that were tested in both datasets
# list genes we tested for (gene names)
our_genes <- unique((toptable$gtf.gene_name))
our_gene_ids <- unique((toptable$gtf.gene_id))
print("n unique gene names - POPS")
length(our_genes)
print("n unique gene ids - POPS")
length(our_gene_ids)
# list genes they tested (gene names)
val_genes <- unique(validation_degs_subset$gtf.gene_name)
val_gene_ids <- unique(validation_degs_subset$gtf.gene_id)
print("n unique gene names - GomezLopez")
length(val_genes)
print("n unique gene ids - GomezLopez")
length(val_gene_ids)
# list overlap
gene_id_intersect <- intersect(our_gene_ids, val_gene_ids) 
print("n genes intersecting in both POPS and GomezLopez")
length(gene_id_intersect)

# subset their gene list to overlapping genes 
validation_degs_sub <- validation_degs_subset %>% filter(gtf.gene_id %in% gene_id_intersect)

# make the same type of df for our top table as the one we read in 
# AND subset our gene list to overlapping theirs 
our_degs_for_testing <- toptable %>% filter(gtf.gene_id %in% gene_id_intersect)

# list time point comparisons
tp_comps_for_validation <- unique(our_degs_for_testing$comp)

# joined validation toptable
joined_val_toptable <- our_degs_for_testing %>% 
  mutate(DEG_pops = (adj.P.Val < p_val & abs(logFC) > log_FC)) %>% 
  full_join(validation_degs_sub %>% 
              mutate(DEG_val = (adj.P.Val < p_val & abs(logFC) > log_FC)), 
            #by = c("gtf.gene_id" = "ensembl_gene_id", "comp" = "comp")) %>% 
            by = c("gtf.gene_id", "comp")) %>% 
  mutate(`Gene differentially\nexpressed in` = factor(ifelse(DEG_pops & DEG_val, "Both", 
                                                             ifelse(DEG_pops, "POPS2", 
                                                                    ifelse(DEG_val, "Gomez-Lopez", "Neither"))), 
                                                      levels = c("Both", "POPS2", "Gomez-Lopez", "Neither"))) 

# Check correlations colored by significance (Supplementary Figure 9B)
joined_val_toptable %>% 
  ggplot(aes(x = logFC.x, y = logFC.y)) + 
  geom_point(aes(color = `Gene differentially\nexpressed in`), size = 0.8) + 
  scale_color_manual(values=c("mediumseagreen", "#32769B", "maroon", "gray")) +
  ggpubr::stat_cor(method="pearson", label.y.npc="top", label.x.npc = "left", size = 3) + 
  xlab("LogFC in POPS2 dataset") + 
  ylab("LogFC in Gomez-Lopez 2019 validation dataset") + 
  ggtitle("Correlation of logFC in POPS2 vs. Gomez-Lopez validation dataset") + 
  facet_wrap(~comp, scales = "free", ncol = 2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  geom_hline(yintercept = 0, size = 0.2) + 
  theme(strip.background = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 3))) 

# calculate the numners 
# n genes detected in both cohorts
dim(joined_val_toptable)/6 # div 6 because each gene appears 6 times, one in each comp
# [1] 12814.000000     3.333333
# n genes DE in gomez-lopez dataset 
joined_val_toptable %>% filter(DEG_val == TRUE) %>% dplyr::select(gtf.gene_id) %>% distinct() %>% nrow()
# [1] 48
# n genes we recapitulate
joined_val_toptable %>% filter(DEG_val == TRUE, DEG_pops == TRUE) %>% dplyr::select(gtf.gene_id) %>% distinct() %>% nrow()
# [1] 43
# prop we re-capitulate
43/48
# [1] 0.8958333

# Venn diagram of unique DEGs adjusted for size (Supplementary Figure 9A)
unique_degs_pops <- our_degs_for_testing %>% 
  dplyr::filter(adj.P.Val < p_val, abs(logFC) > log_FC) %>% 
  pull(gtf.gene_id) %>% unique()
length(unique_degs_pops)
unique_degs_gl <- validation_degs_sub %>% 
  dplyr::filter(adj.P.Val < p_val, abs(logFC) > log_FC) %>% 
  pull(gtf.gene_id) %>% unique()
length(unique_degs_gl)
a <- list(POPs2 = unique_degs_pops, `Gomez-Lopez` = unique_degs_gl)

venn_data <- euler(
  c(
    "POPs2" = length(setdiff(a$POPs2, a$`Gomez-Lopez`)),#889,
    "Gomez-Lopez" = length(setdiff(a$`Gomez-Lopez`, a$POPs2)),#5,
    "POPs2&Gomez-Lopez" = length(intersect(a$POPs2, a$`Gomez-Lopez`)) #43
  )
)
plot(
  venn_data,
  fills = c("#32769B", "maroon"), 
  quantities = TRUE
)

# 07_9_characterize_interaction_clusters.R

################################################################################

# 7.9. Characterize eQTL time interaction clusters

################################################################################

# Aim: Identify pathways and regulons enriched among time interaction clusters 

########################### Output paths ##########################

########################### Input paths ###########################
path <- "genotyping/analysis/eQTL/"
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"
RData.location="/path/to/dir/XGR/bigdata"

########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)
library(XGR)
library(decoupleR)
library(lme4)
library(emmeans)

# set ggplot theme 
theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Load data ###########################
# load significant interactions 
time_interactions <- read.csv(paste0(path,"output_data/interactions_flanders_2/time_interactions_flanders_2.csv"))
neutrophil_interactions <- read.csv(paste0(path,"output_data/interactions_flanders_2/neutrophil_interactions_flanders_2.csv"))
monocyte_interactions <- read.csv(paste0(path,"output_data/interactions_flanders_2/monocyte_interactions_flanders_2.csv"))
# conditional eQTL for background 
conditional_eQTL <- readRDS(paste0(path, "output_data/flanders_conditional_eQTL/all_conditional_signals_sign_aligned_pops.rds"))
# gene id_name converter 
gtf <- rtracklayer::import(gtf_inpath)
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)
# read in time interaction clusters
time_int_clust <- read.csv(paste0(path, "output_data/interactions_flanders_2/time_interaction_clusters_2.csv")) %>%
  mutate(cluster_label = cluster_label + 1) %>%
  separate(gsp, sep = "_", into = c("Gene", "SNP"), remove = FALSE, extra = "merge") %>%
  left_join(id_name, by = c("Gene" = "gtf.gene_id"))
table(time_int_clust$cluster_label)
# counts
counts <- read.table("rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt", header = TRUE)
sample.info <- read.csv("rna-seq/data/sample_covariates_clin_tech.csv")

########################### Analysis ###########################

###############################################################################
# Run XGR pathway enrichment
###############################################################################

background <- id_name %>% filter(gtf.gene_id %in% unique(conditional_eQTL$phenotype_id)) %>% pull(gtf.gene_name) %>% unique()
length(background)
# [1] 14262

#First, are all time interactions enriched for any pathways? 
all_time_eQTL_pathway_enrichment <- xEnricherGenes(data=time_interactions %>% left_join(id_name, by = c("Gene" = "gtf.gene_id")) %>% pull(gtf.gene_name), 
                                                   background=background, ontology="REACTOME", ontology.algorithm="none", RData.location=RData.location) #https://xgr.r-forge.r-project.org/
plt <- xEnrichBarplot(all_time_eQTL_pathway_enrichment)

saveRDS(all_time_eQTL_pathway_enrichment, file = paste0(path, "output_data/interactions_flanders_2/xgr_enrichment_all_interactions.rds"))

# Next, are clusters enriched for different pathways? 
list_eTerm_clust = list()
list_consise_eTerm_clust <- list()

for(i in 1:max(time_int_clust$cluster_label)){
  clust = i
  print(paste("Cluster", clust, "Start"))
  data <- time_int_clust %>% filter(cluster_label == clust) %>% pull(gtf.gene_name)
  GO <- xEnricherGenes(data=data, background=background, ontology="REACTOME", ontology.algorithm="none", RData.location=RData.location) #https://xgr.r-forge.r-project.org/
  list_eTerm_clust[[i]] <- GO
  names(list_eTerm_clust)[i] <- paste0("Cluster", clust)
  if(nrow(GO$term_info) > 1){
    list_consise_eTerm_clust[[i]] <- xEnrichConciser(GO)
  }else{
    list_consise_eTerm_clust[[i]] <- GO
  }
  names(list_consise_eTerm_clust)[i] <- paste("Cluster", clust)
  print(paste(clust, "End"))
}

# write this out 
saveRDS(list_consise_eTerm_clust, file = paste0(path, "output_data/interactions_flanders_2/xgr_enrichment_clusters.rds"))

###############################################################################
# Run dorothea regulon enrichment per cluster 
###############################################################################

net <- decoupleR::get_dorothea(organism='human', levels=c('A', 'B', 'C'))

# function to calculate fold change of enrichment 
fcHyper <- function(genes.group, genes.term, genes.universe) {
  genes.hit <- intersect(genes.group, genes.term)
  X <- length(genes.hit)
  K <- length(genes.group)
  M <- length(genes.term)
  N <- length(genes.universe)
  x.exp <- K * M/N
  fc <- X/x.exp
  return(fc)
}

# get all genes for background
genes.using <- id_name %>% filter(gtf.gene_id %in% unique(conditional_eQTL$phenotype_id))
# limit regulons to measured genes
net_sub <- subset(net, target %in% genes.using$gtf.gene_name)
genes.using$regulon <- genes.using$gtf.gene_name %in% net$target

regulon.enrichment <- list()

# Test each cluster for enrichment of each regulon
for(i in 1:max(time_int_clust$cluster_label)){
  clust = i
  print(paste("Cluster", clust, "Start"))
  genes <- time_int_clust %>% filter(cluster_label == clust) %>%
    select(gtf.gene_name, cluster_label)
  
  regulon.enrichment.cluster <- data.frame(Cluster=clust,
                                           TF=unique(net$source),
                                           OR=NA,
                                           Pval=NA, 
                                           fc=NA)
  
  for(r in 1:length(unique(net$source))){
    tf <- unique(net$source)[r]
    regulon <- subset(net, source == tf)
    #print(tf)
    # only test regulons with at least 5 genes assessed
    if(nrow(regulon)>4 & nrow(genes)>4){
      cont.table <- table(unique(genes.using$gtf.gene_name) %in% genes$gtf.gene_name,
                          unique(genes.using$gtf.gene_name) %in% regulon$target)
      # Add a filter for number of genes
      if(dim(cont.table)[2] == 2){
        if(cont.table[2, 2] >2){
          results <- fisher.test(cont.table, alternative = "greater")
          regulon.enrichment.cluster$OR[r] <- results$estimate
          regulon.enrichment.cluster$Pval[r] <- results$p.value
          regulon.enrichment.cluster$fc[r] <- fcHyper(genes$gtf.gene_name, regulon$target, unique(genes.using$gtf.gene_name))
        }
      }
    }
  }
  regulon.enrichment.cluster <- regulon.enrichment.cluster[complete.cases(regulon.enrichment.cluster),]
  regulon.enrichment.cluster$FDR <- p.adjust(regulon.enrichment.cluster$Pval, method="fdr")
  regulon.enrichment[[clust]] <- regulon.enrichment.cluster
}

enrichment.results <- data.table::rbindlist(regulon.enrichment)

enrichment.results$FDRall <- p.adjust(enrichment.results$Pval, method="fdr")
sig.results <- subset(enrichment.results, FDRall<0.05)

write.csv(sig.results, file = paste0(path, "output_data/interactions_flanders_2/time_int_cluster_regulon_enrichment.csv"), row.names = FALSE, quote = FALSE)

# select eQTL downstream of STAT1 and SMAD3 that have eQTL interaction in cluster 4
# Identify list of eQTL with eGenes downstream of TAL1
tf <- "STAT1"
stat1_regulon <- subset(net, source == tf)
genes_c4 <- time_int_clust %>% filter(cluster_label == 4) %>%
  select(gtf.gene_name, cluster_label)
stat1_genes_with_int <- genes.using %>% # out of all conditionally independent eQTL 
  filter(gtf.gene_name %in% c(genes_c4$gtf.gene_name), # filter to time interaction eGenes
         gtf.gene_name %in% stat1_regulon$target) # filter to genes downstream of TAL1
dim(stat1_genes_with_int)

cont.table <- table(in_c4 = unique(genes.using$gtf.gene_name) %in% genes_c4$gtf.gene_name,
                    in_STAT1_regulon = unique(genes.using$gtf.gene_name) %in% stat1_regulon$target)
results <- fisher.test(cont.table, alternative = "greater")

time_int_clust %>% filter(cluster_label == 4, Gene %in% stat1_genes_with_int$gtf.gene_id) %>% 
  left_join(conditional_eQTL %>% select(phenotype_id, snps, SNPpos, chr), by = c("Gene" = "phenotype_id",
                                                                                 "SNP" = "snps", 
                                                                                 "SNPpos")) %>% 
  write.csv(paste0(path, "output_data/interactions_flanders_2/time_int_cluster_4_STAT1_downstream_eQTL.csv"), row.names = FALSE, quote = FALSE)


tf <- "SMAD3"
SMAD3_regulon <- subset(net, source == tf)
SMAD3_genes_with_int <- genes.using %>% # out of all conditionally independent eQTL 
  filter(gtf.gene_name %in% c(genes_c4$gtf.gene_name), # filter to time interaction eGenes
         gtf.gene_name %in% SMAD3_regulon$target) # filter to genes downstream of TAL1
dim(SMAD3_genes_with_int)

cont.table <- table(in_c4 = unique(genes.using$gtf.gene_name) %in% genes_c4$gtf.gene_name,
                    in_SMAD3_regulon = unique(genes.using$gtf.gene_name) %in% SMAD3_regulon$target)
results <- fisher.test(cont.table, alternative = "greater")

time_int_clust %>% filter(cluster_label == 4, Gene %in% SMAD3_genes_with_int$gtf.gene_id) %>% 
  left_join(conditional_eQTL %>% select(phenotype_id, snps, SNPpos, chr), by = c("Gene" = "phenotype_id",
                                                                                 "SNP" = "snps", 
                                                                                 "SNPpos")) %>% 
  write.csv(paste0(path, "output_data/interactions_flanders_2/time_int_cluster_4_SMAD3_downstream_eQTL.csv"), row.names = FALSE, quote = FALSE)

################ look at regulon enrichment in all interactions ###############
genes <- time_int_clust %>%
  select(gtf.gene_name, cluster_label)

regulon.enrichment.all <- data.frame(#Cluster=clust,
  TF=unique(net$source),
  OR=NA,
  Pval=NA, 
  fc=NA)

for(r in 1:length(unique(net$source))){
  tf <- unique(net$source)[r]
  regulon <- subset(net, source == tf)
  #print(tf)
  # only test regulons with at least 5 genes assessed
  if(nrow(regulon)>4 & nrow(genes)>4){
    cont.table <- table(unique(genes.using$gtf.gene_name) %in% genes$gtf.gene_name,
                        unique(genes.using$gtf.gene_name) %in% regulon$target)
    # Add a filter for number of genes
    if(dim(cont.table)[2] == 2){
      if(cont.table[2, 2] >2){
        results <- fisher.test(cont.table, alternative = "greater")
        regulon.enrichment.all$OR[r] <- results$estimate
        regulon.enrichment.all$Pval[r] <- results$p.value
        regulon.enrichment.all$fc[r] <- fcHyper(genes$gtf.gene_name, regulon$target, unique(genes.using$gtf.gene_name))
      }
    }
  }
}
regulon.enrichment.all <- regulon.enrichment.all[complete.cases(regulon.enrichment.all),]
regulon.enrichment.all$FDR <- p.adjust(regulon.enrichment.all$Pval, method="fdr")
regulon.enrichment.all %>% filter(FDR < 0.05)

write.csv(regulon.enrichment.all %>% filter(FDR < 0.05), 
          file = paste0(path, "output_data/interactions_flanders_2/time_int_all_regulon_enrichment.csv"), row.names = FALSE, quote = FALSE)

###############################################################################
# Look at dorothea regulon activity in SMAD3 and STAT1
###############################################################################

regulon.enrichment.cluster <- read.csv(paste0(path, "output_data/interactions_flanders_2/time_int_cluster_regulon_enrichment.csv"))

# first, change rownames of expression data to gene symbols
counts <- counts %>% rownames_to_column(var = "gtf.gene_id") %>%
  left_join(id_name, by = "gtf.gene_id") %>%
  add_count(gtf.gene_name) %>%   # adds a column n with number of times each gene name occurs
  filter(n == 1) %>%             # keep only unique gene names
  select(-n, -gtf.gene_id) %>%                  # drop the count column
  column_to_rownames(var = "gtf.gene_name")
counts <- as.matrix(counts)

net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
net_sub <- net %>% 
  # subset net to only contain genes in our count matrix
  filter(target %in% rownames(counts),
         # subset source to be only those of interest
         source %in% regulon.enrichment.cluster$TF)

sample_acts_ulm <- decoupleR::run_ulm(
  mat = counts,
  network = net_sub,
  .source = "source",
  .target = "target",
  #consensus_score = TRUE,
  minsize = 5)

# took 6 mins to run for 2 regulons

# save this 
sample_acts_ulm %>% write.csv(paste0(path, "output_data/interactions_flanders_2/SMAD3_STAT1_regulon_activity.csv"), quote = FALSE, row.names = FALSE)
# read it back in 
sample_acts_ulm <- read.csv(paste0(path, "output_data/interactions_flanders_2/SMAD3_STAT1_regulon_activity.csv")) %>%
  left_join(sample.info %>% select(RNA_sanger_sample_id, Sample_taken_at, ANON_ID), by = c("condition" = "RNA_sanger_sample_id"))

# n significant
sample_acts_ulm %>% 
  group_by(source) %>%
  summarize(n_sig = sum(p_value < 0.05)/n())

sample_acts_ulm %>%
  ggplot(aes(x = score, y = p_value, color = Sample_taken_at)) +
  geom_point() + 
  facet_wrap(~source, scales = "free") + 
  geom_hline(yintercept = 0.05)

# test for differences between time points (paired)
model_smad3 <- lmer(score ~ Sample_taken_at + (1 | ANON_ID),
                    data = sample_acts_ulm %>% filter(source == "SMAD3"))


model_stat1 <- lmer(score ~ Sample_taken_at + (1 | ANON_ID),
                    data = sample_acts_ulm %>% filter(source == "STAT1"))

anova(model_smad3)
anova(model_stat1)

emm_smad3 <- emmeans(model_smad3, ~ Sample_taken_at)
emm_stat1 <- emmeans(model_stat1, ~ Sample_taken_at)

smad3_diffs <- pairs(emm_smad3, adjust = "BH")
smad3_diffs %>% write.csv(paste0(path, "output_data/interactions_flanders_2/SMAD3_change_in_regulon_activity.csv"), row.names = FALSE, quote = FALSE)
stat1_diffs <- pairs(emm_stat1, adjust = "BH")
stat1_diffs %>% write.csv(paste0(path, "output_data/interactions_flanders_2/STAT1_change_in_regulon_activity.csv"), row.names = FALSE, quote = FALSE)

########################### Plot ###########################
# xgr enriched pathways 
list_consise_eTerm_clust <- readRDS(paste0(path, "output_data/interactions_flanders_2/xgr_enrichment_clusters.rds"))
all_time_eQTL_pathway_enrichment <- readRDS(paste0(path, "output_data/interactions_flanders_2/xgr_enrichment_all_interactions.rds"))
# dorothea enriched regulons
regulon_enrichment_all_time_interactions <- read.csv(paste0(path, "output_data/interactions_flanders_2/time_int_all_regulon_enrichment.csv"))
regulon_enrichment_cluster_time_interactions <- read.csv(paste0(path, "output_data/interactions_flanders_2/time_int_cluster_regulon_enrichment.csv"))
# regulon activity scores
regulon_activity_scores <- read.csv(paste0(path, "output_data/interactions_flanders_2/SMAD3_STAT1_regulon_activity.csv")) %>%
  left_join(sample.info %>% dplyr::select(RNA_sanger_sample_id, Sample_taken_at, ANON_ID), by = c("condition" = "RNA_sanger_sample_id"))

# Plot pathway enrichment of time interaction clusters
xEnrichBarplot_update <- function(eTerm, top_num = 10, displayBy = c("fc", "adjp", "fdr", 
                                                                     "zscore", "pvalue"), FDR.cutoff = 0.05, bar.label = TRUE, 
                                  bar.label.size = 3, bar.color = "lightyellow-orange", bar.width = 0.8, 
                                  wrap.width = NULL, font.family = "sans", signature = TRUE) 
{
  displayBy <- match.arg(displayBy)
  if (is.null(eTerm)) {
    warnings("There is no enrichment in the 'eTerm' object.\n")
    return(NULL)
  }
  df <- xEnrichViewer(eTerm, top_num = "all")
  if (top_num == "auto") {
    top_num <- sum(df$adjp < FDR.cutoff)
    if (top_num <= 1) {
      top_num <- 10
    }
  }
  df <- xEnrichViewer(eTerm, top_num = top_num, sortBy = "adjp")
  if (!is.null(wrap.width)) {
    width <- as.integer(wrap.width)
    res_list <- lapply(df$name, function(x) {
      x <- gsub("_", " ", x)
      y <- strwrap(x, width = width)
      if (length(y) > 1) {
        paste0(y[1], "...")
      }
      else {
        y
      }
    })
    df$name <- unlist(res_list)
  }
  name <- height <- direction <- hjust <- NULL
  adjp <- zscore <- pvalue <- fc <- NULL
  df <- df %>% dplyr::mutate(direction = ifelse(zscore > 0, 
                                                1, -1))
  if (displayBy == "adjp" | displayBy == "fdr") {
    df <- df %>% dplyr::arrange(direction, desc(adjp), zscore) %>% 
      dplyr::mutate(height = -1 * log10(adjp)) %>% dplyr::mutate(hjust = 1)
    df$name <- factor(df$name, levels = df$name)
    if (length(df$height[!is.infinite(df$height)]) == 0) {
      df$height <- 10
    }
    else {
      df$height[is.infinite(df$height)] <- max(df$height[!is.infinite(df$height)])
    }
    p <- ggplot(df, aes(x = name, y = height))
    p <- p + ylab(expression(paste("Enrichment significance: ", 
                                   -log[10]("FDR"))))
  }
  else if (displayBy == "fc") {
    df <- df %>% dplyr::arrange(direction, fc, desc(adjp)) %>% 
      dplyr::mutate(height = log2(fc)) %>% dplyr::mutate(hjust = ifelse(height >= 
                                                                          0, 1, 0))
    df$name <- factor(df$name, levels = df$name)
    p <- ggplot(df, aes(x = name, y = height))
    p <- p + ylab(expression(paste("Enrichment changes: ", 
                                   log[2]("FC"))))
  }
  else if (displayBy == "pvalue") {
    df <- df %>% dplyr::arrange(direction, desc(pvalue), 
                                zscore) %>% dplyr::mutate(height = -1 * log10(pvalue)) %>% 
      dplyr::mutate(hjust = 1)
    df$name <- factor(df$name, levels = df$name)
    if (length(df$height[!is.infinite(df$height)]) == 0) {
      df$height <- 10
    }
    else {
      df$height[is.infinite(df$height)] <- max(df$height[!is.infinite(df$height)])
    }
    p <- ggplot(df, aes(x = name, y = height))
    p <- p + ylab(expression(paste("Enrichment significance: ", 
                                   -log[10]("p-value"))))
  }
  else if (displayBy == "zscore") {
    df <- df %>% dplyr::arrange(direction, zscore, desc(adjp)) %>% 
      dplyr::mutate(height = zscore) %>% dplyr::mutate(hjust = ifelse(height >= 
                                                                        0, 1, 0))
    df$name <- factor(df$name, levels = df$name)
    p <- ggplot(df, aes(x = name, y = height))
    p <- p + ylab("Enrichment z-scores")
  }
  if (is.null(bar.color)) {
    bp <- p + geom_col(color = "grey80", fill = "transparent", 
                       width = bar.width)
  }
  else {
    bar.color <- unlist(strsplit(bar.color, "-"))
    if (length(bar.color) == 2) {
      bp <- p + geom_col(aes(fill = height), width = bar.width) + 
        scale_fill_gradient(low = bar.color[1], high = bar.color[2])
    }
    else {
      bp <- p + geom_col(color = "grey80", fill = "transparent", 
                         width = bar.width)
    }
  }
  bp <- bp + theme_bw() + theme(legend.position = "none", 
                                axis.title.y = element_blank(), axis.text.y = element_text(size = 10, 
                                                                                           color = "black"), axis.title.x = element_text(size = 12, 
                                                                                                                                         color = "black")) + coord_flip()
  bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (bar.label) {
    to_scientific_notation <- function(x) {
      res <- format(x, scientific = T)
      res <- sub("\\+0?", "", res)
      sub("-0?", "-", res)
    }
    label <- to_scientific_notation(df$adjp)
    df$label <- paste("FDR", as.character(label), sep = "=")
    df_sub <- df %>% dplyr::filter(hjust == 1)
    bp <- bp + geom_text(data = df_sub, aes(label = label), 
                         hjust = 1, size = bar.label.size, family = font.family)
    df_sub <- df %>% dplyr::filter(hjust == 0)
    bp <- bp + geom_text(data = df_sub, aes(label = label), 
                         hjust = 0, size = bar.label.size, family = font.family)
  }
  if (signature) {
    caption <- paste("Created by xEnrichBarplot from XGR version", 
                     utils::packageVersion("XGR"))
    bp <- bp + labs(caption = caption) + theme(plot.caption = element_text(hjust = 1, 
                                                                           face = "bold.italic", size = 8, colour = "#002147"))
  }
  bp <- bp + theme(text = element_text(family = font.family))
  bp <- bp + theme(axis.line.x = element_line(arrow = arrow(angle = 30, 
                                                            length = unit(0.25, "cm"), type = "open")))
  bp <- bp + scale_y_continuous(position = "left")
  invisible(bp)
}
# Pathway enrichment of all time interactions (Supplementary figure 18A)
all_int_XGR <- xEnrichBarplot_update(all_time_eQTL_pathway_enrichment, bar.color = "gray-gray", signature = FALSE)
all_int_XGR

# Dotplot of pathway enrichment of time interaction clusters
# Extract data used for plotting
get_xEnrichCompare_data <- function(list_eTerm, FDR.cutoff = 0.05, displayBy = "fc") {
  list_eTerm <- base::Filter(base::Negate(is.null), list_eTerm)
  list_names <- names(list_eTerm)
  if (is.null(list_names)) {
    list_names <- paste("Enrichment", 1:length(list_eTerm), sep = " ")
  }
  res_ls <- lapply(1:length(list_eTerm), function(i) {
    df <- xEnrichViewer(list_eTerm[[i]], top_num = "all", sortBy = "none")
    if (is.null(df)) return(NULL)
    cbind(group = rep(list_names[i], nrow(df)), id = rownames(df), df, stringsAsFactors = FALSE)
  })
  df_all <- do.call(rbind, res_ls)
  ind <- which(df_all$adjp < FDR.cutoff)
  d <- df_all[ind, c("id", "name", "fc", "adjp", "zscore", "pvalue", "group", "or")]
  d$group <- factor(d$group, levels = rev(list_names))
  d <- d %>% dplyr::mutate(direction = ifelse(zscore > 0, 1, -1))
  d$name <- factor(d$name, levels = unique(d$name))
  return(d)
}

df_dotplot <- get_xEnrichCompare_data(list_consise_eTerm_clust, FDR.cutoff = 0.05) 

# Get original order BEFORE modification
original_order <- levels(df_dotplot$name)

name_order <- levels(df_dotplot$name)
wrapped_names <- stringr::str_wrap(name_order, width = 46)
name_map <- setNames(wrapped_names, name_order)
df_dotplot <- df_dotplot %>%
  mutate(
    name_wrapped = name_map[as.character(name)],
    name_wrapped = factor(name_wrapped, levels = wrapped_names)
  )

# Add names of clusters that don't have any enrichment 1, 5
new_group <- "Cluster t1"
new_group_2 <- "Cluster t5"

# All unique pathways in your current df_dotplot
all_pathways <- unique(df_dotplot$name)

# Create empty rows for the new group
empty_rows <- data.frame(
  name = all_pathways,
  group = factor(rep(new_group, length(all_pathways)), levels = levels(df_dotplot$group)),
  fc = NA_real_,
  adjp = NA_real_,
  pvalue = NA_real_,
  zscore = NA_real_
)

empty_rows_2 <- data.frame(
  name = all_pathways,
  group = factor(rep(new_group_2, length(all_pathways)), levels = levels(df_dotplot$group)),
  fc = NA_real_,
  adjp = NA_real_,
  pvalue = NA_real_,
  zscore = NA_real_
)


# Bind to original df_dotplot
df_dotplot_extended <- dplyr::bind_rows(df_dotplot, empty_rows, empty_rows_2) %>%
  mutate(# Map original name to wrapped label (use as.character to avoid factor issues)
    name_wrapped = name_map[as.character(name)],
    # Convert to factor with wrapped levels to preserve order
    name_wrapped = factor(name_wrapped, levels = wrapped_names))

pathway_order_wrapped <- df_dotplot_extended %>%
  filter(!is.na(fc)) %>%
  count(name_wrapped, name = "n_groups") %>%
  arrange(desc(n_groups)) %>%
  pull(name_wrapped)

# Pathway enrichment dotplot (Main Figure 6C)
df_dotplot_extended %>% 
  mutate(
    name_wrapped = factor(name_wrapped, levels = pathway_order_wrapped), 
    log2fc = log2(fc), 
    #group = factor(group, levels = c("36wkGA", "28wkGA", "20wkGA", "12wkGA")))
    ) %>%
  ggplot(aes(x = name_wrapped, y = group)) +
  geom_point(aes(size = -log10(adjp), fill = log2fc), color = "black", shape = 21) +
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    name = "log2 fold\nenrichment"
  ) +
  scale_size_continuous(
    name = "-log10(FDR)"
  ) +
  labs(
    size = "-log10(FDR)",
    x = "Pathway",
    y = "Pregnancy time-point"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color="black", size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.title = element_text(size = 9),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  # add extra space above#,
    #legend.position = "top"
  ) + 
  geom_vline( xintercept = 2 + 0.5 ) + 
  geom_vline( xintercept = 11 + 0.5 ) + 
  geom_vline( xintercept = 14 + 0.5 ) + 
  annotate(
    "text",
    x = 1.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Shared",
    size = 3
  ) +
  annotate(
    "text",
    x = 2+(9/2)+0.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Cluster t2",
    size = 3
  ) +
  annotate(
    "text",
    x = 11+(3/2)+0.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Cluster t3",
    size = 3
  ) +
  annotate(
    "text",
    x = 14+(2/2)+0.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Cluster t4",
    size = 3
  ) +
  coord_cartesian(clip = "off")

# Regulon enrichments of time interaction clusters
# (Supplementary Figure 20)
regulon_enrichment_cluster_time_interactions %>%
  mutate(Cluster = factor(Cluster, levels = rev(1:5)), 
         log2fc = log2(fc)) %>%
  ggplot(aes(x = Cluster, y = TF)) + 
  geom_point(aes(size = -log10(FDRall), fill = log2fc),
             color = "black", shape = 21) + 
  scale_x_discrete(drop = FALSE) +   # 👈 key line
  scale_size_continuous(breaks = c(1, 1.5, 2, 2.5, 3)) +
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    na.value = "gray",
    #name = "Odds ratio"
    name = "log2 fold\nenrichment"
  ) +
  coord_flip() +
  labs(
    size = "-log10(FDR)",
    x = "Time interaction eQTL cluster",
    y = "Regulon (TF)",
    title = "Enriched regulons\nin eQTL time interaction clusters"
  )

# Regulon enrichments of all time interactions
# (Supplementary Figure 18B)
regulon_enrichment_all_time_interactions %>%
  mutate(Cluster = "All time-point ieQTL", 
         log2fc = log2(fc)) %>%
  #mutate(Cluster = factor(Cluster, levels = 1:5)) %>%
  ggplot(aes(x = Cluster, y = TF)) + 
  geom_point(aes(size = -log10(FDR), fill = log2fc),
             color = "black", shape = 21) + 
  scale_x_discrete(drop = FALSE) +  
  scale_size_continuous(
    range = c(2, 4),   # tweak these numbers
    breaks = c(1.75, 1.8, 1.85)
  ) + 
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    na.value = "gray",
    #name = "Odds ratio"
    name = "log2 fold\nenrichment",
    limits = c(0.1, 0.4)
  ) +
  coord_flip() +
  labs(
    size = "-log10(FDR)",
    x = " ",
    y = "Regulon (TF)",
    title = "Enriched regulons\nin all eQTL time interactions"
  )

# Plot regulon enrichment scores of SMAD3 and STAT1 over time 
# (Supplementary figure 21)
regulon_activity_scores %>%
  mutate(Sample_taken_at = gsub("_weeks", "", Sample_taken_at)) %>%
  ggplot(aes(x = Sample_taken_at, y = score)) +
  geom_violin(aes(fill = Sample_taken_at)) + 
  geom_boxplot(fill = "transparent", width = 0.1, color = "black") + 
  facet_wrap(~source, scales = "free") + 
  scale_fill_manual(values=c("#BCE4D8", "#83C4CB", "#439FB7", "#32769B")) + 
  xlab("Pregnancy time-point (weeks)") + 
  ylab("decoupleR regulon activity score") + 
  theme(legend.position = "none")
# Add significance bars with info from Supplementary Table 9 (SMAD3_change_in_regulon_activity.csv, STAT1_change_in_regulon_activity.csv)


# 06_4_characterizing_pregnancy_magnified_egenes.R

################################################################################

# 6.4. Characterize pregnancy-magnified eGenes

################################################################################

# Aim: run pathway and regulon enrichment on pregnancy-magnified eGenes to understand their functions

########################### Output paths ##########################
dir <- "genotyping/analysis/eQTL/output_data/"

########################### Input paths ###########################
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Load packages ###########################
library(tidyverse)
library(XGR)
library(decoupleR)

########################### Load data and Analysis ###########################

# load data
# main effects for background 
main_effects_12wk <-  readRDS(paste0(dir, "single_tp_analysis/12_weeks/12_weeks_ciseqtl_all.rds"))
# id name conversion
gtf <- rtracklayer::import(gtf_inpath)
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)
# mashR results
mashr.results <- read.table(paste0(dir, "interval_comp/pops_alltps_intervalFRA_mashr_results_sharing_filtered_leads.txt"), sep="\t") %>%
  separate(Pairs_notrsID, into = c("chr", "SNPpos", "gtf.gene_id"), sep = "_", remove = FALSE) %>%
  mutate(shared = factor(shared, levels = c("Bigger effect in POPS2", "Only significant in POPS2", "Bigger effect in Interval", "Opposite direction of effect", "shared", "not significant")),
         SNPpos = as.integer(SNPpos)) %>%
  left_join(id_name, by = "gtf.gene_id") %>%
  left_join(main_effects_12wk %>% select(snps, gene, SNPpos), by = c("gtf.gene_id" = "gene", "SNPpos"))

################################################################################
# First, get groups of mashR signals from the characterization in 6.3
################################################################################

# pregnancy-magnified signals per time-point
pregnancy_magnified_12wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_12wk",
         shared %in% c("Bigger effect in POPS2", "Only significant in POPS2"))
pregnancy_magnified_20wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_20wk",
         shared %in% c("Bigger effect in POPS2", "Only significant in POPS2")) 
pregnancy_magnified_28wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_28wk",
         shared %in% c("Bigger effect in POPS2", "Only significant in POPS2")) 
pregnancy_magnified_36wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_36wk",
         shared %in% c("Bigger effect in POPS2", "Only significant in POPS2")) 

# pregnancy-dampened signals per time-point
pregnancy_dampened_12wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_12wk",
         shared  == "Bigger effect in Interval")
pregnancy_dampened_20wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_20wk",
         shared  == "Bigger effect in Interval")
pregnancy_dampened_28wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_28wk",
         shared  == "Bigger effect in Interval")
pregnancy_dampened_36wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_36wk",
         shared  == "Bigger effect in Interval")

# opposite direction of effect signals per time-point
op_dir_eff_12wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_12wk",
         shared  == "Opposite direction of effect")
op_dir_eff_20wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_20wk",
         shared  == "Opposite direction of effect")
op_dir_eff_28wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_28wk",
         shared  == "Opposite direction of effect") 
op_dir_eff_36wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_36wk",
         shared  == "Opposite direction of effect")

# shared signals per time-point
shared_12wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_12wk",
         shared  == "shared")
shared_20wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_20wk",
         shared  == "shared")
shared_28wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_28wk",
         shared  == "shared")
shared_36wk <- mashr.results %>% 
  filter(POPS2_timepoint == "POPS2_36wk",
         shared  == "shared")

################################################################################
# XGR on pregnancy-enhanced eQTL 
################################################################################

context_signals_list_genes <- list(pregnancy_magnified_12wk = unique(pregnancy_magnified_12wk$gtf.gene_name),
                                   pregnancy_magnified_20wk = unique(pregnancy_magnified_20wk$gtf.gene_name), 
                                   pregnancy_magnified_28wk = unique(pregnancy_magnified_28wk$gtf.gene_name), 
                                   pregnancy_magnified_36wk = unique(pregnancy_magnified_36wk$gtf.gene_name))

tp_background_gene_list = list(bg_12wk = mashr.results %>% filter(POPS2_timepoint == "POPS2_12wk") %>% pull(gtf.gene_name) %>% unique(),
                               bg_20wk = mashr.results %>% filter(POPS2_timepoint == "POPS2_20wk") %>% pull(gtf.gene_name) %>% unique(),
                               bg_28wk = mashr.results %>% filter(POPS2_timepoint == "POPS2_28wk") %>% pull(gtf.gene_name) %>% unique(),
                               bg_36wk = mashr.results %>% filter(POPS2_timepoint == "POPS2_36wk") %>% pull(gtf.gene_name) %>% unique())

# for pregnancy specific genes
list_eTerm_clust = list()
list_consise_eTerm_clust <- list()

for(i in 1:length(context_signals_list_genes)){
  # set name of the signal group
  nm <- names(context_signals_list_genes)[i]
  cat("Start:", i, nm, "\n")
  # set background
  if(grepl("36", nm)){
    bg <- "bg_36wk"
  } else if(grepl("28", nm)) {
    bg <- "bg_28wk"
  } else if(grepl("20", nm)) {
    bg <- "bg_20wk"
  } else {
    bg <- "bg_12wk"
  }
  print(bg)
  print(length(tp_background_gene_list[[bg]]))
  print(length(context_signals_list_genes[[i]]))
  GO <- xEnricherGenes(data=context_signals_list_genes[[i]], 
                       background=tp_background_gene_list[[bg]], 
                       ontology="REACTOME", 
                       ontology.algorithm="none") #https://xgr.r-forge.r-project.org/
  if(!is.null(GO)){
    list_eTerm_clust[[i]] <- GO
    names(list_eTerm_clust)[i] <- nm
    if(nrow(GO$term_info) > 1){
      list_consise_eTerm_clust[[i]] <- xEnrichConciser(GO)
    }else{
      list_consise_eTerm_clust[[i]] <- GO
    }
    names(list_consise_eTerm_clust)[i] <- nm
  }
  print(paste("End", nm))
}

# enrichment of pregnancy-magnified genes
bp_Pathway_smaller_clusters_enhanced <- xEnrichCompare(list_consise_eTerm_clust, 
                                                       displayBy="fc",
                                                       FDR.cutoff=0.05,
                                                       facet = TRUE,
                                                       wrap.width = 60,
                                                       bar.label = TRUE,
                                                       #bar.label.size = 1,
                                                       sharings = c(1:4))
bp_Pathway_smaller_clusters_enhanced

# write out results
saveRDS(list_consise_eTerm_clust, file = paste0(dir, "interval_comp/context_group_XGR_pathway_enrichment_preg_mag.rds"))

################################################################################
# regulon enrichment on pregnancy-enhanced eQTL
################################################################################

net <- decoupleR::get_dorothea(organism='human', levels=c('A', 'B', 'C'))

# limit regulons to measured genes
net_subs <- list(net_sub_12wk = subset(net, target %in% tp_background_gene_list[["bg_12wk"]]),
                 net_sub_20wk = subset(net, target %in% tp_background_gene_list[["bg_20wk"]]),
                 net_sub_28wk = subset(net, target %in% tp_background_gene_list[["bg_28wk"]]),
                 net_sub_36wk = subset(net, target %in% tp_background_gene_list[["bg_36wk"]]))

regulon.enrichment.eQTL_context_all <- list()

preg_magnified_signals_list_genes <- context_signals_list_genes[1:4]

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

# Test each context group for enrichment of each regulon
for(i in 1:length(preg_magnified_signals_list_genes)){
  gene_list = names(preg_magnified_signals_list_genes)[i]
  print(paste("Gene list", gene_list, "Start"))
  # get the genes in this group
  genes <-  preg_magnified_signals_list_genes[[i]] 
  # set background
  if(grepl("36", gene_list)){
    bg <- "bg_36wk"
  } else if(grepl("28", gene_list)) {
    bg <- "bg_28wk"
  } else if(grepl("20", gene_list)) {
    bg <- "bg_20wk"
  } else {
    bg <- "bg_12wk"
  }
  print(bg)
  
  if(grepl("36", gene_list)){
    net_name <- "net_sub_36wk"
  } else if(grepl("28", gene_list)) {
    net_name <- "net_sub_28wk"
  } else if(grepl("20", gene_list)) {
    net_name <- "net_sub_20wk"
  } else {
    net_name <- "net_sub_12wk"
  }
  print(net_name)
  
  regulon.enrichment.eQT_context <- data.frame(eQT_context=gene_list,
                                               TF=unique(net_subs[[net_name]]$source),
                                               OR=NA,
                                               Pval=NA,
                                               fc=NA)
  
  for(r in 1:length(unique(net_subs[[net_name]]$source))){
    tf <- unique(net_subs[[net_name]]$source)[r]
    regulon <- subset(net_subs[[net_name]], source == tf)
    #print(tf)
    # only test regulons with at least 5 genes assessed
    if(nrow(regulon)>4 & length(genes)>4){
      cont.table <- table(unique(tp_background_gene_list[[bg]]) %in% genes, # which genes tested at this tp are in the list of genes in this context group 
                          unique(tp_background_gene_list[[bg]]) %in% regulon$target) # which genes tested in this tp are in the regulon downstream
      # Add a filter for number of genes
      if(dim(cont.table)[2] == 2){
        if(cont.table[2, 2] >2){
          results <- fisher.test(cont.table, alternative = "greater")
          regulon.enrichment.eQT_context$OR[r] <- results$estimate
          regulon.enrichment.eQT_context$Pval[r] <- results$p.value
          regulon.enrichment.eQT_context$fc[r] <- fcHyper(genes, regulon$target, unique(tp_background_gene_list[[bg]]))
        }
      }
    }
  }
  # remove rows where the regulon couldn't be tested 
  regulon.enrichment.eQT_context <- regulon.enrichment.eQT_context[complete.cases(regulon.enrichment.eQT_context),]
  # do multiple testing correction within the context group
  regulon.enrichment.eQT_context$FDR <- p.adjust(regulon.enrichment.eQT_context$Pval, method="fdr")
  # add the results to the gene list 
  regulon.enrichment.eQTL_context_all[[gene_list]] <- regulon.enrichment.eQT_context
}

enrichment.results <- data.table::rbindlist(regulon.enrichment.eQTL_context_all)
# multiple testing correction 
enrichment.results$FDRall <- p.adjust(enrichment.results$Pval, method="fdr")
sig.results <- subset(enrichment.results, FDRall<0.05)

# write out results
sig.results %>% write.csv(paste0(dir, "interval_comp/context_group_regulon_enrichment_magnified_only.csv"), quote = FALSE, row.names = FALSE)

sig.results %>%
  mutate(eQT_context = as.factor(eQT_context)) %>%
  ggplot(aes(x = eQT_context, y = TF)) + 
  geom_point(aes(size = -log10(FDRall), fill = OR), color = "black", shape = 21) + 
  scale_size_continuous(
    breaks = c(1, 1.5, 2, 2.5, 3)
  ) + 
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    na.value = "gray",
    #limits = color_limits,
    name = "Odds ratio"
  ) +
  coord_flip() +
  labs(
    size = "-log10(FDR)",
    x = "",
    y = "Regulon (TF)"
  ) + 
  theme(legend.position = "top")

########################### Plot ###########################

# read in pathway enrichment
list_enhanced <- readRDS(file = paste0(dir, "interval_comp/context_group_XGR_pathway_enrichment_preg_mag.rds"))
names(list_enhanced) <- c("12 weeks", "20 weeks", "28 weeks", "36 weeks")
# read in regulon enrichment
regulon_enrichment_preg_enhanced <- read.csv(paste0(dir, "interval_comp/context_group_regulon_enrichment_magnified_only.csv"))


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

df_dotplot <- get_xEnrichCompare_data(list_enhanced, FDR.cutoff = 0.05) 
# check correlation between fc and OR
df_dotplot %>% ggplot(aes(x = or, y = fc)) + 
  geom_point() + ggpubr::stat_cor()

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


pathway_order_wrapped <- df_dotplot %>%
  filter(!is.na(fc)) %>%
  count(name_wrapped, name = "n_groups") %>%
  arrange(desc(n_groups)) %>%
  pull(name_wrapped)

# set limits so that scales are consistent across pathway and regulon enrichment plots
# Combined limits for fill
fill_limits <- range(
  c(
    log2(regulon_enrichment_preg_enhanced$fc),
    log2(df_dotplot$fc)
  ),
  na.rm = TRUE
)

# Combined limits for size
size_limits <- range(
  c(
    -log10(regulon_enrichment_preg_enhanced$FDRall),
    -log10(df_dotplot$adjp)
  ),
  na.rm = TRUE
)

size_breaks <- c(1, 1.5, 2)    # pick breaks that span the combined data
fill_breaks <- pretty(fill_limits, 5)

# Pathway enrichment dotplot (Main Figure 5D)
df_dotplot %>% 
  mutate(#group = factor(group, levels = c("12v20", "12v28", "12v36", "20v28", "20v36", "28v36")),
    name_wrapped = factor(name_wrapped, levels = pathway_order_wrapped), 
    log2fc = log2(fc), 
    group = gsub(" weeks", "wkGA", group), 
    group = factor(group, levels = c("36wkGA", "28wkGA", "20wkGA", "12wkGA"))) %>%
  ggplot(aes(x = name_wrapped, y = group)) +
  geom_point(aes(size = -log10(adjp), fill = log2fc), color = "black", shape = 21) +
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    limits = fill_limits,
    breaks = fill_breaks,
    name = "log2 fold\nenrichment"
  ) +
  scale_size_continuous(
    limits = size_limits,
    breaks = size_breaks,
    name = "-log10(FDR)"
  ) +
  #coord_flip() +
  labs(
    size = "-log10(FDR)",
    x = "Pathway",
    y = "Pregnancy time-point"
  ) +
  #ggtitle("Pregnancy-enhanced eGene pathway enrichment")  +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color="black", size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.title = element_text(size = 9),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  # add extra space above#,
    #legend.position = "top"
  ) + 
  geom_vline( xintercept = 5 + 0.5 ) + 
  geom_vline( xintercept = 9 + 0.5 ) + 
  geom_vline( xintercept = 14 + 0.5 ) + 
  geom_vline( xintercept = 25 + 0.5 ) + 
  annotate(
    "text",
    x = 2.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Common",
    size = 3
  ) +
  annotate(
    "text",
    x = 7.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "12wkGA",
    size = 3
  ) +
  annotate(
    "text",
    x = 12,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "20wkGA",
    size = 3
  ) +
  annotate(
    "text",
    x = 20,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "28wkGA",
    size = 3
  ) +
  annotate(
    "text",
    x = 26,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "36wkGA",
    size = 3
  ) +
  coord_cartesian(clip = "off")

# now, plot regulon enrichment
# check correlation between fc and OR
regulon_enrichment_preg_enhanced %>% ggplot(aes(x = OR, y = fc)) + 
  geom_point() + ggpubr::stat_cor()

# Regulon enrichment dotplot (Main Figure 5C)
regulon_enrichment_preg_enhanced %>%
  mutate(eQTL_context_long = paste0(gsub("pregnancy_magnified_", "", eQT_context), "GA"), 
         eQTL_context_long = factor(eQTL_context_long, levels = c("36wkGA", "28wkGA", "20wkGA", "12wkGA")),
         log2fc = log2(fc)) %>%
  ggplot(aes(x = eQTL_context_long, y = TF)) + 
  geom_point(aes(size = -log10(FDRall), fill = log2fc), color = "black", shape = 21) + 
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    limits = fill_limits,
    breaks = fill_breaks,
    name = "log2 fold\nenrichment"
  ) +
  scale_size_continuous(
    limits = size_limits,
    breaks = size_breaks,
    name = "-log10(FDR)"
  ) +
  coord_flip() +
  labs(
    size = "-log10(FDR)",
    x = "Pregnancy time-point",
    y = "Regulon (TF)"
  ) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, color="black"))

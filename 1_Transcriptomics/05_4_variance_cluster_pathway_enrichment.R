# 05_4_variance_cluster_pathway_enrichment.R 

################################################################################

# 5.4. Pathway enrichment analysis of variance clusters 

################################################################################

# Aim: Identify pathways enriched among variance clusters and plot variance clusters

########################### Output paths ##########################

########################### Input paths ###########################
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"
# variance cluster inputs
inpath <- "rna-seq/analysis/variance/results/"
# gene variances
gene_var_inpath="rna-seq/analysis/variance/results/gene_variance_mean_pertimepoint.tsv"

########################### Parameters ############################
# Data location for XGR 
RData.location="/software/team282/XGR/bigdata"

########################### Load packages ###########################
library(tidyverse)
library(XGR)

theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Load data ###########################
# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
gene_variance <- read.delim(gene_var_inpath)

########################### Analysis ###########################

# Prepare gene ID key 
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

cluster_genes <- c()
k = 6
# variance clusters
for(i in 0:(k-1)){
  clust_genes_tmp <- read.delim(paste0(inpath, "variance_cluster_genes_sh50/cluster", i, ".txt"), header = FALSE) %>%
    rename("gtf.gene_id" = "V1") %>%
    left_join(id_name, by = c("gtf.gene_id")) %>%
    mutate(cluster = i + 1)
  
  cluster_genes <- rbind(cluster_genes, clust_genes_tmp)
}

# write this out for future use
cluster_genes %>% distinct() %>% write.csv(paste0(inpath, "cluster_gene_summary_sh50.csv"), row.names = FALSE, quote = FALSE)
cluster_genes <- read.csv(paste0(inpath, "cluster_gene_summary_sh50.csv"))

# add cluster variance
cluster_variance <- gene_variance %>%
  # add column with cluster name
  left_join(cluster_genes, by = c("gene" = "gtf.gene_id"))

################################################################################
# Pathway enrichment - variance clusters 
################################################################################
list_eTerm_clust = list()
list_consise_eTerm_clust <- list()

background <- id_name[id_name$gtf.gene_id %in% unique(gene_variance$gene),2]
for(i in 1:6){
  clust = i
  print(paste("Cluster", clust, "Start"))
  data_table <- cluster_genes %>% filter(cluster == clust)
  data <- data_table$gtf.gene_name
  GO <- xEnricherGenes(data=data, background=background, ontology="REACTOME", ontology.algorithm="none", 
                       RData.location=RData.location) #https://xgr.r-forge.r-project.org/
  list_eTerm_clust[[i]] <- GO
  names(list_eTerm_clust)[i] <- paste("Cluster", clust)
  if(nrow(GO$term_info) > 1){
    list_consise_eTerm_clust[[i]] <- xEnrichConciser(GO)
  }else{
    list_consise_eTerm_clust[[i]] <- GO
  }
  names(list_consise_eTerm_clust)[i] <- paste("Cluster", clust)
  print(paste(clust, "End"))
}


bp_Pathway_smaller_clusters <- xEnrichCompare(list_consise_eTerm_clust, 
                                              displayBy="fc",
                                              FDR.cutoff=5e-2,
                                              facet = TRUE,
                                              wrap.width = 60,
                                              bar.label = TRUE,
                                              #bar.label.size = 1,
                                              sharings = c(1, 2, 3, 4, 5, 6, 7))
bp_Pathway_smaller_clusters + theme(axis.text.y=element_text(size=10), plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=rep("gray", 6))


saveRDS(list_eTerm_clust, file = paste0(inpath, "clust_var_XGR_results_REACTOME_sh50.rds"))
saveRDS(list_consise_eTerm_clust, file = paste0(inpath, "clust_var_XGR_concise_results_REACTOME_sh50.rds"))

########################### Plot ###########################

# Cluster trajectory plots 
n_genes_per_cluster <- cluster_genes %>%
  group_by(cluster) %>% 
  summarize(n = n()) 

cluster_variance_plt_df <- cluster_variance %>%
  # add column with cluster name
  dplyr::filter(!is.na(cluster)) %>%
  left_join(n_genes_per_cluster, by = "cluster") %>%
  mutate(cluster = as.factor(paste0("Cluster v", cluster, ", n genes = ", n)),
         timepoint = gsub(" weeks", "", timepoint)) %>%
  group_by(gene) %>%
  mutate(variance_scaled = as.numeric(scale(variance)))

# mean variance per time-point per cluster
mean_per_clust <- cluster_variance_plt_df %>%
  # add column with cluster name
  group_by(cluster, timepoint) %>% 
  summarize(cluster_tp_mean = mean(variance_scaled)) 

# Plot cluster trajectories (Main Figure 3B)
cluster_variance_plt_df %>%
  ggplot(aes(x = timepoint, y = variance_scaled, group = gene)) + 
  geom_line(color = "#B07BB9", alpha = 1, linewidth = 0.05) + 
  geom_line(
    data = mean_per_clust, 
    aes(x = timepoint, y = cluster_tp_mean, group = cluster), 
    color = "black", 
    #linetype = "dashed",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) + 
  geom_point(
    data = mean_per_clust, 
    aes(x = timepoint, y = cluster_tp_mean, group = cluster), 
    color = "black"
  ) + 
  facet_wrap(~cluster, nrow = 2) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)) + 
  labs(#title = "Scaled variance for genes clustered by change in variance over time",
    x = "Sample time-point (weeks)",
    y = "Scaled variance") + 
  scale_x_discrete(expand = c(0.05, 0))


# Plot pathway enrichment 

list_consise_eTerm <- readRDS(paste0(inpath, "clust_var_XGR_concise_results_REACTOME_sh50.rds"))

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
  d <- df_all[ind, c("id", "name", "fc", "adjp", "zscore", "pvalue", "group")]
  d$group <- factor(d$group, levels = rev(list_names))
  d <- d %>% dplyr::mutate(direction = ifelse(zscore > 0, 1, -1))
  d$name <- factor(d$name, levels = unique(d$name))
  return(d)
}

names(list_consise_eTerm) <- paste0("Cluster v", 1:6)
df_dotplot <- get_xEnrichCompare_data(list_consise_eTerm, FDR.cutoff = 0.05)

# Get original order BEFORE modification
original_order <- levels(df_dotplot$name)

short_names <- c(
  "Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template" = 
    "Y-family polymerase translesion synthesis bypasses DNA lesions",
  "Formation of the ternary complex, and subsequently, the 43S complex" = "Ternary and 43S complex formation",
  "Downstream signaling events of B Cell Receptor (BCR)" = "BCR signaling", 
  "APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1" = "APC/C:Cdh1 degrades Cdc20 in late mitosis/G1"
)

# Replace levels using the named vector
levels(df_dotplot$name) <- ifelse(
  levels(df_dotplot$name) %in% names(short_names),
  short_names[levels(df_dotplot$name)],
  levels(df_dotplot$name)
)

name_order <- levels(df_dotplot$name)
wrapped_names <- name_order#stringr::str_wrap(name_order, width = 50)
name_map <- setNames(wrapped_names, name_order)
df_dotplot <- df_dotplot %>%
  mutate(
    name_wrapped = name_map[as.character(name)],
    name_wrapped = factor(name_wrapped, levels = wrapped_names)
  )

# Add names of clusters that weren't enriched for any pathways
new_group <- "Cluster v5"
new_group_2 <- "Cluster v6"

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
  dplyr::filter(!is.na(fc)) %>%
  dplyr::count(name_wrapped, name = "n_groups") %>%
  arrange(desc(n_groups)) %>%
  pull(name_wrapped)

# Pathway enrichment plot (Main Figure 3C)
df_dotplot_extended %>% 
  mutate(
    name_wrapped = factor(name_wrapped, levels = pathway_order_wrapped), 
    log2fc = log2(fc)) %>%
  ggplot(aes(x = name_wrapped, y = group)) +
  geom_point(aes(size = -log10(adjp), fill = log2fc), color = "black", shape = 21) +
  scale_size_continuous(
    breaks = c(2, 5, 10)
  ) + 
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    na.value = "gray",
    #limits = color_limits,
    name = "log2 fold\nenrichment"
  ) +
  #coord_flip() +
  labs(
    size = "-log10(FDR)",
    x = "Pathway",
    y = "Variance cluster"
  ) +
  #ggtitle("Change in variance cluster pawthway enrichment") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
    axis.text.y = element_text(color="black", size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.title = element_text(size = 9),
    legend.position = "bottom"
  )  +
  geom_vline( xintercept = 2 + 0.5 ) + 
  geom_vline( xintercept = 40 + 0.5 ) + 
  geom_vline( xintercept = 55 + 0.5 ) + 
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
    x = 21.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Cluster v1",
    size = 3
  ) +
  annotate(
    "text",
    x = 48.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Cluster v2",
    size = 3
  ) +
  annotate(
    "text",
    x = 57,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Cluster v3",
    size = 3
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  # add extra space above
  )


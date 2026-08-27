# 03_3_deg_pathway_enrichment.R 

################################################################################

# 3.3. Pathway enrichment analysis 

################################################################################

# Aim: identify pathways enriched among differentially expressed genes 

########################### Output paths ##########################
xgr_outpath <- "rna-seq/analysis/dream/outputs/XGR_results_REACTOME.rds"
xgr_concise_outpath <- "rna-seq/analysis/dream/outputs/XGR_concise_results_REACTOME.rds"

########################### Input paths ###########################
# top table
toptable_inpath <- "rna-seq/analysis/dream/outputs/pops_toptable.rds"

########################### Parameters ############################
# define deg cutoffs
p_val = 0.01
log_FC = log2(1.5)
# Data location for XGR 
RData.location="/software/team282/XGR/bigdata"

########################### Load packages ###########################
library(XGR)

########################### Load data ###########################
toptable <- readRDS(toptable_inpath)

########################### Analysis ###########################
list_eTerm = list()
list_consise_eTerm <- list()

comps <- unique(toptable$comp)
for(i in 1:length(comps)){
  print(paste(comps[i], "Start"))
  data_table <- toptable %>% filter(comp == comps[i])
  background <- data_table$gtf.gene_name
  data <- (data_table %>% filter(adj.P.Val < p_val, abs(logFC) > log_FC))$gtf.gene_name
  GO <- xEnricherGenes(data=data, background=background, ontology="REACTOME", ontology.algorithm="none", 
                       RData.location=RData.location) #https://xgr.r-forge.r-project.org/
  list_eTerm[[i]] <- GO
  if(nrow(GO$term_info) > 1){
    list_consise_eTerm[[i]] <- xEnrichConciser(GO)
  }else{
    list_consise_eTerm[[i]] <- GO
  }
  print(paste(comps[i], "End"))
}

names(list_eTerm) <- (comps)
names(list_consise_eTerm) <- (comps)

# write this out as an r object
saveRDS(list_eTerm, file = xgr_outpath)
saveRDS(list_consise_eTerm, xgr_concise_outpath)

########################### Plot ###########################
# read in data
list_consise_eTerm <- readRDS(xgr_concise_outpath)

# Plot with XGR default
bp_Pathway_smaller <- xEnrichCompare(list_consise_eTerm[rev(comps)], #rev(comps)
                                     displayBy="fc",
                                     FDR.cutoff=5e-2,
                                     facet = TRUE,
                                     wrap.width = 60,
                                     bar.label = TRUE,
                                     #bar.label.size = 1,
                                     sharings = c(1, 2, 3, 4, 5, 6))
bp_Pathway_smaller + theme(axis.text.y=element_text(size=10), plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values=rep("gray", 6))


# Plot in format for paper (Main figure 1E)
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

df_dotplot <- get_xEnrichCompare_data(list_consise_eTerm, FDR.cutoff = 0.05) 

# Get original order BEFORE modification
original_order <- levels(df_dotplot$name)

levels(df_dotplot$name) <- ifelse(
  levels(df_dotplot$name) == 
    "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)",
  "Regulation of IGF transport and uptake by IGF binding proteins",
  levels(df_dotplot$name)
)

name_order <- levels(df_dotplot$name)
wrapped_names <- stringr::str_wrap(name_order, width = 50)
name_map <- setNames(wrapped_names, name_order)
df_dotplot <- df_dotplot %>%
  mutate(
    name_wrapped = name_map[as.character(name)],
    name_wrapped = factor(name_wrapped, levels = wrapped_names)
  )


# Include time-point contrast that doesn't have any enriched genes
new_group <- "28v36"

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

# Bind to original df_dotplot
df_dotplot_extended <- dplyr::bind_rows(df_dotplot, empty_rows) %>%
  mutate(# Map original name to wrapped label (use as.character to avoid factor issues)
    name_wrapped = name_map[as.character(name)],
    # Convert to factor with wrapped levels to preserve order
    name_wrapped = factor(name_wrapped, levels = wrapped_names))

pathway_order_wrapped <- df_dotplot_extended %>%
  filter(!is.na(fc)) %>%
  count(name_wrapped, name = "n_groups") %>%
  arrange(desc(n_groups)) %>%
  pull(name_wrapped)

annot_df <- tibble(
  label = c("Metabolism", "Immune", "Signaling"),
  x = c(5, 15, 25),   # midpoint of each chunk (numeric positions)
  y = Inf
)

df_dotplot_extended %>% 
  mutate(#group = factor(group, levels = c("12v20", "12v28", "12v36", "20v28", "20v36", "28v36")),
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
    y = "Time-point contrast"
  ) +
  #ggtitle("DEG pathway enrichment") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.title = element_text(size = 9)#,
    #legend.position = "top"
  ) + 
  geom_vline( xintercept = 10 + 0.5 ) + 
  geom_vline( xintercept = 15 + 0.5 ) + 
  geom_vline( xintercept = 27 + 0.5 ) + 
  annotate(
    "text",
    x = 5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "Shared",
    size = 3
  ) +
  annotate(
    "text",
    x = 13,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "12v20",
    size = 3
  ) +
  annotate(
    "text",
    x = 21.5,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "12v36",
    size = 3
  ) +
  annotate(
    "text",
    x = 29,          # horizontal position (in data units)
    y = Inf,         # places the text at the top edge
    vjust = -1,      # pushes it ABOVE the plot panel
    label = "20v36",
    size = 3
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  # add extra space above
  )


# 04_4_WGCNA_module_annotation.R 

################################################################################

# 4.4. Annotating WGCNA modules

################################################################################

# Aim: Identify the biological function of each WGCNA module based on its genes

########################### Output paths ##########################
outfiles_path <- "rna-seq/analysis/WGCNA/outputs/"

########################### Input paths ###########################
# sample info 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
# eigengenes
eigengenes_inpath <- "rna-seq/analysis/WGCNA/outputs/eigengenes_consensus_signed_hybrid.csv"
# modules
modules_inpath <- "rna-seq/analysis/WGCNA/outputs/consensus_modules_signed_hybrid.csv"
# connectivity
connectivity_inpath <- "rna-seq/analysis/WGCNA/outputs/connectivity_consensus_signed_hybrid.csv"
# cibersort output 
cibersort_inpath <- "rna-seq/analysis/cibersort/outputs/CIBERSORTx_Adjusted.txt"
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Parameters ############################
# Data location for XGR 
RData.location="/software/team282/XGR/bigdata"

########################### Load packages ###########################
library(tidyverse)
library(XGR)
library(xCell)
library(decoupleR)

########################### Load data ###########################
# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
# read in metadata  
sample.info <- read.csv(sample_info_inpath)
# eigengenes
eigengenes <- read.csv(eigengenes_inpath, row.names = 1)
# modules
modules <- read.csv(modules_inpath, row.names = 1)
# connectivity
connectivity <- read.csv(connectivity_inpath, row.names = 1)
# cibersort
cibersort <- read.delim(cibersort_inpath) %>% 
  column_to_rownames(var = "Mixture") %>% select(-c(P.value, Correlation, RMSE))

########################### Analysis ###########################

# Prepare gene ID key 
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

################################################################################
# 1. Module pathway enrichment - Reactome
################################################################################
list_eTerm = list()
list_consise_eTerm = list()
# list background genes which will be the same for each module
gene_ids_using <- rownames(modules)
background <- id_name %>% filter(gtf.gene_id %in% gene_ids_using) %>% pull(gtf.gene_name)
Modules <- unique(modules$Module)[1:(length(unique(modules$Module))-1)]
for(i in 1:length(Modules)){
  print(paste(Modules[i], "Start"))
  # convert these to gene NAMES
  data <- (modules %>% filter(Module == Modules[i])  %>%
             rownames_to_column(var = "Gene") %>% 
             left_join(id_name, by = c("Gene" = "gtf.gene_id")))$gtf.gene_name
  GO <- xEnricherGenes(data=data, background=background, ontology="REACTOME", ontology.algorithm="none", 
                       RData.location=RData.location)
  list_eTerm[[Modules[i]]] <- GO
  if(nrow(GO$term_info) > 1){
    list_consise_eTerm[[Modules[i]]] <- xEnrichConciser(GO)
  }else{
    list_consise_eTerm[[Modules[i]]] <- GO
  }
  print(paste(Modules[i], "End"))
}
length(list_eTerm)
names(list_eTerm)
length(list_consise_eTerm)
names(list_consise_eTerm)

# write this out as an r object
saveRDS(list_eTerm, file = paste0(outfiles_path,"module_annotation/module_XGR_results_reactome.rds"))
saveRDS(list_consise_eTerm, file = paste0(outfiles_path,"module_annotation/module_XGR_results_concise_reactome.rds"))
# read it back in 
list_concise_eTerm_reactome <- readRDS(paste0(outfiles_path,"module_annotation/module_XGR_results_concise_reactome.rds"))
# make a df
module_pathway_annot_reactome_concise <- c()
for(i in 1:length(list_concise_eTerm_reactome)){
  pathways <- xEnrichViewer(list_concise_eTerm_reactome[[i]], 
                            top_num = length(list_concise_eTerm_reactome[[i]]$annotation)) %>% 
    filter(adjp < 0.05) %>% 
    arrange(desc(fc)) %>% 
    select(name, adjp, fc) %>% 
    mutate(Module = names(list_concise_eTerm_reactome)[[i]])
  module_pathway_annot_reactome_concise <- rbind(module_pathway_annot_reactome_concise, pathways)
}
# Write this out 
module_pathway_annot_reactome_concise %>% 
  write.csv(file = paste0(outfiles_path,"module_annotation/module_XGR_results_concise_signif_reactome.txt"), row.names = FALSE)

################################################################################
# Module association with cibersort cell proportions
################################################################################
# Associate modules with cell proportions only based on module eigengene values at time point 1 
# remove cell types with 0 variation 
cells_0_variance <- names(which(apply(X = cibersort, MARGIN = 2, FUN = var) == 0))
cibersort_cells_var <- cibersort %>% select(-cells_0_variance)
cibersort.cells <- colnames(cibersort_cells_var)
miscarried_samples_to_remove <- c("pops212970284", "pops212971619")

# subset eigengenes to time point 1 
eigengenes_tp1 <- eigengenes %>% 
  rownames_to_column(var = "RNA_sanger_sample_id") %>%
  left_join(sample.info %>% select(RNA_sanger_sample_id, Sample_taken_at)) %>%
  # remove the possible double samples
  filter(Sample_taken_at == "12 weeks",
         !(RNA_sanger_sample_id %in% miscarried_samples_to_remove)) %>%
  column_to_rownames(var = "RNA_sanger_sample_id") %>%
  select(-Sample_taken_at)
# for each eigengene 
estimates <- lapply(colnames(eigengenes_tp1), function(eigengene.name) {
  # for each cell type 
  lapply(cibersort.cells, function(cell.type) {
    # test the correlation between the the cell type and the eigengene values
    test.result = cor.test(cibersort_cells_var[, cell.type], eigengenes[, eigengene.name], alternative="two.sided", method="spearman", exact=FALSE)
    return(c(test.result$estimate, test.result$p.value, cell.type))
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(Eigengene=eigengene.name)
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::select(Rho=1, P.Value=2, Association.Variable=3, Eigengene) %>%
  dplyr::mutate(Rho=as.numeric(Rho), P.Value=as.numeric(P.Value)) %>%
  dplyr::arrange(desc(abs(Rho))) %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="BH")) %>%
  dplyr::mutate(Association.Variable.Type="CIBERSORTx Cell Proportion", Statistic.Type="Rho") %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Rho, Statistic.Type, P.Value, Adjusted.P.Value)

# write this out 
write.csv(estimates, paste0(outfiles_path, "module_annotation/estimates.eigengene.cibersort.cell.proportion.association.csv"))
# Plot so you can see where to make cutoff for how to filter for cell props
# Cutoff at +0.3
cor_cel_prop <- read.csv(paste0(outfiles_path, "module_annotation/estimates.eigengene.cibersort.cell.proportion.association.csv"), row.names = 1)
cor_cel_prop %>% select(Eigengene, Statistic, Adjusted.P.Value, Association.Variable) %>%
  ggplot(aes(x = Eigengene, y = Statistic, color = Adjusted.P.Value < 0.05)) + 
  geom_point() + geom_hline(yintercept = 0.3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

################################################################################
# xCell cell type enrichment 
################################################################################
# Prep gene annotations + count matrix 
# get list of gene IDs in our count df 
gene.map <- data.frame(Ensembl.ID=rownames(modules)) %>% 
  # joun to gene names 
  left_join(id_name, by = c("Ensembl.ID" = "gtf.gene_id")) %>%
  #rename("gene_name" = "gtf.gene_name") %>%
  column_to_rownames(var = "Ensembl.ID")
colnames(gene.map)[1] <- "gene_name"

# Remove genes with duplicate gene names
duplicated.genes <- gene.map$gene_name[duplicated(gene.map$gene_name)]
gene.map <- gene.map %>%
  filter(!(gene_name %in% duplicated.genes))

# Overrepresentation Analysis with xCell Signatures
# Reduce Modules to Genes with Gene Names
modules.list <- modules %>%
  rownames_to_column(var = "Gene") %>% 
  dplyr::filter(Gene %in% rownames(gene.map)) %>%
  left_join(gene.map %>% rownames_to_column(var = "Gene")) %>%
  dplyr::select(Gene=gene_name, Module) %>%
  # make into a list 
  split(., .$Module)

# Perform Overrepresentation Test

# get cell unique cell types from xCell 
signature.cell.types <- sapply(strsplit(names(xCell.data$signatures), "%"), function(x) { x[1] })
unique.cell.types <- unique(signature.cell.types)
# subset these to overlapping with LM22 cell types for consistency with the previous analysis 
LM22_celltypes <- colnames(cibersort)
unique.cell.types_LM22 <- c("Memory B-cells", "naive B-cells", "B-cells", "Plasma cells", "CD8+ T-cells", "CD4+ naive T-cells", "CD4+ memory T-cells", "Tregs", "Tgd cells", "NK cells", "Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2", "DC", "aDC", "cDC", "Mast cells", "Eosinophils", "Neutrophils")

# get list of marker genes for each cell type
gene.lists <- list()
for (cell.type in unique.cell.types_LM22) {
  n.signatures = sum(signature.cell.types == cell.type)
  gene.lists[[cell.type]] = lapply(1:n.signatures, function(i) {
    xCell.data$signatures[signature.cell.types == cell.type][[i]]@geneIds
  })
}

# set up output df for overrepresentation analysis results
ora.results <- list()
ora.results[["Module"]] <- list()
ora.results[["Cell.Type"]] <- list()
ora.results[["Overlap"]] <- list()
ora.results[["N.Gene.Set"]] <- list()
ora.results[["N.Module"]] <- list()
ora.results[["N.Not.Gene.Set"]] <- list()
ora.results[["P.Value"]] <- list()
ora.results[["xCell.Genes"]] <- list()
ora.results[["Module.Overlap.Genes"]] <- list()

counter <- 0

# for each module, for each cell type, for each list of marker genes, use phyper to do overrepresentation analysis 
for (module in names(modules.list)) {
  
  for (cell.type in unique.cell.types_LM22) {
    
    for (i in 1:length(gene.lists[[cell.type]])) {
      
      counter <- counter + 1
      
      gene.list = gene.lists[[cell.type]][[i]]
      module.list = modules.list[[module]]$Gene
      
      within.module.and.gene.set = length(intersect(gene.list, module.list))
      within.gene.set = length(intersect(gene.list, gene.map$gene_name))
      not.within.gene.set = length(setdiff(gene.map$gene_name, gene.list))
      within.module = length(module.list)
      
      p.value = phyper(
        within.module.and.gene.set - 1,
        within.gene.set,
        not.within.gene.set,
        within.module,
        lower.tail=F
      )
      
      ora.results[["Module"]][[counter]] <- module
      ora.results[["Cell.Type"]][[counter]] <- cell.type
      ora.results[["Overlap"]][[counter]] <- within.module.and.gene.set
      ora.results[["N.Gene.Set"]][[counter]] <- within.gene.set
      ora.results[["N.Module"]][[counter]] <- within.module
      ora.results[["N.Not.Gene.Set"]][[counter]] <- not.within.gene.set
      ora.results[["P.Value"]][[counter]] <- p.value
      ora.results[["xCell.Genes"]][[counter]] <- paste0(gene.list, collapse="|")
      ora.results[["Module.Overlap.Genes"]][[counter]] <- paste0(intersect(gene.list, module.list), collapse="|")
    }
  }
}

# make the ora results into a df 
ora.results.df <- as.data.frame(lapply(ora.results, unlist))

# update results 
filtered.cell.type.enrichment <- ora.results.df %>%
  # adjust p vals for multiple testing correction 
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="BH")) %>%
  # filter to p-vals that are significant 
  dplyr::filter(Adjusted.P.Value < 0.05) %>%
  # make enrichment score 
  dplyr::mutate(Enrichment=(Overlap / N.Module) / (N.Gene.Set / N.Not.Gene.Set)) %>%
  # group by module and cell type (across the different marker gene sets per cell type)
  dplyr::group_by(Module, Cell.Type) %>%
  # get mean, median, min, max enrichment scores
  dplyr::summarize(
    Mean.Enrichment=mean(Enrichment), Median.Enrichment=median(Enrichment),
    Enrichment.Min=min(Enrichment), Enrichment.Max=max(Enrichment),
    .groups="drop"
  ) %>%
  dplyr::arrange(Module, desc(Median.Enrichment))

# write this out 
write.csv(filtered.cell.type.enrichment, paste0(outfiles_path,"module_annotation/module.xcell.signature.enrichment.csv"))

################################################################################
# Regulon enrichment with dorothea
################################################################################
# Code from https://github.com/davenportlab/sepsis_eqtl/blob/main/3.modules/3.modQTL_TFs.R
# Dorothea regulons
net <- decoupleR::get_dorothea(organism='human', levels=c('A', 'B', 'C'))
# get all genes for background
genes.using <- id_name %>% filter(gtf.gene_id %in% rownames(modules))
# limit regulons to measured genes
net <- subset(net, target %in% genes.using$gtf.gene_name)
genes.using$regulon <- genes.using$gtf.gene_name %in% net$target

regulon.enrichment <- list()

# Test each module for enrichment of each regulon
for(i in 1:(length(unique(modules$Module))-1)){
  module <- paste0("Module_", i)
  print(module)
  genes <- subset(modules %>% rownames_to_column(var = "Gene") %>% left_join(id_name, by = c("Gene" = "gtf.gene_id")), Module == module)
  
  regulon.enrichment.module <- data.frame(Module=module,
                                          TF=unique(net$source),
                                          OR=NA,
                                          Pval=NA)
  
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
          regulon.enrichment.module$OR[r] <- results$estimate
          regulon.enrichment.module$Pval[r] <- results$p.value
        }
      }
    }
  }
  regulon.enrichment.module <- regulon.enrichment.module[complete.cases(regulon.enrichment.module),]
  regulon.enrichment.module$FDR <- p.adjust(regulon.enrichment.module$Pval, method="fdr")
  regulon.enrichment[[module]] <- regulon.enrichment.module
}

enrichment.results <- data.table::rbindlist(regulon.enrichment)
enrichment.results$FDRall <- p.adjust(enrichment.results$Pval, method="fdr")
sig.results <- subset(enrichment.results, FDRall<0.05)

# write these out 
sig.results %>% write.csv(file = paste0(outfiles_path, "outputs/module_annotation/regulons_enriched_per_module.csv"), row.names = FALSE)

################################################################################
# Identify hub genes for each module
################################################################################
# Plot rank of connecticity score for each gene within a module to find elbow of cutoff
# Not a clear one, but 10 hub genes looks good enough 
connectivity %>% group_by(Module) %>% 
  mutate(within_mod_rank = rank(-kWithin)) %>% 
  ggplot(aes(x = within_mod_rank, y = kWithin, color = within_mod_rank < 10)) + 
  geom_point() + 
  facet_wrap(~Module, scales = "free") + 
  theme(strip.background = element_blank()) + 
  ggtitle("Ranked within-module connectivity for genes in module")

################################################################################
# Make a table of module annotation information 
################################################################################
# n genes per mod 
n_genes_per_mod = as.data.frame(table(modules$Module))
colnames(n_genes_per_mod) = c("Module", "Number_of_genes")
# top 10 hub genes per mod (as a string)
top_hub_genes_per_module_cat <- connectivity %>%
  rownames_to_column(var = "gtf.gene_id") %>% 
  left_join(id_name) %>%
  # find top kWithin for each module 
  group_by(Module) %>%
  top_n(n = 10, kWithin) %>%
  select(Module, gtf.gene_name) %>%
  group_by(Module) %>%
  mutate(top_10_hub_genes = paste(gtf.gene_name, collapse = ", ")) %>%
  select(-gtf.gene_name) %>% distinct() 
# pathways (all signif reactome)
module_pathway_annot_reactome_concise_cat <- read.csv(paste0(outfiles_path,"module_annotation/module_XGR_results_concise_signif_reactome.txt")) %>%
  select(Module, name) %>%
  group_by(Module) %>%
  mutate(all_signif_enriched_reactome_pathways = paste(name, collapse = ", ")) %>% 
  select(-name) %>% distinct()
# regulon enrichment per module (all signif)
enriched_regulons_per_mod <- read.csv(paste0(outfiles_path, "module_annotation/regulons_enriched_per_module.csv")) %>%
  select(Module, TF) %>%
  group_by(Module) %>%
  mutate(enriched_regulons = paste(TF, collapse = ", ")) %>%
  select(-TF) %>% distinct()
# cell enrichment per module (all signif)
enriched_celltypes_per_mod <- read.csv(paste0(outfiles_path,"module_annotation/module.xcell.signature.enrichment.csv")) %>% 
  select(Module, Cell.Type) %>%
  group_by(Module) %>%
  mutate(enriched_celltypes = paste(Cell.Type, collapse = ", ")) %>%
  select(-Cell.Type) %>% distinct()
# top correlated cell type per module
correlated_celltypes_per_mod <- read.csv(paste0(outfiles_path, "module_annotation/estimates.eigengene.cibersort.cell.proportion.association.csv"), row.names = 1) %>%
  # filter to significant and strong positive correlation (> 0.5)
  filter(Adjusted.P.Value < 0.05, 
         Statistic > 0.3) %>%
  select(Eigengene, Association.Variable) %>%
  group_by(Eigengene) %>% 
  mutate(correlated_celltypes = paste(Association.Variable, collapse = ", ")) %>%
  select(-Association.Variable) %>% distinct() %>%
  #top_n(n = 1, Adjusted.P.Value) %>%
  ungroup() %>%
  mutate(Module = paste0("Module_", gsub("ME_", "", Eigengene))) %>%
  select(-Eigengene)

# join them together 
module_annotation_df <- n_genes_per_mod %>%
  filter(Module != "Unassigned") %>%
  # join to hub genes
  left_join(top_hub_genes_per_module_cat, by = "Module") %>%
  # join to enriched pathways 
  left_join(module_pathway_annot_reactome_concise_cat, by = "Module") %>%
  mutate(all_signif_enriched_reactome_pathways = ifelse(is.na(all_signif_enriched_reactome_pathways), "None", all_signif_enriched_reactome_pathways)) %>%
  # join to enriched cell types per module
  left_join(enriched_celltypes_per_mod, by = "Module") %>%
  mutate(enriched_celltypes = ifelse(is.na(enriched_celltypes), "None", enriched_celltypes)) %>%
  # join to top correlated celltype per mod
  left_join(correlated_celltypes_per_mod, by = "Module") %>%
  mutate(correlated_celltypes = ifelse(is.na(correlated_celltypes), "None", correlated_celltypes)) %>%
  left_join(enriched_regulons_per_mod, by = "Module") %>%
  mutate(enriched_regulons = ifelse(is.na(enriched_regulons), "None", enriched_regulons))
# write this out (Supplementary Table)
module_annotation_df %>% write.table(file = paste0(outfiles_path, "module_annotation/module_annotation_df.txt"), row.names = FALSE, sep = "\t")
# test reading it back in 
mod_annot_df_in <- read.delim(paste0(outfiles_path, "module_annotation/module_annotation_df.txt"), sep = "\t")

# Based on this integrated evidence, modules were assigned the following functions: 
# innate immune response (m1)
# cell proliferation (m2) cell motility (m3)
# T cell immune response (m4)
# protein synthesis and translation (m5, m9)
# B cell immune response (m6)
# eosinophil allergic-like response (m7)
# cytotoxic immune response (m8)
# coagulation and clotting (m10) 
# histones (m11)
# mitochondria and energy production (m12)
# small molecule transport (m13)
# antiviral interferon response (m14) 
# inflammatory response regulation (m15)

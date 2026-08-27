# 04_2_make_consensus_WGCNA_modules.R 

################################################################################

# 4.3. Make WGCNA consensus modules 

################################################################################

# Aim: Create WGCNA consensus modules from separate TOM matrices
# Note: This WGCNA analysis was guided majorly by Nikhil Milind's WGCNA analysis. 
#       https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/04_Expression/scripts/wgcna_01_gene_expression.ipynb
#       Consensus module creation was guided by consensus module WGCNA tutorial. 
#       https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?e=4&preview=Consensus-NetworkConstruction-man.pdf&rlkey=3za224679kv8ulbpl5s9mwrkt&dl=0
#       Longitudinal modeling was guided by Jack Gisby's WGCNA analysis.
#       https://github.com/jackgisby/covid-longitudinal-multi-omics/blob/v1.0/notebooks/4_module_analysis.md

########################### Output paths ##########################
outfiles_path <- "rna-seq/analysis/WGCNA/outputs/"

########################### Input paths ###########################
infiles_path <- "rna-seq/analysis/WGCNA/outputs/sep_preprocessing/"

########################### Parameters ############################

########################### Load packages ###########################
library(WGCNA)
library(tidyverse)

# set ggplot theme 
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             text=element_text(size=15)))

########################### Load data ###########################
# read in tom matricies 
TOM_tp1 <- readRDS(paste0(outfiles_path, "sep_preprocessing/tom.mtx_st13_tp", 1, ".rds"))
TOM_tp2 <- readRDS(paste0(outfiles_path, "sep_preprocessing/tom.mtx_st13_tp", 2, ".rds"))
TOM_tp3 <- readRDS(paste0(outfiles_path, "sep_preprocessing/tom.mtx_st13_tp", 3, ".rds"))
TOM_tp4 <- readRDS(paste0(outfiles_path, "sep_preprocessing/tom.mtx_st13_tp", 4, ".rds"))

########################### Analysis ###########################

# Scale to 95th percentile: This ensures that the scaling adjusts for systematic 
#                           differences while maintaining the relative relationships 
#                           within each TOM. For computational efficiency, 
#                           sample a sufficient number of TOM entries from 
#                           each dataset to estimate the percentile.
scaleP <- 0.95  # Reference percentile
set.seed(12345)  # Reproducibility
# Determine the number of samples needed for accurate percentile estimation
nGenes <- nrow(TOM_tp1)  # Number of genes (rows and columns of the TOMs)
nSamples <- as.integer(1 / (1 - scaleP) * 1000)
# Randomly sample from the upper triangle of each TOM (excluding diagonal)
scaleSample <- sample(nGenes * (nGenes - 1) / 2, size = nSamples)
TOMScalingSamples <- list()
# Store sampled TOM entries for each dataset
TOMScalingSamples[[1]] <- as.dist(TOM_tp1)[scaleSample]
TOMScalingSamples[[2]] <- as.dist(TOM_tp2)[scaleSample]
TOMScalingSamples[[3]] <- as.dist(TOM_tp3)[scaleSample]
TOMScalingSamples[[4]] <- as.dist(TOM_tp4)[scaleSample]
# calculate the 95th percentile
scaleQuant <- sapply(TOMScalingSamples, function(sample) {
  quantile(sample, probs = scaleP, type = 8)
})
# Scale each TOM so that the selected percentile matches the reference percentile. 
# Use the first dataset as the reference.
scalePowers <- log(scaleQuant[1]) / log(scaleQuant)  # Compute scaling powers

# Scale each TOM (apply scaling to all except the reference TOM)
TOM1_scaled <- TOM_tp1  # Reference TOM, no scaling needed
TOM2_scaled <- TOM_tp2^scalePowers[2]
TOM3_scaled <- TOM_tp3^scalePowers[3]
TOM4_scaled <- TOM_tp4^scalePowers[4]

# The array TOM now contains the scaled TOMs. To see what the scaling achieved, 
# we form a quantile-quantile plot of topological overlaps before and after scaling:
# For plotting, also scale the sampled TOM entries
scaledTOMSamples_WGCNA = list()
scaledTOMSamples_WGCNA[[1]] = TOMScalingSamples[[1]]^scalePowers[1]
for(set in 2:4){
  scaledTOMSamples_WGCNA[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
  # qq plot of unscaled samples 
  qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[set]], plot.it = TRUE, cex = 0.6,
                      xlab=paste("TOM in time point 1"), ylab=paste("TOM in time point", set), main = paste0("Q-Q plot of TOM tp1 vs tp", set), pch = 20)
  # qq plot of scaled samples 
  qqScaled = qqplot(scaledTOMSamples_WGCNA[[1]], scaledTOMSamples_WGCNA[[set]], plot.it = FALSE)
  points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
  abline(a = 0, b = 1, col = "blue")
  legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
  dev.off()
}
# All plots saved toegther to outputs qq.plt_all.pdf

# Calculation of consensus Topological Overlap
consensusTOM <- pmin(TOM1_scaled, TOM2_scaled, TOM3_scaled, TOM4_scaled)
# write this out 
saveRDS(consensusTOM, paste0(outfiles_path, "consensus_TOM.rds"))

# Make tom distance matrix
TOM.dist = 1 - consensusTOM

# Clustering TOM Matrix
# The gene co-expression modules are generated by clustering the TOM matrix.
# Make the dendrogram
dendrogram = hclust(as.dist(TOM.dist), method="average")
# Plot dendrogram 
plot(dendrogram, labels=FALSE, main="Gene Expression TOM Dendrogram")

# Do tree cutting 
dynamic.mods = cutreeDynamic(
  dendro=dendrogram, distM=TOM.dist, 
  pamRespectsDendro=FALSE, # parameter from jack and Nikhil 
  method = "hybrid", # default parameter, jack 
  # Adjust for cluster size adjustment
  deepSplit = 4, # 1 by default, can be 0-4. Higher values = more + smaller clusters
  minClusterSize = 50
)
# Show number of modules and their sizes
table(dynamic.mods)

# Read in the expression data
multiExpr = vector(mode="list",length = 4)
for(i in 1:4){
  # Read in one file at a time
  multiExpr[[i]] <- list(data = readRDS(paste0(outfiles_path, "sep_preprocessing/gene.exp_tp", i, ".rds")))
  names(multiExpr[[i]]$data) = paste0("tp", i)
}
# Check this looks good 
exprSize = checkSets(multiExpr)

# merge similar modules based on their module eigengenes. 
# Merging all modules with a distance of less than 0.1 
# (a correlation between module eigengenes of greater than 0.9)
me.data = multiSetMEs(multiExpr, colors = NULL, universalColors = dynamic.mods)
consMEDiss = consensusMEDissimilarity(me.data)
me.tree = hclust(as.dist(consMEDiss), method ="average")
# Plot the clustered module eigengenes with a dist threshold of 0.1
plot(me.tree, main="Clustering of Module Eigengenes")
me.dist.threshold = 0.1
abline(h=me.dist.threshold, col="firebrick1")
# Merge similar clusters (we don't have any that are below threshold)
# This plot is saved to initial_module_tree.pdf
merge = mergeCloseModules(
  multiExpr,#gene.exp, 
  colors=dynamic.mods, 
  cutHeight=me.dist.threshold,
  corOptions=list(use="p", method="spearman")
)
merge.mods = merge$colors
# The merged modules are not labeled by frequency. 
# Refactor the module names to name them in decreasing order. 
# The 0-th module represents unassigned genes and remains 0 
# regardless of the number of genes in the module.
merge.mods.no.zero = merge.mods[merge.mods != 0]
freq.mod.labels = 0:length(table(merge.mods.no.zero))
names(freq.mod.labels) = c("0", names(table(merge.mods.no.zero))[order(table(merge.mods.no.zero), decreasing=T)])
merge.mods = as.numeric(as.character(plyr::revalue(factor(merge.mods), replace=freq.mod.labels)))
# Show the results of this 
table(merge.mods)

# Plot dendrogram with module colors (Supplementary Figure 25)
mergedColors = labels2colors(dynamic.mods)
plotDendroAndColors(dendrogram, 
                    mergedColors, 
                    "Module colors", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE,
                    guideHang = 0.05)
# This plot is output to colored_dendrogram.pdf

# Intra- and Inter-Module Connectivity
connectivity <- intramodularConnectivity(consensusTOM, colors=merge.mods)
rownames(connectivity) <- colnames(multiExpr[[1]]$data)

# Co-Expression Module Size Distribution (Supplementary Figure 10)
plot.data <- as.data.frame(table(merge.mods)) %>%
  dplyr::select(Module=1, Frequency=Freq) %>%
  dplyr::filter(Module != 0)
mod_size_plt <- ggplot(plot.data) +
  geom_col(aes(x=Module, y=Frequency), fill="black", color="black", width=1, position=position_dodge(1)) +
  xlab("Module") + ylab("Number of Genes") + ggtitle("Number of genes in each module") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
# This plot is output to mod_size.pdf

### Module Assignment
gene.list <- colnames(multiExpr[[1]]$data)
mod.assignment <- merge.mods
mod.labels <- paste0("Module_", mod.assignment)
mod.labels[mod.labels == "Module_0"] <- "Unassigned"

modules <- data.frame(
  Gene=gene.list,
  Module=mod.labels
) %>%
  dplyr::arrange(Module, Gene)
# write it out 
write.csv(modules, paste0(outfiles_path, "consensus_modules_signed_hybrid.csv"), row.names=F)

# load data
# normalized counts 
print("Reading in count data...")
counts <- read.table("outputs/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt", header = TRUE) 
dim(counts)
counts[1:5, 1:5]
print("Count data loaded.")

### Eigengenes
# Calculate eigengenes on all samples
eigengenes.result <- moduleEigengenes(t(counts)[,mod.labels != "Unassigned"], colors=mod.labels[mod.labels != "Unassigned"])
eigengene.order <- order(as.numeric(sapply(strsplit(colnames(eigengenes.result$eigengenes), "_"), function(x) x[2])))
eigengenes <- eigengenes.result$eigengenes
colnames(eigengenes) <- sapply(strsplit(colnames(eigengenes), "_"), function(x) paste0("ME_", x[2]))

eigengenes <- eigengenes[,eigengene.order]
rownames(eigengenes) <- colnames(counts)
dim(eigengenes)
eigengenes[1:5, 1:5]
# Write out eigengenes 
write.csv(eigengenes, paste0(outfiles_path, "eigengenes_consensus_signed_hybrid.csv"))

# variance explained by each eigengene
variance.explained <- data.frame(
  Module=colnames(eigengenes),
  Var.Explained=t(eigengenes.result$varExplained[eigengene.order])[,1]
)
# variance.explained %>% ggplot(aes(x = Var.Explained)) + geom_histogram() + ggtitle("Variance explaiend by eigengenes")
# Write out variance explained
write.csv(variance.explained, paste0(outfiles_path, "var_explained_consensus_signed_hybrid.csv"), row.names = FALSE)

### Connectivity 
connectivity <- connectivity[modules$Gene,]
connectivity <- cbind(connectivity, modules$Module)
colnames(connectivity)[ncol(connectivity)] <- "Module"
# Write out connectivity 
write.csv(connectivity, file = paste0(outfiles_path, "connectivity_consensus_signed_hybrid.csv"))

# Check signs of eigengenes 
# Get module eigengene correlation with gene expression 
module_indexes <- 1:length(unique(modules$Module)[! unique(modules$Module) == "Unassigned"])
eigengene_corrs <- list()
eigengene_corr_sign <- list()
counts_t <- counts %>% 
  t() %>% 
  as.data.frame()
for(i in 1:length(module_indexes)){
  print(i)
  # For each module
  # 1. Find member correlations
  # Get eigengene value
  eigengene_val <- eigengenes %>% 
    select(paste0("ME_", module_indexes[i])) %>% 
    rownames_to_column(var = "RNA_sanger_sample_id")
  # Get gene expression values of genes in this module
  genes_in_module <- (modules %>% filter(Module == paste0("Module_", module_indexes[i])))$Gene
  counts_in_module <- counts_t %>%
    select(all_of(genes_in_module)) %>% 
    rownames_to_column(var = "RNA_sanger_sample_id") 
  member_cors <- cor(eigengene_val[,-1], counts_in_module[,-1])
  eigengene_corrs[[i]] <- member_cors
  names(eigengene_corrs)[[i]] <- paste0("Module_", module_indexes[i])
  eigengene_corr_sign[[i]] <- sign(mean(member_cors))
  names(eigengene_corr_sign)[[i]] <- paste0("Module_", module_indexes[i])
}

eigengene_corr_sign_vec <- unlist(eigengene_corr_sign, use.names=TRUE)

# Change signs - in this case we don't need to switch any - all are positive
eigengenes_signed = sweep(eigengenes, 2, eigengene_corr_sign_vec,`*`) 
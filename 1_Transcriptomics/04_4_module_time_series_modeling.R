# 04_4_module_time_series_modeling.R

################################################################################

# 4.4. Time-series modeling of WGCNA modules

################################################################################

# Aim: Model WGCNA module eigengene values over time

########################### Output paths ##########################
outfiles_path <- "rna-seq/analysis/WGCNA/outputs/"

########################### Input paths ###########################
# sample info 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
# eigengenes
eigengenes_inpath <- "rna-seq/analysis/WGCNA/outputs/eigengenes_consensus_signed_hybrid.csv"
# connectivity
connectivity_inpath <- "rna-seq/analysis/WGCNA/outputs/connectivity_consensus_signed_hybrid.csv"
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"
# module pathway enrichment
pathway_enrichment_inpath <- "rna-seq/analysis/WGCNA/outputs/module_annotation/module_XGR_results_concise_reactome.rds"


########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)
library(rstatix)
library(variancePartition)
library(ggeffects)
library(interactions)
library(dendextend)

# set ggplot theme 
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             text=element_text(size=15)))


########################### Load data ###########################
# make sample info df
sample.info <- read.csv(sample_info_inpath)
# eigengenes
eigengenes <- read.csv(eigengenes_inpath, row.names = 1)
# connectivity
connectivity <- read.csv(connectivity_inpath, row.names = 1)
# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
# read in module pathway enrichment
list_concise_eTerm_reactome <- readRDS(pathway_enrichment_inpath)

########################### Analysis ###########################

################ Module association with time point ############################
# Nikhil used a repeated measures ANOVA (which is an extension of the paired 
#                                        -test) to test for association between eigengenes and time points.
# For this analysis to work, we can only have one sample per time point per person. 
# We have two individuals who miscarried and re-enrolled. 
# SO we will remove pops212970284 and pops212971619 from the following analysis
test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["GES"]] <- list()
miscarried_samples_to_remove <- c("pops212970284", "pops212971619")

# merge sample info to eigengene info
eigen.time.data <- sample.info %>% left_join(eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id"), by = "RNA_sanger_sample_id") %>%
  mutate(ANON_ID = as.factor(ANON_ID), 
         Sample_taken_at = as.factor(Sample_taken_at)) %>%
  filter(!(RNA_sanger_sample_id %in% miscarried_samples_to_remove))

for (eigengene in colnames(eigengenes)) {
  
  res <- eigen.time.data %>%
    anova_test(dv=lazyeval::interp(eigengene), wid=ANON_ID, within=Sample_taken_at) %>%
    get_anova_table()
  
  test.results[["P.Value"]][[eigengene]] <- res$p
  test.results[["GES"]][[eigengene]] <- res$ges
}

associated.eigengenes <- lapply(test.results, unlist) %>%
  as.data.frame() %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="BH")) %>%
  dplyr::arrange(desc(GES)) %>%
  dplyr::mutate(Association.Variable="D1/D3/D5", Association.Variable.Type="Time Point", Statistic.Type="Generalized Effect Size") %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=GES, Statistic.Type, P.Value, Adjusted.P.Value)

print(paste("Number of modules associated with eigengenes is", associated.eigengenes %>% filter(Adjusted.P.Value < 0.05) %>% nrow(), "out of", length(colnames(eigengenes))))

write.csv(associated.eigengenes, paste0(outfiles_path, "estimates.eigengene.time.point.association.csv"))

# plot this 
plot_data <- eigengenes %>% 
  rownames_to_column(var = "RNA_sanger_sample_id") %>%
  # join to discrete time point
  left_join(sample.info %>%
              select(RNA_sanger_sample_id, Sample_taken_at, SampleGA), by = "RNA_sanger_sample_id") %>%
  pivot_longer(cols = -c(RNA_sanger_sample_id, Sample_taken_at, SampleGA), names_to = "Module", values_to = "Eigengene_expression") %>%
  mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module))), 
         `Sample time-point (weeks)` = gsub("_weeks", "", Sample_taken_at)) 

module.pval.annotation <- associated.eigengenes %>% 
  mutate(round_adj.p = paste0("p=",signif(Adjusted.P.Value, digits=3))) %>%
  select(round_adj.p) %>% rownames_to_column(var = "Module_name") %>%
  mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module_name)))) %>%
  left_join(plot_data %>% select(`Sample time-point (weeks)`, Module, Eigengene_expression)) %>%
  group_by(Module) %>% 
  mutate(max_val = max(Eigengene_expression)) %>%
  top_n(1, Eigengene_expression) %>%
  select(-Module_name)
# Plot module eigengene expression over time (Supplementary Figure 26)
plt_discrete_tp_long <- plot_data %>%
  ggplot(aes(x = `Sample time-point (weeks)`, 
             y = Eigengene_expression 
             #fill = `Sample time-point (weeks)`, 
             )) + 
  geom_jitter(aes(color = `Sample time-point (weeks)`), 
              size = 0.1, 
              width = 0.3) + 
  facet_wrap(~Module, scales = "free", ncol = 3) + 
  ggtitle("Signed hybrid consensus module eigengene values by time-point") + 
  geom_boxplot(fill = "transparent", outlier.alpha = 0) + ylab("Module eigengene value") + 
  #scale_fill_manual(values=c("#BCE4D8", "#83C4CB", "#439FB7", "#32769B")) +
  scale_color_manual(values=c("#BCE4D8", "#83C4CB", "#439FB7", "#32769B")) + 
  theme(legend.position = "none",
        strip.background = element_blank()) + 
  geom_text(data = module.pval.annotation, aes(label = round_adj.p, x = 1, y = max_val+.01))


################ Continuous temporal modeling of modules ############################

# Mixed effects modeling of continuous time point (GA)
# Run variance partition to see which covariates are contributing to the module eigengene values 
# Transpose eigengenes for input into variance partition
var_part_eig <- eigengenes %>% 
  #select(-sample_id) %>% 
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(.)))
# make variance partition formula
variance_partition_formula <- ~(1|Sample_taken_at)+(1|id_run_position)+(1|ANON_ID)+(1|an_alcohol_status)+an_est_age+an_scan1_bmi+(1|an_level_edu_str)+(1|an_diab_tf)+pn_monthdel+(1|an_eth_cat_str)+(1|an_smokstat)
# make sample info df with clinical covariates 
sample_info_vp <- sample.info %>%
  mutate(an_diab_tf = as.logical(an_diab)) %>% 
  arrange(factor(RNA_sanger_sample_id, levels = colnames(var_part_eig))) %>%
  column_to_rownames("RNA_sanger_sample_id")

# run variance partition (takes ~one minute)
varPart <- fitExtractVarPartModel(var_part_eig, variance_partition_formula, sample_info_vp)
saveRDS(varPart, file=paste0(outfiles_path, "module_variance_partition.RDS"))

varPart <- readRDS(paste0(outfiles_path, "module_variance_partition.RDS"))

# Plot vp (Supplementary Figure 27)
plotVarPart(sortCols(varPart), , col = c(rep("white", 11), "grey85")) + 
  scale_x_discrete(labels=c("ANON_ID" = "Individual", 
                            "id_run_position" = "Sequencing\nbatch", 
                            "Sample_taken_at" = "Sample\ntime-point", 
                            "an_eth_cat_str" = "Self-reported\nethnicity", 
                            "an_scan1_bmi" = "BMI", 
                            "pn_monthdel" = "Month of\ndelivery", 
                            "an_est_age" = "Age", 
                            "an_level_edu_str" = "Education\nlevel", 
                            "an_alcohol_status" = "Alcohol\nstatus", 
                            "an_smokstat" = "Smoking\nstatus", 
                            "an_diab_tf" = "Diabetes\nstatus")) + 
  xlab("Covariate") + 
  ggtitle("Variance partition plot for eigengenes") + 
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), 
        text=element_text(size=15), 
        plot.title = element_text(hjust=0))

# Check if the eigengene values are normally distributed 
eigengenes %>%
  rownames_to_column(var = "RNA_sanger_sample_id") %>%
  pivot_longer(cols = -RNA_sanger_sample_id, names_to = "Module", values_to = "Eigengene_value") %>%
  ggplot(aes(x = Eigengene_value)) + geom_histogram() + facet_wrap(~Module)
# Good enough 

# Fit the following model based on the variancePartition
# mixed.lmer <- lmer(Eigengene_value ~ SampleGA + (1|ANON_ID) + (1|id_run_position), data = eigen.time.data)
# merge sample info to eigengene info
eigen.time.data <- sample.info %>%
  left_join(eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id"), by = "RNA_sanger_sample_id") %>%
  mutate(ANON_ID = as.factor(ANON_ID), 
         Sample_taken_at = as.factor(Sample_taken_at))

#somewhere to store results
lme.results <- list()
lme_prediction_df <- c()
lme.partial.residuals <- c()

for (eigengene in colnames(eigengenes)) {
  res <- lmerTest::lmer(as.formula(get(eigengene) ~ splines::bs(x = SampleGA, degree = 3) + (1|ANON_ID) + (1|id_run_position)), 
                        data = eigen.time.data, 
                        REML = FALSE)
  
  lme.results[[eigengene]] <- res
  
  # based on the model, make a prediction 
  # use ggemmeans to estimate confidence intervals/fits
  prediction <- as.data.frame(ggemmeans(lme.results[[eigengene]], c("SampleGA [all]"))) %>% mutate(Mod = eigengene)
  colnames(prediction) <- c("SampleGA", "pred_eig_val", "std.err", "conf.low", "conf.high", "group", "Module") 
  lme_prediction_df <- rbind(lme_prediction_df, prediction)
  
  # extract partial residuals 
  # Using the argument partial.residuals = TRUE, what is plotted instead is the observed data with the effects of all the control variables accounted for
  plt <- interact_plot(model = lme.results[[eigengene]], pred = SampleGA, modx = id_run_position, plot.points = TRUE, partial.residuals = TRUE, data = eigen.time.data, colors = rep("black", length(unique(eigen.time.data$id_run_position))))
  partial_residuals <- plt$layers[[2]]$data %>% select(`get(eigengene)`)
  colnames(partial_residuals) <- "partial_residuals"
  SampleGA <- plt$layers[[2]]$data$SampleGA
  RNA_sanger_sample_id <- plt$layers[[2]]$data$RNA_sanger_sample_id
  partial_residuals_plt_dat <- as.data.frame(cbind(cbind(partial_residuals, SampleGA), RNA_sanger_sample_id)) %>% mutate(Mod = eigengene)
  lme.partial.residuals <- rbind(lme.partial.residuals, partial_residuals_plt_dat)
  
}

# Based on these models, make time series predictions for each module 
# This is an initial view of Main Figure 2, but the final figure will be created later 
lme_prediction_df %>% 
  mutate(Module = factor(as.numeric(gsub("ME_", "", Module)))) %>% 
  ggplot(aes(x = SampleGA)) + 
  geom_line(aes(SampleGA, pred_eig_val), size = 1.7, color = "red") +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high), alpha = 0.055, fill = "red") + 
  ylab("Module eigengene value") + xlab("Sample gestational age (weeks)") + 
  facet_wrap(~Module, ncol = 3) + 
  ggtitle("Signed hybrid consensus module eigengene trajectories") + 
  theme(legend.position = "none", strip.background = element_blank())

# cluster modules for displaying similar modules together systematically 
module_clustering = eigengenes %>% as.matrix() %>% t() %>% dist() %>% hclust(method ="average")
dendrogram <- as.dendrogram(module_clustering)
dend_rotated <- dendextend::rotate(dendrogram, order = c(1, 3, 7, 14, 15, 8, 4, 6, 5, 9, 12, 13, 11, 10,2))
module_order <- gsub("ME_", "", labels(dend_rotated))
# Plot the dendrogram (Supplementary Figure 11)
plot(dend_rotated, main = "Clustered module eigengenes")

################################################################################
# Create module figure (Main Figure 2)
################################################################################

# Prepare gene ID key 
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)

# Prep module annotations
# top 5 hub genes per mod (as a string)
top_5_hub_genes_per_module_cat <- connectivity %>%
  rownames_to_column(var = "gtf.gene_id") %>%
  left_join(id_name) %>%
  # find top kWithin for each module 
  group_by(Module) %>%
  top_n(n = 5, kWithin) %>%
  select(Module, gtf.gene_name) %>%
  group_by(Module) %>%
  mutate(top_5_hub_genes = paste(gtf.gene_name, collapse = ", ")) %>%
  select(-gtf.gene_name) %>% distinct() %>%
  mutate(Module = gsub("Module_", "", Module))

# top pathways 
# make a df
top_3_pathways_per_module_cat <- c()
for(i in 1:length(list_concise_eTerm_reactome)){
  pathways <- xEnrichViewer(list_concise_eTerm_reactome[[i]], 
                            top_num = length(list_concise_eTerm_reactome[[i]]$annotation)) %>% 
    filter(adjp < 0.05) %>% 
    arrange(desc(fc)) %>% 
    select(name, adjp, fc) %>% 
    top_n(n = 3, wt = fc) %>%
    mutate(Module = names(list_concise_eTerm_reactome)[[i]]) 
  top_3_pathways_per_module_cat <- rbind(top_3_pathways_per_module_cat, pathways)
}
top_3_pathways_per_module_cat_2 <- top_3_pathways_per_module_cat %>%
  select(Module, name) %>%
  # remove trailing spaces
  mutate(name = stringr::str_trim(name), 
         name = ifelse(name == "Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways", "NLR signaling pathways", name),
         name = ifelse(name == "Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.", "Mitochondrial respiration and heat production", name),
         name = ifelse(name == "RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function", "RUNX1 regulation of megakaryocyte and platelet genes", name),
         name = ifelse(name == "Hh mutants that don't undergo autocatalytic processing are degraded by ERAD", "ERAD degrades unprocessed Hh mutants", name),
         name = ifelse(name == "Antigen activates B Cell Receptor (BCR) leading to generation of second messengers", "Antigen activates BCR triggering second messengers", name),
         name = ifelse(name == "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell", "Lymphoid-non-Lymphoid immunoregulatory interactions", name),
         name = ifelse(name == "Formation of the ternary complex, and subsequently, the 43S complex", "Formation of ternary and 43S complexes", name)) %>%
  group_by(Module) %>%
  mutate(top_3_pathways = paste(name, collapse = "\n")) %>% 
  select(-name) %>% distinct() %>% 
  mutate(Module = gsub("Module_", "", Module)) 

mod_functions <- c("Innate immune response (m1)", 
                   "Cell proliferation (m2)", 
                   "Cell motility (m3)", 
                   "T cell immune response (m4)", 
                   "Protein synthesis/translation (m5)",
                   "B cell immune response (m6)",
                   "Eosinophil allergic-like response (m7)",
                   "Cytotoxic immune response (m8)",
                   "Protein sythesis/translation (m9)",
                   "Coagulation/clotting (m10)",
                   "Histones (m11)",
                   "Mitochondria/energy production (m12)",
                   "Small molecule transport (m13)",
                   "Antiviral interferon response (m14)",
                   "Inflammatory response regulation (m15)", 
                   "Unassigned")

mod_functions_df <- as.data.frame(cbind(Module = c(1:15, "Unassigned"), mod_functions))

# join these annotations into one string
annotations <- top_5_hub_genes_per_module_cat %>% 
  left_join(top_3_pathways_per_module_cat_2 , by = "Module") %>%
  mutate(top_3_pathways = ifelse(is.na(top_3_pathways), "None significant", top_3_pathways),
         joined_annot = paste0(top_5_hub_genes, "\n", top_3_pathways)) %>%
  left_join(mod_functions_df, by = "Module")

custom_labels <- setNames(annotations$mod_functions, annotations$Module)


module_name_labels <- as.data.frame(cbind(Module = annotations$Module[1:15], 
                                          xpos = rep(11.6, 15), 
                                          ypos = rep(0.015, 15),#rep(-0.02, 15), 
                                          label = annotations$joined_annot[1:15])) %>%
  mutate(xpos = as.numeric(xpos), ypos = as.numeric(ypos),
         Module = factor(Module, levels = module_order))

lme_prediction_df %>% 
  mutate(Module = gsub("ME_", "", Module), 
        Module = factor(Module, levels = module_order)) %>%
  ggplot(aes(x = SampleGA)) + 
  geom_line(aes(SampleGA, pred_eig_val), size = 1.7, color = "#32769B") +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high), alpha = 0.055, fill = "#32769B") + 
  ylab("Module eigengene") + xlab("Sample gestational age (weeks)") + 
  facet_wrap(~Module, ncol = 3, labeller = labeller(Module = custom_labels)) + 
  #ggtitle("Module trajectories") + 
  theme(legend.position = "none",
        strip.background = element_blank(),#element_rect(fill = "white", color = "black"),
        axis.text.y=element_blank(),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 14)) + 
  scale_x_continuous(expand = c(0.02,0),
                     breaks = c(12, 20, 28, 36)) +
  geom_text(
    data = module_name_labels,
    aes(x = xpos, y = ypos, label = label),
    hjust = 0, vjust = 0.95,
    #fontface = "bold",
    inherit.aes = FALSE,
    size = 3
  )

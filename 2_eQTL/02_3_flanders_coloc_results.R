# 02_3_flanders_coloc_results.R 

################################################################################

# 2.3. Evaluate colocalization results from flanders

################################################################################

# Aim: Identidfy colocalizations between eQTL and pregnancy GWAS

########################### Output paths ##########################
# specify paths 
path_lustre <- "genotyping/analysis/eQTL/"
path <- paste0(path_lustre, "output_data/flanders_5/flanders_output/results/")
interval_path <- "INTERVAL_RNAseq/"

########################### Input paths ###########################
# gtf file for gene id to gene name conversion
gtf_inpath <- "rna-seq/data/ref_data/Homo_sapiens.GRCh38.99.gtf"

########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)
library(locuszoomr)
library(AnnotationHub) 
library(GenomicRanges)

theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Load data ###########################
# read in gtf file
gtf <- rtracklayer::import(gtf_inpath)
# edit since we added two genes that are not in the reference 
HLA_update <- as.data.frame(rbind(c("HLA-DRB3", "HLA-DRB3"), c("HLA-DRB4", "HLA-DRB4")))
colnames(HLA_update) <- c("gtf.gene_id", "gtf.gene_name")
id_name <- rbind(data.frame(gtf$gene_id, gtf$gene_name) %>% distinct(), HLA_update)
genes_gtf <- gtf[gtf$type == "gene"]

# coloc results
coloc.results <- read.delim(paste0(path, "coloc/coloc_run_colocalization.table.all.tsv")) %>% arrange(desc(PP.H4.abf))
# mashR results
mashr.results_Interval <- read.table(paste0(path_lustre, "output_data/interval_comp/pops_alltps_intervalFRA_mashr_results_sharing_filtered_leads.txt"), sep="\t") %>%
  mutate(shared = factor(shared, levels = c("Bigger effect in POPS2", "Only significant in POPS2", "Bigger effect in Interval", "Opposite direction of effect", "shared", "not significant")))
# Interval finemapping results 
finemapped.loci_interval <- readRDS(paste0(path_lustre, "output_data/flanders_finemap/finemapped_loci_interval_full.rds"))
# Interval main effects results
Interval_pval_thresholds <- read.delim(paste0(interval_path, "results_summary_stats/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_phenotype_summary.tsv"), header = TRUE)
# Interval annotations
interval_annotations <- read.delim(paste0(interval_path, "Feature_Annotation_Ensembl_gene_ids_ensembl99.txt"))
# GWAS finemapped loci (for plotting)
finemapped.loci_GWAS <- readRDS(paste0(path_lustre, "output_data/flanders_finemap_2/finemapped_loci_gwas.rds"))
# POPs2 finemapped loci (for plotting loci)
finemapped.loci_pops <- readRDS(paste0(path_lustre, "output_data/flanders_finemap_2/finemapped_loci_pops.rds"))

########################### Analysis ###########################

# First, check the PE GWAS - do signals across the two studies coloc with each other? 
pe_gwas_colocs_tested <- coloc.results %>% 
  # filter to PE only 
  filter(grepl("PE", t1_study_id) & grepl("PE", t2_study_id),
         # filter to different studies only 
         t1_study_id != t2_study_id)
print("n signals tested for coloc between the two PE studies:")
pe_gwas_colocs_tested %>% nrow()
# [1] 18

print("n signals with significant coloc (PPH4 > 0.8) between the two PE studies:")
PE_gwas_overlap <- pe_gwas_colocs_tested %>% 
  # filter to significant colocs
  filter(PP.H4.abf > 0.8)

PE_gwas_overlap %>% nrow()
# [1] 10

# of these signals that coloc between the two GWAS studies, do they coloc with eQTL? 
PE_overlap_colocs <- coloc.results %>% 
  dplyr::filter(
    # filter to cross-dataset colocs
    t1_study_id != t2_study_id,
    # filter to colocs between a GWAS and an eQTL study 
    !grepl("GWAS", t1_study_id) | !grepl("GWAS", t2_study_id),
    grepl("GWAS", t1_study_id) | grepl("GWAS", t2_study_id)
  ) %>%
  # set info on GWAS vs eQTL datasets 
  mutate(GWAS_dataset = ifelse(grepl("GWAS", t1_study_id), t1_study_id, t2_study_id), 
         GWAS_hit = ifelse(grepl("GWAS", t1_study_id), hit1, hit2), 
         eQTL_dataset = ifelse(!grepl("GWAS", t1_study_id), t1_study_id, t2_study_id),
         eQTL_hit = ifelse(!grepl("GWAS", t1_study_id), hit1, hit2)) %>%
  dplyr::filter(
    # filter to colocs with GWAS hits that coloc between the two PE GWAS
    GWAS_hit %in% PE_gwas_overlap$hit1 | GWAS_hit %in% PE_gwas_overlap$hit2, 
    # filter to significant colocs
    PP.H4.abf > 0.8
  ) %>%
  dplyr::select(-t1_study_id, -t2_study_id, -hit1, -hit2)

print("PE GWAS signals that coloc with each other and also eQTL signals:")
unique(PE_overlap_colocs$GWAS_hit)
#[1] "chr15::PE_FINNGEN_GWAS::full::chr15:90887387:A:G::L1"            "chr16::PE_FINNGEN_GWAS::full::chr16:28516005:A:G::L1"           
#[3] "chr15::PE_GWAS_multiancestry_meta::full::chr15:90885406:A:G::L1" "chr16::PE_GWAS_multiancestry_meta::full::chr16:28516005:A:G::L1"

# two overlapping signals
PE_overlap_colocs %>% group_by(GWAS_hit, eQTL_dataset) %>% summarize(n = n()) %>% arrange(eQTL_dataset, GWAS_hit)
# GWAS_hit                                                        eQTL_dataset                  n
# <chr>                                                           <chr>                     <int>
# 1 chr16::PE_FINNGEN_GWAS::full::chr16:28516005:A:G::L1            Interval_females_u45_eQTL     2
# 2 chr16::PE_GWAS_multiancestry_meta::full::chr16:28516005:A:G::L1 Interval_females_u45_eQTL     3
# 3 chr15::PE_FINNGEN_GWAS::full::chr15:90887387:A:G::L1            Interval_full_eQTL            1
# 4 chr15::PE_GWAS_multiancestry_meta::full::chr15:90885406:A:G::L1 Interval_full_eQTL            2
# 5 chr16::PE_GWAS_multiancestry_meta::full::chr16:28516005:A:G::L1 Interval_full_eQTL            1
# 6 chr16::PE_FINNGEN_GWAS::full::chr16:28516005:A:G::L1            POPs2_bulk_all                1
# 7 chr16::PE_GWAS_multiancestry_meta::full::chr16:28516005:A:G::L1 POPs2_bulk_all                2

# based on these results, we will remove the following GWAS signals from further analysis
# to avoid double-counting the same signal twice
# chr16::PE_FINNGEN_GWAS::full::chr16:28516005:A:G::L1
# because chr16::PE_GWAS_multiancestry_meta::full::chr16:28516005:A:G::L1 colocs with 2 POPs2 signals (vs 1) and 3 Interval FRA signals (vs 2) and 1 INTERVAL signal (vs 0)

# chr15::PE_FINNGEN_GWAS::full::chr15:90887387:A:G::L1 
# because chr15::PE_GWAS_multiancestry_meta::full::chr15:90885406:A:G::L1 colocs with 2 INTERVAL signals (vs 1)

# to summarize: 
# we removed one of each GWAS signals that overlapped between the two PE GWAS datasets (PPH4 > 0.8)
# for each, we removed the GWAS signal that colocalized with fewer eQTL signals 

# now, remove these colocs from the coloc results and format 
coloc.results_crossdataset <-  coloc.results %>% 
  dplyr::filter(
    # retain colocs across datasets
    t1_study_id != t2_study_id,
    # retain colocs across GWAS and eQTL studies 
    !grepl("GWAS", t1_study_id) | !grepl("GWAS", t2_study_id),
    grepl("GWAS", t1_study_id) | grepl("GWAS", t2_study_id)
  ) %>%
  # make explicit GWAS vs eQTL datasets + hits 
  mutate(GWAS_dataset = ifelse(grepl("GWAS", t1_study_id), t1_study_id, t2_study_id), 
         GWAS_hit = ifelse(grepl("GWAS", t1_study_id), hit1, hit2), 
         eQTL_dataset = ifelse(!grepl("GWAS", t1_study_id), t1_study_id, t2_study_id),
         eQTL_hit = ifelse(!grepl("GWAS", t1_study_id), hit1, hit2),
         eGene_id = ifelse(!grepl("GWAS", t1_study_id), t1, t2)) %>%
  # remove the GWAS signals that overlap, decided above 
  filter(GWAS_hit != "chr16::PE_FINNGEN_GWAS::full::chr16:28516005:A:G::L1", 
         GWAS_hit != "chr15::PE_FINNGEN_GWAS::full::chr15:90887387:A:G::L1") %>%
  # remove old columns where GWAS and eQTL were not specified
  select(-t1, -t2, -hit1, -hit2, -t1_study_id, -t2_study_id) %>%
  # add eGene name 
  left_join(id_name %>% rename("eGene_name" = "gtf.gene_name"), by = c("eGene_id" = "gtf.gene_id"))

# write this out 
saveRDS(coloc.results_crossdataset, paste0(path_lustre, "output_data/flanders_5_coloc_outputs/coloc_results_crossdataset_all_noHG.rds"))
# save all coloc results for paper (only POPs2 and INTERVAL)
write.table(coloc.results_crossdataset %>%
              filter(eQTL_dataset == "Interval_full_eQTL" | eQTL_dataset == "POPs2_bulk_all"), 
            paste0(path_lustre, "output_data/flanders_5_coloc_outputs/coloc_results_crossdataset_paper.txt"), 
            quote = FALSE, row.names = FALSE)

# subset to significant colocalizations across datasets
sig.coloc.results_crossdataset <- coloc.results_crossdataset %>% 
  filter(PP.H4.abf > 0.8) 

saveRDS(sig.coloc.results_crossdataset, 
        paste0(path_lustre, "output_data/flanders_5_coloc_outputs/coloc_results_crossdataset_sig_flanders_5.rds"))

# How many GWAS were tested for coloc with POPs2? 
coloc.results_crossdataset %>% dplyr::filter(
  eQTL_dataset == "POPs2_bulk_all"
) %>% 
  nrow()
# [1] 1339

# how many were tested per GWAS study? 
# per gwas hit
coloc.results_crossdataset %>% dplyr::filter(
  eQTL_dataset == "POPs2_bulk_all"
) %>% 
  group_by(GWAS_dataset) %>%
  summarize(n_colocs_tested = n(), 
            n_unique_gwas_hits = length(unique(GWAS_hit)))
# GWAS_dataset                  n_colocs_tested n_unique_gwas_hits
# <chr>                                   <int>              <int>
#   1 GDM_GWAS_EUR_meta                    544                129
# 2 GestDura_GWAS_EUR_META                  97                 35
# 3 PE_FINNGEN_GWAS                        310                 68
# 4 PE_GWAS_multiancestry_meta             304                 59
# 5 sPTB_GWAS_EUR_META                      83                 22

# how many significant colocs? 
coloc.results_crossdataset %>% dplyr::filter(
  eQTL_dataset == "POPs2_bulk_all",
  PP.H4.abf > 0.8
) %>% 
  nrow()
# [1] 78

# significant colocalizations (PPH4 > 0.8) by GWAS
coloc.results_crossdataset %>% dplyr::filter(
  eQTL_dataset == "POPs2_bulk_all",
  PP.H4.abf > 0.8
) %>% 
  group_by(GWAS_dataset) %>%
  summarize(n_colocs_tested = n(), 
            n_unique_gwas_hits = length(unique(GWAS_hit)))
# GWAS_dataset               n_colocs_tested n_unique_gwas_hits
# <chr>                                <int>              <int>
# 1 GDM_GWAS_EUR_meta                       43                 28
# 2 GestDura_GWAS_EUR_META                   9                  6
# 3 PE_FINNGEN_GWAS                         12                  9
# 4 PE_GWAS_multiancestry_meta              12                  9
# 5 sPTB_GWAS_EUR_META                       2                  2

################################################################################
# make sup table for paper
gwas_name_key <- as.data.frame(rbind(c("GDM_GWAS_EUR_meta", "Gestational diabetes"),
                                     c("PE_FINNGEN_GWAS", "Pre-eclampsia"),
                                     c("PE_GWAS_multiancestry_meta", "Pre-eclampsia"),
                                     c("sPTB_GWAS_EUR_META", "Preterm birth"),
                                     c("GestDura_GWAS_EUR_META", "Gestational duration"))) %>%
  dplyr::rename("GWAS_name" = "V1", "GWAS_trait" = "V2")

# first, list all of the signals 
# list all of the signals
base_coloc_sup_table <- coloc.results_crossdataset %>% 
  left_join(gwas_name_key, by = c("GWAS_dataset" = "GWAS_name")) %>%
  # pick top signal per gene per gwas
  group_by(GWAS_hit, eGene_id, eQTL_dataset) %>%
  slice_max(n = 1, PP.H4.abf, with_ties = FALSE) %>% 
  pivot_wider(
    id_cols = c(GWAS_hit, eGene_name, eGene_id, GWAS_dataset, GWAS_trait),
    names_from = eQTL_dataset,
    values_from = c(
      PP.H4.abf,
      eQTL_hit
    ),
    names_sep = "_"
  ) %>%
  # filter to only rows where there is a sig hit in pops
  dplyr::filter(!is.na(PP.H4.abf_POPs2_bulk_all) & PP.H4.abf_POPs2_bulk_all > 0.8) %>% 
  arrange(desc(PP.H4.abf_POPs2_bulk_all)) %>%
  mutate(
    
    sig_POPs2    = !is.na(PP.H4.abf_POPs2_bulk_all) & PP.H4.abf_POPs2_bulk_all > 0.8,
    sig_Interval = !is.na(PP.H4.abf_Interval_full_eQTL) & PP.H4.abf_Interval_full_eQTL > 0.8,
    sig_Interval_FRA = !is.na(PP.H4.abf_Interval_females_u45_eQTL) & PP.H4.abf_Interval_females_u45_eQTL > 0.8,
    sig_Placenta = !is.na(PP.H4.abf_Placenta_eQTL) & PP.H4.abf_Placenta_eQTL > 0.8,
    
    significance_thesis = case_when(
      sig_POPs2 & sig_Interval & sig_Placenta ~ "Colocalised in all 3 datasets",
      sig_POPs2 & sig_Interval ~ "POPs2 and Interval",
      sig_POPs2 & sig_Placenta ~ "POPs2 and Placenta",
      sig_Interval & sig_Placenta ~ "Interval and Placenta",
      sig_POPs2 ~ "POPs2 only",
      sig_Interval ~ "Interval only",
      sig_Placenta ~ "Placenta only"
    ),
    
    significance_paper = case_when(
      sig_POPs2 & sig_Interval  ~ "Colocalised in both datasets",
      sig_POPs2 ~ "POPs2 only"
    ), 
    
    significance_all = case_when(
      sig_POPs2 & sig_Interval & sig_Placenta & sig_Interval_FRA ~ "Colocalised in all 4 datasets",
      
      sig_POPs2 & sig_Interval & sig_Placenta ~ "POPs2 and Interval and Placenta",
      sig_POPs2 & sig_Interval & sig_Interval_FRA ~ "POPs2 and Interval and IntervalFRA",
      sig_POPs2 & sig_Placenta & sig_Interval_FRA ~ "POPs2 and Placenta and IntervalFRA",
      sig_Interval & sig_Placenta & sig_Interval_FRA ~ "Interval and Placenta and IntervalFRA",
      
      sig_POPs2 & sig_Interval ~ "POPs2 and Interval",
      sig_POPs2 & sig_Placenta ~ "POPs2 and Placenta",
      sig_POPs2 & sig_Interval_FRA ~ "POPs2 and IntervalFRA",
      sig_Interval & sig_Placenta ~ "Interval and Placenta",
      sig_Interval & sig_Interval_FRA ~ "Interval and IntervalFRA",
      sig_Placenta & sig_Interval_FRA ~ "Placenta and IntervalFRA",
      
      sig_POPs2 ~ "POPs2 only",
      sig_Interval ~ "Interval only",
      sig_Placenta ~ "Placenta only",
      sig_Interval_FRA ~ "IntervalFRA only"
    ),
    
  )  %>% 
  ungroup() 

# join this to info about genes 
interval_eQTL_info <- base_coloc_sup_table %>%
  dplyr::filter(sig_Interval == FALSE) %>%
  left_join(interval_annotations %>%
              dplyr::select(feature_id, gene_biotype) %>%
              distinct(), 
            by = c("eGene_id" = "feature_id")) %>%
  mutate(tested_for_coloc = ifelse(is.na(PP.H4.abf_Interval_full_eQTL), FALSE, TRUE),
         tested_for_eQTL_not_sig = eGene_id %in% (Interval_pval_thresholds %>% dplyr::filter(qval_sig == FALSE) %>% pull(phenotype_id)),
         sig_eQTL_not_finemapped = eGene_id %in% (Interval_pval_thresholds %>% dplyr::filter(qval_sig == TRUE) %>% pull(phenotype_id)),
         eGene_finemapped = eGene_id %in% unique(finemapped.loci_interval$phenotype_id),
         reason_not_coloc = case_when(
           tested_for_coloc ~ "tested for coloc, pph4 < 0.8",
           eGene_finemapped ~ "finemapped, not tested for coloc",
           sig_eQTL_not_finemapped ~ "sig eQTL, not finemapped",
           tested_for_eQTL_not_sig ~ "tested for eQTL, not sig",
           grepl("pseudogene", gene_biotype) ~ "pseudogene filtered from interval",
           !grepl("pseudogene", gene_biotype) ~ "not expressed in interval"
         ), 
         reason_not_coloc = factor(reason_not_coloc, levels = c("tested for coloc, pph4 < 0.8", 
                                                                "finemapped, not tested for coloc", 
                                                                "sig eQTL, not finemapped", 
                                                                "tested for eQTL, not sig", 
                                                                "pseudogene filtered from interval",
                                                                "not expressed in interval"))
         
  )
table(interval_eQTL_info$reason_not_coloc, useNA = "always")

# join this info back to the coloc table
base_coloc_sup_table <- base_coloc_sup_table %>% 
  left_join(interval_eQTL_info %>%
              dplyr::select(GWAS_hit, eGene_name, eGene_id, 
                     GWAS_dataset, GWAS_trait, PP.H4.abf_POPs2_bulk_all, eQTL_hit_POPs2_bulk_all, 
                     reason_not_coloc)) %>%
  # add mashR info
  left_join(mashr.results_Interval %>% separate(Pairs_notrsID, into = c("chr", "pos", "eGene")) %>% 
              dplyr::select(-Interval, -POPS2_effect, -pos, -chr) %>% 
              pivot_wider(names_from = POPS2_timepoint, values_from = shared), 
            by = c("eGene_id" = "eGene")) %>%
  # re-arrange
  relocate(POPS2_28wk, .after = POPS2_20wk)

base_coloc_sup_table %>% write.table(paste0(path_lustre, "output_data/flanders_5_coloc_outputs/coloc_summary_table.csv"), 
                                     quote = FALSE, row.names = FALSE) 

# thesis sup 
paper_sup <- base_coloc_sup_table %>% 
  separate(GWAS_hit, into = c("Chr", "GWAS_dataset", "Trait", "GWAS_hit_SNP", "Locus_signal"), sep = "::") %>%
  separate(eQTL_hit_Interval_full_eQTL, into = c("Chr", "eQTL_dataset_Interval", "eGene_id_interval", "Interval_eSNP", "Locus_signal"), sep = "::") %>%
  separate(eQTL_hit_POPs2_bulk_all, into = c("Chr", "eQTL_dataset_POPS", "eGene_id_POPS", "POPs2_eSNP", "Locus_signal"), sep = "::") %>%
  separate(GWAS_hit_SNP, into = c("Chr", "Pos", "minor", "major"), sep = ":", remove = FALSE) %>%
  mutate(Pregnancy_enhanced = paste0(POPS2_12wk, POPS2_20wk, POPS2_28wk, POPS2_36wk), 
         Chr = as.numeric(gsub("chr", "", Chr)), 
         PP.H4.abf_Interval_full_eQTL = round(PP.H4.abf_Interval_full_eQTL, digits = 2), 
         PP.H4.abf_POPs2_bulk_all = round(PP.H4.abf_POPs2_bulk_all, digits = 2), 
         GWAS_study = case_when(
           GWAS_dataset == "PE_FINNGEN_GWAS" ~ "Pre-eclampsia (FinnGen)",
           GWAS_dataset == "PE_GWAS_multiancestry_meta" ~ "Pre-eclampsia (multi-ancestry)",
           GWAS_dataset == "sPTB_GWAS_EUR_META" ~ "Preterm birth (European)", 
           GWAS_dataset == "GestDura_GWAS_EUR_META" ~ "Gestational duration (European)", 
           GWAS_dataset == "GDM_GWAS_EUR_meta" ~ "Gestational diabetes mellitus (European)", 
         )) %>%
  arrange(GWAS_study, desc(significance_paper), Chr, Pos, reason_not_coloc) %>%
  dplyr::select(GWAS_study, GWAS_hit_SNP, eGene_name, eGene_id, PP.H4.abf_POPs2_bulk_all, PP.H4.abf_Interval_full_eQTL, POPs2_eSNP, Interval_eSNP, significance_paper, reason_not_coloc) %>%
  rename("reason_not_coloc" = "Interval_coloc_info", "significance_paper" = "significant_coloc") 

# save this 
paper_sup %>% write.csv(paste0(path_lustre, "output_data/flanders_5_coloc_outputs/coloc_summary_table_paper_sup.csv"), 
                        row.names = FALSE, na = "") 

########################### Plot ###########################

# Plot n finemapped GWAS loci vs. colocs (Supplementary Figure 14)
gwas_loci <- finemapped.loci_GWAS %>%
  group_by(study_id) %>%
  summarize(n_loci = n())
coloc.results_crossdataset %>%
  dplyr::filter(eQTL_dataset == "POPs2_bulk_all", 
         PP.H4.abf > 0.8) %>%
  group_by(GWAS_dataset) %>%
  summarize(n_colocs = n()) %>%
  left_join(gwas_loci, by = c("GWAS_dataset" = "study_id")) %>%
  mutate(dataset = case_when(
    GWAS_dataset == "GDM_GWAS_EUR_meta" ~ "GDM GWAS",
    GWAS_dataset == "PE_FINNGEN_GWAS" ~ "PE GWAS (FinnGen)",
    GWAS_dataset == "PE_GWAS_multiancestry_meta" ~ "PE GWAS (multi-ancestry)",
    GWAS_dataset == "sPTB_GWAS_EUR_META" ~ "PTB GWAS",
    GWAS_dataset == "GestDura_GWAS_EUR_META" ~ "Gestational duration GWAS")) %>%
  ggplot(aes(x = n_loci, y = n_colocs)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = dataset)) +
  ggpubr::stat_cor() + 
  labs(title = "Correlation between n GWAS loci and n colocalizations with pregnancy eQTL", 
       x = "N GWAS loci (p < 1e-5)", 
       y = "N significant colocalizations (PPH4 > 0.8)")

# Barchart summarizing n colocs with POPs and Interval
colocs.summary_df <- 
  coloc.results_crossdataset %>%
  mutate(GWAS_trait = case_when(
    GWAS_dataset == "GDM_GWAS_EUR_meta" ~ "Gestational diabetes",
    GWAS_dataset == "PE_FINNGEN_GWAS" ~ "Pre-eclampsia",
    GWAS_dataset == "PE_GWAS_multiancestry_meta" ~ "Pre-eclampsia",
    GWAS_dataset == "sPTB_GWAS_EUR_META" ~ "Preterm birth",
    GWAS_dataset == "GestDura_GWAS_EUR_META" ~ "Gestational duration"
  ), 
  GWAS_trait = factor(GWAS_trait, levels = c("Gestational diabetes", "Pre-eclampsia", "Gestational duration", "Preterm birth"))) %>%
  dplyr::filter(eQTL_dataset != "Interval_females_u45_eQTL") %>%
  # pick top signal per gene per gwas
  group_by(GWAS_hit, eGene_id, eQTL_dataset) %>%
  slice_max(n = 1, PP.H4.abf, with_ties = FALSE) %>% 
  pivot_wider(
    id_cols = c(GWAS_hit, eGene_name, eGene_id, GWAS_dataset, GWAS_trait),
    names_from = eQTL_dataset,
    values_from = c(
      PP.H4.abf,
      eQTL_hit
    ),
    names_sep = "_"
  ) %>%
  mutate(
    
    sig_POPs2    = !is.na(PP.H4.abf_POPs2_bulk_all) & PP.H4.abf_POPs2_bulk_all > 0.8,
    sig_Interval = !is.na(PP.H4.abf_Interval_full_eQTL) & PP.H4.abf_Interval_full_eQTL > 0.8,
    
    significance_paper = case_when(
      sig_POPs2 & sig_Interval  ~ "POPs2 and INTERVAL",
      sig_POPs2 ~ "POPs2 only",
      sig_Interval ~ "INTERVAL only"
    ), 
    gwas_trait_forplt = case_when(
      GWAS_trait == "Pre-eclampsia" ~ "Pre-\neclampsia",
      GWAS_trait == "Gestational diabetes" ~ "Gestational\ndiabetes",
      GWAS_trait == "Gestational duration" ~ "Gestational\nduration",
      GWAS_trait == "Preterm birth" ~ "Preterm\nbirth"
    ), 
    gwas_trait_forplt = factor(gwas_trait_forplt, levels = c("Gestational\ndiabetes", "Pre-\neclampsia", "Gestational\nduration", "Preterm\nbirth"))
  ) 

# Barplot including Interval 
colocs.summary_df %>% 
  # filter to just rows where there is a sig eQTL
  dplyr::filter(sig_POPs2 | sig_Interval) %>% 
  group_by(significance_paper, gwas_trait_forplt) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(significance_paper = factor(significance_paper, levels = c("POPs2 only", "POPs2 and INTERVAL", "INTERVAL only"))) %>%
  ggplot(aes(x = gwas_trait_forplt, y = n, fill = significance_paper)) +
  geom_col() +
  scale_fill_manual(values = c(
    "POPs2 and INTERVAL" = "#f3a257",
    "POPs2 only" = "#598eca",
    "INTERVAL only" = "#cb2f43"
    
  )) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "white") +
  labs(
    x = "GWAS trait",
    y = "Number of GWAS locus–eGene pairs",
    fill = "Colocalization category"
  ) +
  #ggtitle("Summary of colocalizations in POPs2 and INTERVAL by GWAS trait") + 
  #coord_flip() + 
  theme(legend.position = "top", 
        legend.title.position = "top", 
        legend.title = element_text(hjust = 0.5))

# Barplot excluding Interval (Main Figure 4A)
colocs.summary_df %>% 
  # filter to just rows where there is a sig eQTL
  dplyr::filter(sig_POPs2) %>% 
  group_by(significance_paper, gwas_trait_forplt) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(significance_paper = factor(significance_paper, levels = c("POPs2 only", "POPs2 and INTERVAL", "INTERVAL only"))) %>%
  ggplot(aes(x = gwas_trait_forplt, y = n, fill = significance_paper)) +
  geom_col() +
  scale_fill_manual(values = c(
    "POPs2 and INTERVAL" = "#9baa77",
    "POPs2 only" = "#598eca"#,
    #"INTERVAL only" = "#cb2f43"
    
  )) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "white") +
  labs(
    x = "GWAS trait",
    y = "Number of GWAS locus–eGene pairs",
    fill = "Colocalization category"
  ) +
  #ggtitle("Summary of colocalizations in POPs2 and INTERVAL by GWAS trait") + 
  #coord_flip() + 
  theme(legend.position = "top", 
        legend.title.position = "top", 
        legend.title = element_text(hjust = 0.5))

# Locus plots 
# Prepare annotations for locuszoom
hub <- AnnotationHub()
query(hub, c("EnsDb", "sapiens", "GRCh38", "99"))
edb <- hub[["AH78783"]]


# find nearest gene to the variants with coloc for abstract
nearest_gene_lookup <- function(hit1, genes_gtf) {
  parts1 <- strsplit(hit1, "::", fixed = TRUE)[[1]]
  parts <- strsplit(parts1[4], ":", fixed = TRUE)[[1]]
  chr <- gsub("^chr", "", parts[1])   # remove "chr" prefix
  pos <- as.integer(parts[2])
  
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(start = pos, end = pos))
  
  # keep only seqlevels that exist in genes_gtf
  common <- intersect(seqlevels(gr), seqlevels(genes_gtf))
  gr <- keepSeqlevels(gr, common, pruning.mode = "coarse")
  genes_sub <- keepSeqlevels(genes_gtf, common, pruning.mode = "coarse")
  
  nearest_idx <- nearest(gr, genes_sub)
  nearest_gene <- genes_sub[nearest_idx]
  
  return(mcols(nearest_gene)$gene_name)
}

# function 
plot_interval_pops_coloc <- function(coloc.results_crossdataset_pops, i, finemapped.loci_pops, finemapped.loci_GWAS, finemapped.loci_interval, genes_gtf, interval_eQTL_dataset = "Interval_full_eQTL"){
  gene_name = coloc.results_crossdataset_pops$eGene_name[i]
  gwas_hit = coloc.results_crossdataset_pops$GWAS_hit[i]
  gwas_dataset <- strsplit(gwas_hit, "::")[[1]][2]
  chr <- gsub("chr", "", strsplit(gwas_hit, "::")[[1]][1])
  lead_pos <- strsplit(strsplit(gwas_hit, "::")[[1]][4], ":")[[1]][2]
  pops_eQTL_hit <- coloc.results_crossdataset_pops$eQTL_hit[i]
  # this removes second hits...
  #pops_eQTL_hit <- coloc.results_crossdataset %>% dplyr::filter(eGene_name == gene_name, GWAS_hit == gwas_hit, eQTL_dataset == "POPs2_bulk_all", eQTL_hit == eQTL_hit_i) %>% pull(eQTL_hit)
  #pops_eQTL_hit <- coloc.results_crossdataset %>% dplyr::filter(eGene_name == gene_name, GWAS_hit == gwas_hit, eQTL_dataset == "POPs2_bulk_all") %>% slice_max(PP.H4.abf, n = 1) %>% pull(eQTL_hit)
  pops_pph4 = coloc.results_crossdataset %>% dplyr::filter(eGene_name == gene_name, GWAS_hit == gwas_hit, eQTL_dataset == "POPs2_bulk_all", eQTL_hit == pops_eQTL_hit) %>% pull(PP.H4.abf)
  ld_inpath <- paste0("/lustre/scratch125/humgen/projects_v2/pops2/genotyping/analysis/eQTL/output_data/coloc_LD/", gwas_dataset, "_chr", chr, "_", lead_pos, ".ld")
  
  
  pops_eQTL_susie_path <- finemapped.loci_pops %>% dplyr::filter(credible_set_name == pops_eQTL_hit) %>% pull(path_rds)
  GWAS_susie_path <- finemapped.loci_GWAS %>% dplyr::filter(credible_set_name == gwas_hit) %>% pull(path_rds)
  # read in susie files 
  pops_eQTL_susie <- readRDS(pops_eQTL_susie_path)
  # find the correct locus
  j1 = which(names(pops_eQTL_susie) == pops_eQTL_hit)
  pops_eQTL_susie <- pops_eQTL_susie[[j1]]
  GWAS_susie <- readRDS(GWAS_susie_path)
  # find the correct locus
  j3 = which(names(GWAS_susie) == gwas_hit)
  GWAS_susie <- GWAS_susie[[j3]]
  
  
  # check if there is an interval coloc tested at this gene
  interval_eQTL_hit <- coloc.results_crossdataset %>% dplyr::filter(eGene_name == gene_name, GWAS_hit == gwas_hit, eQTL_dataset == interval_eQTL_dataset) %>% slice_max(PP.H4.abf, n = 1) %>% pull(eQTL_hit)
  if(length(interval_eQTL_hit) == 0){
    common_snps = intersect(
      pops_eQTL_susie$finemapping_lABFs$position,
      GWAS_susie$finemapping_lABFs$position)
    
    # pick lead SNP from POPs and GWAS 
    common_max_snp <- full_join(pops_eQTL_susie$finemapping_lABFs %>% rename("lABF" = "lABF_eQTL", "is_cs" = "is_cs_eQTL") %>% dplyr::select(-bC, -bC_se), 
                                GWAS_susie$finemapping_lABFs %>% rename("lABF" = "lABF_GWAS", "is_cs" = "is_cs_GWAS")%>% dplyr::select(-bC, -bC_se)) %>%
      dplyr::filter(!is.na(lABF_GWAS),
                    !is.na(lABF_eQTL),
                    is_cs_GWAS,
                    is_cs_eQTL) %>%
      mutate(
        rank_eQTL = rank(-lABF_eQTL, ties.method = "min"),
        rank_GWAS = rank(-lABF_GWAS, ties.method = "min"),
        joint_rank = rank_eQTL + rank_GWAS
      ) %>%
      slice_min(n = 1, order_by = joint_rank, with_ties = FALSE) %>%
      pull(snp)
    
    coloc_lead_pos <- strsplit(common_max_snp, ":")[[1]][2]
    
    # # get the GWAS hit 
    # coloc_lead_pos <- sub(".*:(\\d+):.*", "\\1", coloc.results_crossdataset$GWAS_hit[i])
    # # check if the gwas lead pos is in common snps
    # coloc_lead_pos %in% common_snps
    
    # read in LD
    ld <- read.table(paste0(ld_inpath),header=T)
    ld1 <- ld[ld$BP_A ==coloc_lead_pos,]
    ld2 <- ld[ld$BP_B ==coloc_lead_pos,]
    ld2$BP_B <- ld2$BP_A
    ld <- rbind(ld1,ld2,
                # add a row that is that SNP with itself
                c(ld1[1, 1], ld1[1, 2], ld1[1,3], ld1[1,4], ld1[1, 1], ld1[1, 2], ld1[1,3], ld1[1,4], 1))
    
    coloc_plt_df <- rbind(pops_eQTL_susie$finemapping_lABFs %>% dplyr::filter(position %in% common_snps) %>% mutate(phenotype = gene_name, dataset = pops_eQTL_susie$metadata$study_id),
                          #interval_eQTL_susie$finemapping_lABFs %>% dplyr::filter(position %in% common_snps) %>% mutate(phenotype = gene_name, dataset = interval_eQTL_susie$metadata$study_id),
                          GWAS_susie$finemapping_lABFs %>% dplyr::filter(position %in% common_snps) %>% mutate(phenotype = nearest_gene_lookup(gwas_hit, genes_gtf), dataset = GWAS_susie$metadata$study_id)) %>%
      left_join(ld %>% mutate(BP_B = as.numeric(BP_B),
                              R2 = as.numeric(R2)), by = c("position" = "BP_B")) %>%
      mutate(
        dataset = case_when(
          dataset == "GDM_GWAS_EUR_meta" ~ "Gestational diabetes GWAS",
          dataset == "PE_FINNGEN_GWAS" ~ "Pre-eclampsia GWAS (FinnGen)",
          dataset == "PE_GWAS_multiancestry_meta" ~ "Pre-eclampsia GWAS (multi-ancestry)",
          dataset == "sPTB_GWAS_EUR_META" ~ "Preterm birth GWAS",
          dataset == "GestDura_GWAS_EUR_META" ~ "Gestational duration GWAS",
          dataset == "POPs2_bulk_all" ~ "POPs2 eQTL",
          dataset == "Interval_full_eQTL" ~ "INTERVAL eQTL"
        ),
        gene_dataset = ifelse(grepl("GWAS", dataset), 
                              paste0(dataset, " locus"),#paste0(dataset, ": ", phenotype, " locus"), 
                              paste0(dataset, ": ", phenotype, " expression (", "PPH4: ", round(pops_pph4, 3), ")")),
        chrom = chr, 
        pos = position, 
        rsid = snp,
        p = lABF,
        r2 = R2, 
        pch = ifelse(snp == common_max_snp, 24, 21),
        # R2_group = cut(
        #   R2,
        #   breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
        #   labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
        #   include.lowest = TRUE,
        #   right = FALSE
        # ), 
        #R2_group = factor(R2_group, levels = rev(levels(R2_group)))
        bg = as.character(cut(
          R2,
          breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
          labels = c("#00117F", "#4E9CFF", "#90BE6D", "#F79256", "#D33F49"),
          include.lowest = TRUE,
          right = FALSE
        ))) %>%
      dplyr::select(chrom, pos, p, snp, r2, lABF, bg, pch, dataset, gene_dataset)
    
    loc_eQTL <- locus(
      data = coloc_plt_df %>%
        dplyr::filter(grepl("POPs2", dataset)),
      seqname = chr,                     # chromosome
      xrange = c(min(coloc_plt_df$pos), 
                 max(coloc_plt_df$pos)),  # bp range
      ens_db = edb, #"EnsDb.Hsapiens.v86", 
      yvar = "lABF", LD = "r2",
      index_snp = common_max_snp
    )
    
    loc_GWAS <- locus(
      data = coloc_plt_df %>%
        dplyr::filter(grepl("GWAS", dataset)),
      seqname = chr,                     # chromosome
      xrange = c(min(coloc_plt_df$pos), max(coloc_plt_df$pos)),  # bp range
      ens_db = edb, #"EnsDb.Hsapiens.v86", 
      yvar = "lABF", LD = "r2",
      index_snp = common_max_snp
    )
    
    GWAS_name = unique(coloc_plt_df$gene_dataset)[grepl("GWAS", unique(coloc_plt_df$gene_dataset))]
    eQTL_name = unique(coloc_plt_df$gene_dataset)[!grepl("GWAS", unique(coloc_plt_df$gene_dataset))]
    
    oldpar <- set_layers(2)
    scatter_plot(loc_GWAS, xticks = FALSE, use_layout = FALSE, legend_pos = "topleft", col = NA, showExons = TRUE, border = TRUE)
    mtext(
      GWAS_name,
      side = 3,
      line = 2.6,
      outer = FALSE,
      cex = 1,
      font = 1
    )
    scatter_plot(loc_eQTL, xticks = FALSE, use_layout = FALSE, legend_pos = NULL, col = NA, showExons = TRUE, border = TRUE)
    mtext(
      eQTL_name,
      side = 3,     # top
      line = 2.6,     # move outside panel
      outer = FALSE,
      cex = 1,
      font = 1      # bold
    )
    genetracks(loc_GWAS, border = TRUE)
    par(oldpar)  # revert par() settings
    
  } else {
    interval_pph4 = coloc.results_crossdataset %>% dplyr::filter(eGene_name == gene_name, GWAS_hit == gwas_hit, eQTL_dataset == interval_eQTL_dataset, eQTL_hit == interval_eQTL_hit) %>% pull(PP.H4.abf)
    interval_eQTL_susie_path <- finemapped.loci_interval %>% dplyr::filter(credible_set_name == interval_eQTL_hit) %>% pull(path_rds)
    interval_eQTL_susie <- readRDS(interval_eQTL_susie_path)
    # find the correct locus
    j2 = which(names(interval_eQTL_susie) == interval_eQTL_hit)
    interval_eQTL_susie <- interval_eQTL_susie[[j2]]
    
    common_snps = intersect(intersect(
      pops_eQTL_susie$finemapping_lABFs$position,
      interval_eQTL_susie$finemapping_lABFs$position),
      GWAS_susie$finemapping_lABFs$position)
    
    # pick lead SNP from POPs and GWAS 
    common_max_snp <- full_join(pops_eQTL_susie$finemapping_lABFs %>% rename("lABF" = "lABF_eQTL", "is_cs" = "is_cs_eQTL") %>% dplyr::select(-bC, -bC_se), 
                                GWAS_susie$finemapping_lABFs %>% rename("lABF" = "lABF_GWAS", "is_cs" = "is_cs_GWAS")%>% dplyr::select(-bC, -bC_se)) %>%
      dplyr::filter(!is.na(lABF_GWAS),
                    !is.na(lABF_eQTL),
                    is_cs_GWAS,
                    is_cs_eQTL) %>%
      mutate(
        rank_eQTL = rank(-lABF_eQTL, ties.method = "min"),
        rank_GWAS = rank(-lABF_GWAS, ties.method = "min"),
        joint_rank = rank_eQTL + rank_GWAS
      ) %>%
      slice_min(n = 1, order_by = joint_rank, with_ties = FALSE) %>%
      pull(snp)
    
    coloc_lead_pos <- strsplit(common_max_snp, ":")[[1]][2]
    
    # read in LD
    ld <- read.table(paste0(ld_inpath),header=T)
    ld1 <- ld[ld$BP_A ==coloc_lead_pos,]
    ld2 <- ld[ld$BP_B ==coloc_lead_pos,]
    ld2$BP_B <- ld2$BP_A
    ld <- rbind(ld1,ld2,
                # add a row that is that SNP with itself
                c(ld1[1, 1], ld1[1, 2], ld1[1,3], ld1[1,4], ld1[1, 1], ld1[1, 2], ld1[1,3], ld1[1,4], 1))
    
    coloc_plt_df <- rbind(pops_eQTL_susie$finemapping_lABFs %>% dplyr::filter(position %in% common_snps) %>% mutate(phenotype = gene_name, dataset = pops_eQTL_susie$metadata$study_id),
                          interval_eQTL_susie$finemapping_lABFs %>% dplyr::filter(position %in% common_snps) %>% mutate(phenotype = gene_name, dataset = interval_eQTL_susie$metadata$study_id),
                          GWAS_susie$finemapping_lABFs %>% dplyr::filter(position %in% common_snps) %>% mutate(phenotype = nearest_gene_lookup(gwas_hit, genes_gtf), dataset = GWAS_susie$metadata$study_id)) %>%
      left_join(ld %>% mutate(BP_B = as.numeric(BP_B),
                              R2 = as.numeric(R2)), by = c("position" = "BP_B")) %>%
      mutate(
        dataset = case_when(
          dataset == "GDM_GWAS_EUR_meta" ~ "Gestational diabetes GWAS",
          dataset == "PE_FINNGEN_GWAS" ~ "Pre-eclampsia GWAS (FinnGen)",
          dataset == "PE_GWAS_multiancestry_meta" ~ "Pre-eclampsia GWAS (multi-ancestry)",
          dataset == "sPTB_GWAS_EUR_META" ~ "Preterm birth GWAS",
          dataset == "GestDura_GWAS_EUR_META" ~ "Gestational duration GWAS",
          dataset == "POPs2_bulk_all" ~ "POPs2 eQTL",
          dataset == "Interval_full_eQTL" ~ "INTERVAL eQTL"
        ),
        gene_dataset = ifelse(grepl("GWAS", dataset), 
                              paste0(dataset, " locus"), #paste0(dataset, ": ", phenotype, " locus"), 
                              ifelse(grepl("POPs2", dataset), 
                                     paste0(dataset, ": ", phenotype, " expression (", "PPH4: ", round(pops_pph4, 3), ")"), 
                                     paste0(dataset, ": ", phenotype, " expression (", "PPH4: ", round(interval_pph4, 3), ")"))),
        chrom = chr, 
        pos = position, 
        rsid = snp,
        p = lABF,
        r2 = R2, 
        pch = ifelse(snp == common_max_snp, 24, 21),
        # R2_group = cut(
        #   R2,
        #   breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
        #   labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
        #   include.lowest = TRUE,
        #   right = FALSE
        # ),
        # R2_group = factor(R2_group, levels = rev(levels(R2_group)))
        bg = as.character(cut(
          R2,
          breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
          labels = c("#00117F", "#4E9CFF", "#90BE6D", "#F79256", "#D33F49"),
          include.lowest = TRUE,
          right = FALSE
        ))) %>%
      dplyr::select(chrom, pos, p, snp, r2, lABF, bg, pch, dataset, gene_dataset)
    
    loc_POPs_eQTL <- locus(
      data = coloc_plt_df %>%
        dplyr::filter(grepl("POPs2", dataset)),
      seqname = chr,                     # chromosome
      xrange = c(min(coloc_plt_df$pos), 
                 max(coloc_plt_df$pos)),  # bp range
      ens_db = edb, #"EnsDb.Hsapiens.v86", 
      yvar = "lABF", LD = "r2",
      index_snp = common_max_snp
    )
    
    loc_Interval_eQTL <- locus(
      data = coloc_plt_df %>%
        dplyr::filter(grepl("INTERVAL", dataset)),
      seqname = chr,                     # chromosome
      xrange = c(min(coloc_plt_df$pos), 
                 max(coloc_plt_df$pos)),  # bp range
      ens_db = edb, #"EnsDb.Hsapiens.v86", 
      yvar = "lABF", LD = "r2",
      index_snp = common_max_snp
    )
    
    loc_GWAS <- locus(
      data = coloc_plt_df %>%
        dplyr::filter(grepl("GWAS", dataset)),
      seqname = chr,                     # chromosome
      xrange = c(min(coloc_plt_df$pos), max(coloc_plt_df$pos)),  # bp range
      ens_db = edb, #"EnsDb.Hsapiens.v86", 
      yvar = "lABF", LD = "r2",
      index_snp = common_max_snp
    )
    
    GWAS_name = unique(coloc_plt_df$gene_dataset)[grepl("GWAS", unique(coloc_plt_df$gene_dataset))]
    eQTL_POPS_name = unique(coloc_plt_df$gene_dataset)[grepl("POPs2", unique(coloc_plt_df$gene_dataset))]
    eQTL_Interval_name = unique(coloc_plt_df$gene_dataset)[grepl("INTERVAL", unique(coloc_plt_df$gene_dataset))]
    
    oldpar <- set_layers(3)
    scatter_plot(loc_GWAS, xticks = FALSE, use_layout = FALSE, legend_pos = "topleft", col = NA, showExons = TRUE, border = TRUE)
    mtext(
      GWAS_name,
      side = 3,
      line = 2.6,
      outer = FALSE,
      cex = 1,
      font = 1
    )
    scatter_plot(loc_POPs_eQTL, xticks = FALSE, use_layout = FALSE, legend_pos = NULL, col = NA, showExons = TRUE, border = TRUE)
    mtext(
      eQTL_POPS_name,
      side = 3,     # top
      line = 2.6,     # move outside panel
      outer = FALSE,
      cex = 1,
      font = 1      # bold
    )
    scatter_plot(loc_Interval_eQTL, xticks = FALSE, use_layout = FALSE, legend_pos = NULL, col = NA, showExons = TRUE, border = TRUE)
    mtext(
      eQTL_Interval_name,
      side = 3,     # top
      line = 2.6,     # move outside panel
      outer = FALSE,
      cex = 1,
      font = 1      # bold
    )
    genetracks(loc_GWAS, border = TRUE)
    par(oldpar)  # revert par() settings
    
  }
  
}

# subset the coloc df to colocs with POPS
coloc.results_crossdataset_pops <- coloc.results_crossdataset %>% 
  dplyr::filter(eQTL_dataset == "POPs2_bulk_all", 
                PP.H4.abf > 0.8)

# plot loci for paper
# NOTCH2 (Main Figure 4B)
gene_name = "NOTCH2"
plot_interval_pops_coloc(coloc.results_crossdataset_pops, 
                         i = which(coloc.results_crossdataset_pops$eGene_name == gene_name), 
                         finemapped.loci_pops, finemapped.loci_GWAS, finemapped.loci_interval, genes_gtf, interval_eQTL_dataset = "Interval_full_eQTL")

# OPRL1 (Supplementary Figure 15)
gene_name = "OPRL1"
plot_interval_pops_coloc(coloc.results_crossdataset_pops, 
                         i = which(coloc.results_crossdataset_pops$eGene_name == gene_name)[1], 
                         finemapped.loci_pops, finemapped.loci_GWAS, finemapped.loci_interval, genes_gtf, interval_eQTL_dataset = "Interval_full_eQTL")

# VEGFA (Supplementary Figure 16)
gene_name = "VEGFA"
plot_interval_pops_coloc(coloc.results_crossdataset_pops, 
                         i = which(coloc.results_crossdataset_pops$eGene_name == gene_name)[1], 
                         finemapped.loci_pops, finemapped.loci_GWAS, finemapped.loci_interval, genes_gtf, interval_eQTL_dataset = "Interval_full_eQTL")

# MSMO1
gene_name = "MSMO1"
plot_interval_pops_coloc(coloc.results_crossdataset_pops, 
                         i = which(coloc.results_crossdataset_pops$eGene_name == gene_name), 
                         finemapped.loci_pops, finemapped.loci_GWAS, finemapped.loci_interval, genes_gtf, interval_eQTL_dataset = "Interval_full_eQTL")

# To plot all loci to a pdf, run a loop with a pdf open like this
# open a pdf device
pdf(paste0(path, "locuszoom_all_sig_pops_w_interval_plots.pdf"), width = 9)

# Loop over coloc results
for (i in 1:nrow(coloc.results_crossdataset_pops)) {
  
  print(i)
  
  plot_interval_pops_coloc(coloc.results_crossdataset_pops, i, 
                           finemapped.loci_pops, finemapped.loci_GWAS, finemapped.loci_interval, genes_gtf, interval_eQTL_dataset = "Interval_full_eQTL")
  
}
# close pdf device
dev.off()

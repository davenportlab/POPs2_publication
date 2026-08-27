# 04_5_module_validation.R 

################################################################################

# 4.5. Validate relevant modules in Bar et al data

################################################################################

# Aim: Check for correlation of module eigengene values and measured lab tests from Bar et al 2025. 

########################### Output paths ##########################

########################### Input paths ###########################
# Read in coagulation and red blood cell measures from Bar et al 2025 github 
Bar_path <- "rna-seq/analysis/WGCNA/Bar_coagulation_tests/"
# sample info 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
# eigengenes
eigengenes_inpath <- "rna-seq/analysis/WGCNA/outputs/eigengenes_consensus_signed_hybrid.csv"

########################### Parameters ############################

########################### Load packages ###########################
library(tidyverse)

########################### Load data ###########################
# read in metadata  
sample.info <- read.csv(sample_info_inpath)
# eigengenes
eigengenes <- read.csv(eigengenes_inpath, row.names = 1)
# read in coagulation measures
Bar_coag <- as.data.frame(rbind(read.csv(paste0(Bar_path, "FIBRINOGEN.csv")) %>% mutate(test = "Fibrogen"),
                                read.csv(paste0(Bar_path, "PLT.csv")) %>% mutate(test = "Platelet count"),
                                read.csv(paste0(Bar_path, "PDW.csv")) %>% mutate(test = "Platelet distribution width"),
                                read.csv(paste0(Bar_path, "APTT_R.csv")) %>% mutate(test = "APTT ratio"),
                                read.csv(paste0(Bar_path, "MPV.csv")) %>% mutate(test = "Mean platelet volume"),
                                read.csv(paste0(Bar_path, "APTT_sec.csv")) %>% mutate(test = "APTT time (sec)"),
                                read.csv(paste0(Bar_path, "PT_SEC.csv")) %>% mutate(test = "Prothrombin time (sec)"),
                                read.csv(paste0(Bar_path, "PT_INR.csv")) %>% mutate(test = "Prothrombin INR")))

# read in RBC measures
Bar_RBC <- as.data.frame(rbind(read.csv(paste0(Bar_path, "HCT.csv")) %>% mutate(test = "Hematocrit"),
                               read.csv(paste0(Bar_path, "HGB.csv")) %>% mutate(test = "Hemoglobin"),
                               read.csv(paste0(Bar_path, "IRON.csv")) %>% mutate(test = "Iron"),
                               read.csv(paste0(Bar_path, "RBC.csv")) %>% mutate(test = "Red blood cell count"),
                               read.csv(paste0(Bar_path, "HEMOGLOBIN_A1C_CALCULATED.csv")) %>% mutate(test = "HbA1c"),
                               read.csv(paste0(Bar_path, "FERRITIN.csv")) %>% mutate(test = "Ferritin"),
                               read.csv(paste0(Bar_path, "FOLIC_ACID.csv")) %>% mutate(test = "Folic acid"),
                               read.csv(paste0(Bar_path, "HCT_HGB_RATIO.csv")) %>% mutate(test = "Hct/Hgb ratio"),
                               read.csv(paste0(Bar_path, "HDW.csv")) %>% mutate(test = "HDW"),
                               read.csv(paste0(Bar_path, "HYPERperc.csv")) %>% mutate(test = "Hyperchromic %"),
                               read.csv(paste0(Bar_path, "HYPO_perc.csv")) %>% mutate(test = "Hypochromic %"),
                               read.csv(paste0(Bar_path, "MACROperc.csv")) %>% mutate(test = "Macrocytic %"),
                               read.csv(paste0(Bar_path, "MCH.csv")) %>% mutate(test = "MCH"),
                               read.csv(paste0(Bar_path, "MCHC.csv")) %>% mutate(test = "MCHC"),
                               read.csv(paste0(Bar_path, "MCV.csv")) %>% mutate(test = "MCV"),
                               read.csv(paste0(Bar_path, "MICR_perc.csv")) %>% mutate(test = "Microcytic %"),
                               read.csv(paste0(Bar_path, "MICROperc_HYPOperc.csv")) %>% mutate(test = "Microcytic/Hypochromic ratio"),
                               read.csv(paste0(Bar_path, "RDW.csv")) %>% mutate(test = "RDW"),
                               read.csv(paste0(Bar_path, "RDW_SD.csv")) %>% mutate(test = "RDW SD"),
                               read.csv(paste0(Bar_path, "TRANSFERRIN.csv")) %>% mutate(test = "Transferrin"),
                               read.csv(paste0(Bar_path, "VITAMIN_B12.csv")) %>% mutate(test = "Vitiman B12")
))


########################### Analysis ###########################
# week: Week of the test relative to delivery (=0).
box_width <- 1

eigengene_for_comp <- eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id") %>%
  left_join(sample.info, by = "RNA_sanger_sample_id") %>%
  select(ME_10, ME_13, ME_2, ME_5, ME_9, ME_12, SampleGA, RNA_sanger_sample_id)

# forest plot of correlation between ME 10 and coagulation tests (Supplementary Figure 13A)
Bar_coag %>%
  mutate(
    week_start = as.numeric(str_extract(week, "-?\\d+")),
    week_end = as.numeric(str_extract(week, "(?<=,)-?\\d+")),
    week_mid = (week_start + week_end) / 2 + 0.5,
    gest_age = 40 + week_mid,
    test = factor(test, levels = c("APTT ratio", "APTT time (sec)", 
                                   "Platelet count", "Prothrombin INR", "Prothrombin time (sec)",
                                   "Mean platelet volume", "Platelet distribution width", "Fibrogen"))
  ) %>%
  select(gest_age, val_mean, test) %>%
  filter(gest_age >= min(sample.info$SampleGA), 
         gest_age <= max(sample.info$SampleGA)) %>%
  rename("val_mean" = "Bar_mean") %>%
  left_join(eigengene_for_comp %>%
              mutate(gest_age = floor(SampleGA)) %>%
              group_by(gest_age) %>%
              summarize(POPS_mean = mean(ME_10)), by = c("gest_age")) %>%
  group_by(test) %>%
  summarize(cor_test = list(cor.test(Bar_mean, POPS_mean, method = "pearson"))) %>%
  mutate(estimate = map_dbl(cor_test, ~ .x$estimate),
         conf_low = map_dbl(cor_test, ~ .x$conf.int[1]),
         conf_high = map_dbl(cor_test, ~ .x$conf.int[2]),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  ggplot(aes(x = estimate, y = reorder(test, estimate))) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(
    x = "Pearson correlation coefficient",
    y = "",
    title = "Correlation between coagulation tests and module 10"#,
    #caption = "Error bars = 95% CI"
  ) 

# forest plot of correlation between ME 13 and RBC lab tests (Supplementary Figure 13B)
Bar_RBC %>%
  mutate(
    week_start = as.numeric(str_extract(week, "-?\\d+")),
    week_end = as.numeric(str_extract(week, "(?<=,)-?\\d+")),
    week_mid = (week_start + week_end) / 2 + 0.5,
    gest_age = 40 + week_mid
  ) %>%
  select(gest_age, val_mean, test) %>%
  filter(gest_age >= min(sample.info$SampleGA), 
         gest_age <= max(sample.info$SampleGA)) %>%
  rename("val_mean" = "Bar_mean") %>%
  left_join(eigengene_for_comp %>%
              mutate(gest_age = floor(SampleGA)) %>%
              group_by(gest_age) %>%
              summarize(POPS_mean = mean(ME_13)), by = c("gest_age")) %>%
  group_by(test) %>%
  summarize(cor_test = list(cor.test(Bar_mean, POPS_mean, method = "pearson"))) %>%
  mutate(estimate = map_dbl(cor_test, ~ .x$estimate),
         conf_low = map_dbl(cor_test, ~ .x$conf.int[1]),
         conf_high = map_dbl(cor_test, ~ .x$conf.int[2]),
         p_value = map_dbl(cor_test, ~ .x$p.value)) %>%
  ggplot(aes(x = estimate, y = reorder(test, estimate))) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(
    x = "Pearson correlation coefficient",
    y = "",
    title = "Correlation between RBC tests and module 13"#,
    #caption = "Error bars = 95% CI"
  ) + 
  expand_limits(x = c(-1, 1)) +
  scale_x_continuous(#expand = c(0.02,0),
    breaks = c(-1, -0.5, 0, 0.5, 1))


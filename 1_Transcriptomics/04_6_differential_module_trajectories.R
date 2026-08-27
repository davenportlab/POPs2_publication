# 04_6_differential_module_trajectories.R 

################################################################################

# 4.6. Differential module trajectories across clinical outcome groups

################################################################################

# Aim: Identify modules with differential trajectories acorss clinical subgroups

########################### Output paths ##########################
outfiles_path <- "rna-seq/analysis/WGCNA/outputs/outcome_trajectories/"

########################### Input paths ###########################
# sample info 
sample_info_inpath = "rna-seq/data/sample_covariates_clin_tech.csv"
# eigengenes
eigengenes_inpath <- "rna-seq/analysis/WGCNA/outputs/eigengenes_consensus_signed_hybrid.csv"

########################### Parameters ############################

########################### Load packages ###########################
library(lme4)
library(tidyverse)
library(ggeffects)

# set ggplot theme 
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             text=element_text(size=15)))

########################### Load data ###########################
# make sample info df
sample.info <- read.csv(sample_info_inpath)
# eigengenes
eigengenes <- read.csv(eigengenes_inpath, row.names = 1)

########################### Analysis ###########################
# sample info 
sample.info <- sample.info %>%
  mutate(pn_gdm = as.factor(pn_diabetes_3cat == 2), 
         sga = as.factor(BW_Centile_Br1990 < 10),
         pn_ptd_sp = as.factor(pn_ptd_sp == 1))

################################################################################
# GDM, 30 individuals, 0 module with differential trajectories 
################################################################################
# merge sample info to eigengene info
eigen.time.data_gdm <- sample.info %>% 
  select(RNA_sanger_sample_id, SampleGA, id_run_position, an_eth_cat_str, pn_gdm, ANON_ID, Sample_taken_at) %>%
  filter(!is.na(pn_gdm)) %>%
  left_join(eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id"), by = "RNA_sanger_sample_id") %>%
  mutate(ANON_ID = as.factor(ANON_ID), 
         Sample_taken_at = as.factor(Sample_taken_at)) 
# Results storage
lme.results_gdm <- list()
lme_prediction_df_gdm <- c()
lme_interaction_signif_gdm <- c()
# loop through modules 
for (eigengene in colnames(eigengenes)) {
  print(eigengene)
  res <- lmerTest::lmer(as.formula(get(eigengene) ~ splines::bs(x = SampleGA, degree = 3) * pn_gdm + (1|ANON_ID) + (1|id_run_position)), # the way the interaction term is written with * means it will create separate terms also for splines::bs(x = SampleGA, degree = 3) and pn_gdm as well as splines::bs(x = SampleGA, degree = 3) * pn_gdm. If it were written with : i.e. splines::bs(x = SampleGA, degree = 3):pn_gdm, it would only include that term. So we are happy with this. 
                        data = eigen.time.data_gdm, 
                        REML = FALSE) # check if this is the correct format 
  
  lme.results_gdm[[eigengene]] <- res
  
  # based on the model, make a prediction 
  # use ggemmeans to estimate confidence intervals/fits
  prediction <- as.data.frame(ggemmeans(lme.results_gdm[[eigengene]], c("SampleGA [all]", "pn_gdm"))) %>% mutate(Mod = eigengene)
  colnames(prediction) <- c("SampleGA", "pred_eig_val", "std.err", "conf.low", "conf.high", "group", "Module") 
  lme_prediction_df_gdm <- rbind(lme_prediction_df_gdm, prediction)
  
  # extract whether the interaction is significant 
  interaction_signif <- as.data.frame(t(c(anova(lme.results_gdm[[eigengene]])[,6], eigengene)))
  colnames(interaction_signif) <- c("time_pval", "group_pval", "time.group.interaction_pval", "Module")
  lme_interaction_signif_gdm <- rbind(lme_interaction_signif_gdm, interaction_signif)
}
# plot significant interactions
if(nrow(lme_interaction_signif_gdm %>% filter(time.group.interaction_pval < 0.05)) > 0){
  print("There are significant interactions with gestational diabetes")
  lme_prediction_df_gdm %>% 
    filter(Module %in% (lme_interaction_signif_gdm %>% filter(time.group.interaction_pval < 0.05))$Module) %>%
    # join to interaction p vals 
    mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module)))) %>%
    ggplot(aes(x = SampleGA, group = group, fill = group)) + 
    geom_line(aes(x = SampleGA, y = pred_eig_val, color = group), size = 1.7) +
    geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.055) + 
    ylab("Module eigengene value") + xlab("Sample gestational age (weeks)") + 
    facet_wrap(~Module, scales = "free", ncol = 5) + 
    ggtitle("WGCNA module eigengene values modeled over time\nSigned hybrid consensus modules\nGestational diabetes differential continuous trajectory") 
}else{
  print("There are no significant interactions with gestational diabetes")
}

################################################################################
# Pre eclampsia, 50 individuals 1 modules with differential trajectories  
################################################################################
n_ind_PE <- sample.info %>% 
  select(pn_petACOG13, ANON_ID) %>%
  distinct() %>% 
  filter(pn_petACOG13 == 1) %>% nrow()
# Model
# Merge sample info to eigengene info
eigen.time.data_pe <- sample.info %>% 
  select(RNA_sanger_sample_id, SampleGA, id_run_position, an_eth_cat_str, pn_petACOG13, ANON_ID, Sample_taken_at) %>%
  filter(!is.na(pn_petACOG13)) %>%
  left_join(eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id"), by = "RNA_sanger_sample_id") %>%
  mutate(ANON_ID = as.factor(ANON_ID), 
         Sample_taken_at = as.factor(Sample_taken_at)) 
# Results storage
lme.results_pe <- list()
lme_prediction_df_pe <- c()
lme_interaction_signif_pe <- c()

for (eigengene in colnames(eigengenes)) {
  print(eigengene)
  res <- lmerTest::lmer(as.formula(get(eigengene) ~ splines::bs(x = SampleGA, degree = 3) * pn_petACOG13 + (1|ANON_ID) + (1|id_run_position)), 
                        data = eigen.time.data_pe, 
                        REML = FALSE) # check if this is the correct format 
  
  lme.results_pe[[eigengene]] <- res
  
  # based on the model, make a prediction 
  # use ggemmeans to estimate confidence intervals/fits
  prediction <- as.data.frame(ggemmeans(lme.results_pe[[eigengene]], c("SampleGA [all]", "pn_petACOG13"))) %>% mutate(Mod = eigengene)
  colnames(prediction) <- c("SampleGA", "pred_eig_val", "std.err", "conf.low", "conf.high", "group", "Module") 
  lme_prediction_df_pe <- rbind(lme_prediction_df_pe, prediction)
  
  # extract whether the interaction is significant 
  interaction_signif <- as.data.frame(t(c(anova(lme.results_pe[[eigengene]])[,6], eigengene)))
  colnames(interaction_signif) <- c("time_pval", "group_pval", "time.group.interaction_pval", "Module")
  lme_interaction_signif_pe <- rbind(lme_interaction_signif_pe, interaction_signif)
}

#Plot significant interactions
if(nrow((lme_interaction_signif_pe %>% filter(time.group.interaction_pval < 0.05))) > 0){
  print("There are significant interactions with pre-eclampsia")
  lme_prediction_df_pe %>% 
    filter(Module %in% (lme_interaction_signif_pe %>% filter(time.group.interaction_pval < 0.05))$Module) %>%
    # join to interaction p vals 
    mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module))), 
           # edit the group
           `Pre-eclampsia` = ifelse(group == 1, TRUE, FALSE)) %>%
    ggplot(aes(x = SampleGA, group = `Pre-eclampsia`, fill = `Pre-eclampsia`)) + 
    scale_fill_manual(values = c("#73A2C6", "#F4777F")) + 
    scale_color_manual(values = c("#73A2C6", "#F4777F")) + 
    geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Pre-eclampsia`), size = 1.7) +
    geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Pre-eclampsia`), alpha = 0.055) + 
    ylab("Module eigengene value") + xlab("Sample gestational age (weeks)") + 
    facet_wrap(~Module, scales = "free", ncol = 3) + 
    theme(strip.background = element_blank()) + 
    ggtitle(paste0("WGCNA module eigengene values modeled over time\nPre-eclampsia (n=", n_ind_PE, ") differential continuous trajectory")) 
  
} else{
  print("There are no significant interactions with pre-eclampsia")
}

################################################################################
# Small for gestational age (SGA), 50 individuals 2 modules with differential trajectories
################################################################################
n_ind_SGA <- sample.info %>% 
  select(sga, ANON_ID) %>%
  distinct() %>% 
  filter(sga == TRUE) %>% nrow()
# Model
# merge sample info to eigengene info
eigen.time.data_sga <- sample.info %>% 
  select(RNA_sanger_sample_id, SampleGA, id_run_position, an_eth_cat_str, sga, ANON_ID, Sample_taken_at) %>%
  filter(!is.na(sga)) %>%
  left_join(eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id"), by = "RNA_sanger_sample_id") %>%
  mutate(ANON_ID = as.factor(ANON_ID), 
         Sample_taken_at = as.factor(Sample_taken_at),
         sga = as.numeric(sga)) 

# Results storage
lme.results_sga <- list()
lme_prediction_df_sga <- c()
lme_interaction_signif_sga <- c()

for (eigengene in colnames(eigengenes)) {
  print(eigengene)
  res <- lmerTest::lmer(as.formula(get(eigengene) ~ splines::bs(x = SampleGA, degree = 3) * sga + (1|ANON_ID) + (1|id_run_position)), 
                        data = eigen.time.data_sga, 
                        REML = FALSE) # check if this is the correct format 
  
  lme.results_sga[[eigengene]] <- res
  
  # based on the model, make a prediction 
  # use ggemmeans to estimate confidence intervals/fits
  prediction <- as.data.frame(ggemmeans(lme.results_sga[[eigengene]], c("SampleGA [all]", "sga"))) %>% mutate(Mod = eigengene)
  colnames(prediction) <- c("SampleGA", "pred_eig_val", "std.err", "conf.low", "conf.high", "group", "Module") 
  lme_prediction_df_sga <- rbind(lme_prediction_df_sga, prediction)
  
  # extract whether the interaction is significant 
  interaction_signif <- as.data.frame(t(c(anova(lme.results_sga[[eigengene]])[,6], eigengene)))
  colnames(interaction_signif) <- c("time_pval", "group_pval", "time.group.interaction_pval", "Module")
  lme_interaction_signif_sga <- rbind(lme_interaction_signif_sga, interaction_signif)
  
}

# Plot significant interactions
if(nrow((lme_interaction_signif_sga %>% filter(time.group.interaction_pval < 0.05))) > 0){
  print("There are significant interactions with small for gestational age")
  lme_prediction_df_sga %>% 
    filter(Module %in% (lme_interaction_signif_sga %>% filter(time.group.interaction_pval < 0.05))$Module) %>%
    # join to interaction p vals 
    mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module))),
           `Small for gestational age` = ifelse(group == 1, TRUE, FALSE)) %>%
    ggplot(aes(x = SampleGA, group = `Small for gestational age`, fill = `Small for gestational age`)) + 
    geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Small for gestational age`), size = 1.7) +
    geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Small for gestational age`), alpha = 0.055) + 
    ylab("Module eigengene value") + xlab("Sample gestational age (weeks)") + 
    facet_wrap(~Module, scales = "free", ncol = 3) + 
    theme(strip.background = element_blank()) + 
    scale_fill_manual(values = c("#73A2C6", "#F4777F")) + 
    scale_color_manual(values = c("#73A2C6", "#F4777F")) + 
    ggtitle(paste0("WGCNA module eigengene values modeled over time\nSGA (n=", n_ind_SGA, ") differential continuous trajectory")) 
  
} else{
  print("There are no significant interactions with SGA")
}

################################################################################
# Spontaneous preterm delivery (SPD) 18 individuals, 6 differential trajectories 
################################################################################
n_ind_sptd <- sample.info %>% 
  select(pn_ptd_sp, ANON_ID) %>%
  distinct() %>% 
  filter(pn_ptd_sp == 1) %>% nrow()
# Model
# Merge sample info to eigengene info
eigen.time.data_sptd <- sample.info %>% 
  select(RNA_sanger_sample_id, SampleGA, id_run_position, an_eth_cat_str, pn_ptd_sp, ANON_ID, Sample_taken_at)%>%
  filter(!is.na(pn_ptd_sp)) %>%
  left_join(eigengenes %>% rownames_to_column(var = "RNA_sanger_sample_id"), by = "RNA_sanger_sample_id") %>%
  mutate(ANON_ID = as.factor(ANON_ID), 
         Sample_taken_at = as.factor(Sample_taken_at),
         pn_ptd_sp = as.numeric(pn_ptd_sp)) 

# Results storage
lme.results_sptd <- list()
lme_prediction_df_sptd <- c()
lme_interaction_signif_sptd <- c()

for (eigengene in colnames(eigengenes)) {
  print(eigengene)
  res <- lmerTest::lmer(as.formula(get(eigengene) ~ splines::bs(x = SampleGA, degree = 3) * pn_ptd_sp + (1|ANON_ID) + (1|id_run_position)), 
                        data = eigen.time.data_sptd,
                        REML = FALSE) # check if this is the correct format 
  
  lme.results_sptd[[eigengene]] <- res
  
  # based on the model, make a prediction 
  # use ggemmeans to estimate confidence intervals/fits
  prediction <- as.data.frame(ggemmeans(lme.results_sptd[[eigengene]], c("SampleGA [all]", "pn_ptd_sp"))) %>% mutate(Mod = eigengene)
  colnames(prediction) <- c("SampleGA", "pred_eig_val", "std.err", "conf.low", "conf.high", "group", "Module") 
  lme_prediction_df_sptd <- rbind(lme_prediction_df_sptd, prediction)
  
  # extract whether the interaction is significant 
  interaction_signif <- as.data.frame(t(c(anova(lme.results_sptd[[eigengene]])[,6], eigengene)))
  colnames(interaction_signif) <- c("time_pval", "group_pval", "time.group.interaction_pval", "Module")
  lme_interaction_signif_sptd <- rbind(lme_interaction_signif_sptd, interaction_signif)
  
}

#  Plot significant interactions
if(nrow(lme_interaction_signif_sptd %>% filter(time.group.interaction_pval < 0.05)) > 1){
  print("There are significant interactions with spontaneous preterm delivery")
  lme_prediction_df_sptd %>% 
    filter(Module %in% (lme_interaction_signif_sptd %>% filter(time.group.interaction_pval < 0.05))$Module) %>%
    # join to interaction p vals 
    mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module))),
           `Spontaneous preterm delivery` = ifelse(group == 1, TRUE, FALSE)) %>%
    ggplot(aes(x = SampleGA, group = `Spontaneous preterm delivery`, fill = `Spontaneous preterm delivery`)) + 
    geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Spontaneous preterm delivery`), size = 1.7) +
    geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Spontaneous preterm delivery`), alpha = 0.055) + 
    ylab("Module eigengene value") + xlab("Sample gestational age (weeks)") + 
    facet_wrap(~Module, scales = "free", ncol = 3) + 
    theme(strip.background = element_blank()) + 
    scale_fill_manual(values = c("#F4777F", "#73A2C6")) + 
    scale_color_manual(values = c("#F4777F", "#73A2C6")) + 
    ggtitle(paste0("WGCNA module eigengene values modeled over time\nSpontaneous preterm delivery (n=", n_ind_sptd, ") differential continuous trajectory")) 
}else{
  print("There are no significant interactions with small for gestational age")
}

# PLOT

# make module titles
custom_labels <- setNames(c("1: Innate immune response", 
                            "2: Cell proliferation", 
                            "3: Cell motility", 
                            "4: T cell immune response", 
                            "5: Protein synthesis/translation",
                            "6: B cell immune response",
                            "7: Eosinophil allergic-like response",
                            "8: Cytotoxic immune response",
                            "9: Protein sythesis/translation",
                            "10: Coagulation/clotting",
                            "11: Histones",
                            "12: Mitochondria/energy production",
                            "13: Small molecule transport",
                            "14: Antiviral interferon response",
                            "15: Inflammatory response regulation", 
                            "Unassigned"), 
                          1:15)

# pe 
pe_plot_df <- lme_prediction_df_pe %>%
  left_join(
    lme_interaction_signif_pe %>%
      select(Module, time.group.interaction_pval),
    by = "Module"
  ) %>%
  filter(time.group.interaction_pval < 0.05) %>%
  mutate(
    Module = as.factor(as.numeric(gsub("ME_", "", Module))),
    `Pre-eclampsia` = ifelse(group == 1, TRUE, FALSE)
  )

pe_label_df <- pe_plot_df %>%
  group_by(Module) %>%
  summarise(
    x = max(SampleGA),
    y = max(conf.high),
    p = first(time.group.interaction_pval),
    .groups = "drop"
  ) %>%
  mutate(
    p = as.numeric(p),
    label = paste0("pint = ", signif(p, 3))
  )

pe_mod_plt <- pe_plot_df %>%
  ggplot(aes(x = SampleGA, group = `Pre-eclampsia`, fill = `Pre-eclampsia`)) + 
  scale_fill_manual(values = c("#73A2C6", "#F4777F")) + 
  scale_color_manual(values = c("#73A2C6", "#F4777F")) + 
  geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Pre-eclampsia`), size = 1.7) +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Pre-eclampsia`), alpha = 0.055) + 
  ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
  facet_wrap(~Module, scales = "free", ncol = 3, labeller = labeller(Module = custom_labels)) + 
  theme(legend.position = "inside", 
        legend.position.inside =  c(0.5, 0.9),
        legend.direction = "horizontal",
        axis.text.y=element_blank(),
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) + 
  ggtitle(paste0("Modules with different trajectory")) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  geom_text(
    data = pe_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  )

pe_mod_plt_paper_sup <- pe_plot_df %>%
  ggplot(aes(x = SampleGA, group = `Pre-eclampsia`, fill = `Pre-eclampsia`)) + 
  scale_fill_manual(values = c("#73A2C6", "#F4777F")) + 
  scale_color_manual(values = c("#73A2C6", "#F4777F")) + 
  geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Pre-eclampsia`), size = 1.7) +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Pre-eclampsia`), alpha = 0.055) + 
  ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
  facet_wrap(~Module, scales = "free", ncol = 3, labeller = labeller(Module = custom_labels)) + 
  theme(legend.direction = "horizontal",
        legend.position = "top",
        axis.text.y=element_blank(),
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) + 
  #ggtitle(paste0("Modules with different trajectory")) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  geom_text(
    data = pe_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  )


# sga
sga_plot_df <- lme_prediction_df_sga %>%
  left_join(
    lme_interaction_signif_sga %>%
      select(Module, time.group.interaction_pval),
    by = "Module"
  ) %>%
  filter(time.group.interaction_pval < 0.05) %>%
  mutate(
    Module = as.factor(as.numeric(gsub("ME_", "", Module))),
    SGA = group == 1
  )

sga_label_df <- sga_plot_df %>%
  group_by(Module) %>%
  summarise(
    x = max(SampleGA),
    y = max(conf.high),
    p = first(time.group.interaction_pval),
    .groups = "drop"
  ) %>%
  mutate(
    p = as.numeric(p),
    label = paste0("pint = ", signif(p, 2))
  ) 

sga_mod_plt <- sga_plot_df %>%
  ggplot(aes(x = SampleGA, group = SGA, fill = SGA)) + 
  geom_line(aes(x = SampleGA, y = pred_eig_val, color = SGA), size = 1.7) +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = SGA), alpha = 0.055) + 
  ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
  facet_wrap(~Module, scales = "free", ncol = 3, labeller = labeller(Module = custom_labels)) + 
  theme(legend.position = "top",
        axis.text.y=element_blank(),
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) + 
  scale_fill_manual(values = c("#73A2C6", "#F4777F")) + 
  scale_color_manual(values = c("#73A2C6", "#F4777F")) + 
  ggtitle(paste0("Modules with differential trajectory")) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  geom_text(
    data = sga_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  )

sga_mod_plt_paper_sup <- sga_plot_df %>%
  ggplot(aes(x = SampleGA, group = SGA, fill = SGA)) + 
  geom_line(aes(x = SampleGA, y = pred_eig_val, color = SGA), size = 1.7) +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = SGA), alpha = 0.055) + 
  ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
  facet_wrap(~Module, scales = "free", ncol = 3, labeller = labeller(Module = custom_labels)) + 
  theme(legend.position = "top",
        axis.text.y=element_blank(),
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) + 
  scale_fill_manual(values = c("#73A2C6", "#F4777F"), 
                    name = "Small for gestational age") + 
  scale_color_manual(values = c("#73A2C6", "#F4777F"), 
                     name = "Small for gestational age") + 
  #ggtitle(paste0("Modules with differential trajectory")) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  geom_text(
    data = sga_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  )

# sptd 
sptd_plot_df <- lme_prediction_df_sptd %>%
  left_join(
    lme_interaction_signif_sptd %>%
      select(Module, time.group.interaction_pval),
    by = "Module"
  ) %>%
  filter(time.group.interaction_pval < 0.05) %>%
  mutate(
    Module = as.factor(as.numeric(gsub("ME_", "", Module))),
    `Spontaneous preterm delivery` = ifelse(group == 1, TRUE, FALSE)
  )

sptd_label_df <- sptd_plot_df %>%
  group_by(Module) %>%
  summarise(
    x = max(SampleGA),
    y = max(conf.high),
    p = first(time.group.interaction_pval),
    .groups = "drop"
  ) %>%
  mutate(
    p = as.numeric(p),
    label = paste0("pint = ", signif(p, 2))
  ) 

sptd_mod_plt <- sptd_plot_df %>%
  ggplot(aes(x = SampleGA, group = `Spontaneous preterm delivery`, fill = `Spontaneous preterm delivery`)) + 
  geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Spontaneous preterm delivery`), size = 1.7) +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Spontaneous preterm delivery`), alpha = 0.055) + 
  ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
  facet_wrap(~Module, scales = "free_y", ncol = 1, labeller = labeller(Module = custom_labels)) + 
  theme(legend.position = "top",
        axis.text.y=element_blank(),
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 12)) + 
  scale_fill_manual(values = c("#73A2C6", "#F4777F"),
                    name = "sPTB") + 
  scale_color_manual(values = c("#73A2C6", "#F4777F"),
                     name = "sPTB") + 
  ggtitle(paste0("Modules w different trajectory")) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  geom_text(
    data = sptd_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  )

sptd_mod_plt_paper_sup <- sptd_plot_df %>%
  ggplot(aes(x = SampleGA, group = `Spontaneous preterm delivery`, fill = `Spontaneous preterm delivery`)) + 
  geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Spontaneous preterm delivery`), size = 1.7) +
  geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Spontaneous preterm delivery`), alpha = 0.055) + 
  ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
  facet_wrap(~Module, scales = "free_y", nrow = 2, labeller = labeller(Module = custom_labels)) + 
  theme(legend.position = "top",
        axis.text.y=element_blank(),
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 12)) + 
  scale_fill_manual(values = c("#73A2C6", "#F4777F"),
                    name = "Spontaneous preterm birth") + 
  scale_color_manual(values = c("#73A2C6", "#F4777F"),
                     name = "Spontaneous preterm birth") + 
  #ggtitle(paste0("Modules w different trajectory")) + 
  scale_x_continuous(breaks = c(12, 20, 28, 36)) + 
  geom_text(
    data = sptd_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  )

# sptd_mod_single_main <- lme_prediction_df_sptd %>% 
#   filter(Module == "ME_6") %>%
#   # join to interaction p vals 
#   mutate(Module = as.factor(as.numeric(gsub("ME_", "", Module))),
#          `Spontaneous preterm delivery` = ifelse(group == 1, TRUE, FALSE)) %>%
#   ggplot(aes(x = SampleGA, group = `Spontaneous preterm delivery`, fill = `Spontaneous preterm delivery`)) + 
#   geom_line(aes(x = SampleGA, y = pred_eig_val, color = `Spontaneous preterm delivery`), size = 1.7) +
#   geom_ribbon(aes(SampleGA, pred_eig_val, ymin = conf.low, ymax = conf.high, fill = `Spontaneous preterm delivery`), alpha = 0.055) + 
#   ylab("Module eigengene value") + xlab("Gestational age (weeks)") + 
#   facet_wrap(~Module, scales = "free_y", ncol = 1, labeller = labeller(Module = custom_labels)) + 
#   theme(
#         axis.text.y=element_blank(),
#         plot.title = element_text(size = 15),
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15),
#         strip.text = element_text(size = 15)) + 
#   scale_fill_manual(values = c("#73A2C6", "#F4777F"),
#                     name = "sPTB") + 
#   scale_color_manual(values = c("#73A2C6", "#F4777F"),
#                      name = "sPTB") + 
#   ggtitle(paste0("Module with different trajectory in sPTB")) + 
#   scale_x_continuous(breaks = c(12, 20, 28, 36))
# ggsave(paste0(plt_outpath, 'sPTB_mod6_main.pdf'), sptd_mod_single_main, width = 5, height = 4, useDingbats = FALSE)


# paste all together to make sup figure for paper
all_mod_traj_sup <- gridExtra::grid.arrange(pe_mod_plt_paper_sup + labs(tag = "A"), 
                                            sga_mod_plt_paper_sup + labs(tag = "B"), 
                                            sptd_mod_plt_paper_sup + labs(tag = "C"), 
                                            ncol = 2, nrow = 2, 
                                            widths = c(1, 2), 
                                            heights = c(1.25, 2),
                                            layout_matrix = rbind(c(1, 2),
                                                                  c(3, 3)))
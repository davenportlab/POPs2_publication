# 02_4_plot_cibersort.R 

################################################################################

# 2.4. Plot imputed cell proportions

################################################################################

# Aim: Plot cibersort outputs 

########################### Input paths ###########################
cibersort_inpath <- "rna-seq/analysis/cibersort/outputs/CIBERSORTx_Adjusted.txt"
sample_info_inpath <- "rna-seq/data/sample_covariates_clin_tech.csv"

########################### Load packages ###########################
library(tidyverse)
library(ggrastr)

# set ggplot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Load data ###########################
# read in cibersort output imputed cell proportions
cell_prop <- read.delim(cibersort_inpath)
# read in sample info 
sample.info <- read.csv(sample_info_inpath) 

########################### Analysis ###########################
# Join Metadata patients and time points with the cell proportions
cell_prop <- cell_prop %>% left_join(sample.info %>% dplyr::select(RNA_sanger_sample_id, Sample_taken_at, ANON_ID),  by=c("Mixture" ="RNA_sanger_sample_id")) 
cell_prop2 <- cell_prop
# remove dots from cell type colnames
colnames(cell_prop) <- gsub("\\.", " ", colnames(cell_prop))
colnames(cell_prop) <- gsub("  ", " ", colnames(cell_prop)) 

########################### Plot ###########################

# 1. Stacked bar chart (Supplementary figure 2)

# make distinct color pallette 
col_vector <- rev(c("#0075DC",  "#FF5000", "#2BCE48", "#C20088",  "#E0FF66", "#FFA8BE", "#100AFF", "#FF0010", "#00998F", "#F0A3FF", "#FFFF80", "#9DCC00", "#990000", "#FFE100", "#5EF1F2", "#FDB462", "#4C005C", "#E78AC3", "#FFFFCC", "#426600", "#993F00", "#94FFB5"))
# plot 
cell_prop %>% 
  mutate(Sample_taken_at = gsub("_", " ", Sample_taken_at)) %>%
  # remove columns we don't want to plot 
  dplyr::select(-`P value`, -Correlation, -RMSE) %>% 
  # arrange samples based on neutrophil abundance 
  arrange(Neutrophils) %>% 
  mutate(order = c(1:nrow(cell_prop))) %>%
  gather(Cell_Type, Proportion, `B cells naive`:`Neutrophils`, factor_key=TRUE) %>% 
  group_by(Cell_Type) %>% mutate(mean = mean(Proportion)) %>% ungroup() %>% arrange((mean)) %>% mutate(`Cell type` = factor(Cell_Type, levels = unique(Cell_Type))) %>%
  ggplot(aes(fill=`Cell type`, y=Proportion, x=reorder(Mixture, order))) + 
  facet_grid(~ Sample_taken_at, scale="free_x",space="free_x") + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        text=element_text(size=15),
        legend.position = "bottom") + 
  xlab("Sample") + 
  ggtitle("Relative abundance of CIBERSORTx imputed cell proportions by sample time-point") + 
  scale_fill_manual(values=col_vector) + 
  guides(fill = guide_legend(reverse = TRUE))

# 2. Boxplot (Supplementary Figure 3)
cell_prop %>% mutate(`Sample time-point (weeks)` = gsub("_weeks", "", Sample_taken_at)) %>% 
  dplyr::select(-`P value`, -Correlation, -RMSE, -ANON_ID, -Sample_taken_at) %>%  
  pivot_longer(cols = !c(Mixture, `Sample time-point (weeks)`), names_to = "Cell Type", values_to = "Proportion") %>% 
  ggplot(aes(x = `Sample time-point (weeks)`, y = Proportion)) +
  geom_jitter(aes(color = `Sample time-point (weeks)`), 
              size = 0.01, position = position_jitter(width = 0.25, height = 0)) + 
  geom_boxplot(fill = "transparent", outlier.alpha = 0) +  
  scale_color_manual(values=c("#BCE4D8", "#83C4CB", "#439FB7", "#32769B")) +
  #geom_boxplot(outlier.alpha = 0.4)  + 
  facet_wrap(~`Cell Type`, scales = "free", ncol = 4) + 
  theme(legend.position = "none") +
  ggtitle("Boxplot of distribution of CIBERSORTx imputed cell proportions by sample time-point")

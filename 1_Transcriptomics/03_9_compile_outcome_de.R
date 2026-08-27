# 03_9_compile_outcome_de.R 

################################################################################

# 3.9. Compile outcome differential expression results

################################################################################

# Aim: compile and plot results of outcome DE

########################### Output paths ##########################
sig_toptable_outpath <- "rna-seq/analysis/outcome_pred/outputs/outcome_all_sig_toptable.csv"

########################### Input paths ###########################
# DE dirs
dir_alltp <- "rna-seq/analysis/outcome_pred/outputs/de_all_tps/"
dir_12wk <- "rna-seq/analysis/outcome_pred/outputs/de_12wk/"
dir_20wk <- "rna-seq/analysis/outcome_pred/outputs/de_20wk/"
dir_28wk <- "rna-seq/analysis/outcome_pred/outputs/de_28wk/"
dir_36wk <- "rna-seq/analysis/outcome_pred/outputs/de_36wk/"

########################### Parameters ############################
# define deg cutoffs
p_val = 0.05
log_FC = log2(1.5)

########################### Load packages ###########################
library(tidyverse)
library(edgeR)

########################### Load data ###########################

########################### Analysis ###########################
# make function to read in DE and make top tables
make_toptable <- function(path, outcome, tp, comp){
  print(paste0(path, outcome, "_", tp, "dream_fit.rds"))
  # read in dream fit
  fit <- readRDS(paste0(path, outcome, "_", tp, "dream_fit.rds"))
  # get top table
  toptable <- topTable(fit, coef = comp, number = 18826)
  #return(toptable)
}
# list all that need to be read in 
DE_list <- as.data.frame(cbind(outcome = c(rep("GDM", 5), rep("PE", 5), rep("SGA", 5), rep("PTD", 5)),
                               comp = c(rep("compare_GDM_vs_noGDM", 5), 
                                        rep("compare_PE_vs_noPE", 5), 
                                        rep("compare_SGA_vs_noSGA", 5), 
                                        rep("compare_PTD_vs_noPTD", 5)),
                               tp = c("", "12wk_", "20wk_", "28wk_", "36wk_"),
                               dir = c(dir_alltp, dir_12wk, dir_20wk, dir_28wk, dir_36wk))) %>%
  mutate(tp_str = ifelse(tp == "", "all_tps", gsub("_", "", tp)))

long_toptable <- c()
# loop through all toptables 
for(i in 1:nrow(DE_list)){
  toptable <- make_toptable(path = DE_list$dir[i],
                            outcome = DE_list$outcome[i],
                            comp = DE_list$comp[i],
                            tp = DE_list$tp[i])
  long_toptable <- rbind(long_toptable, 
                         toptable %>% 
                           rownames_to_column(var = "gtf.gene_id") %>%
                           mutate(DE = paste0(DE_list$outcome[i], "_", DE_list$tp_str[i]),
                                  Outcome = DE_list$outcome[i],
                                  Timepoint = DE_list$tp_str[i]))
}

dim(long_toptable)
head(long_toptable)

all_sig_toptable <- long_toptable %>% filter((abs(logFC) > log_FC) & (adj.P.Val < p_val))
dim(all_sig_toptable)
head(all_sig_toptable)
# write out all sig outcome DEGs to one table (Supplementary table)
all_sig_toptable %>% write.csv(sig_toptable_outpath, row.names = FALSE)


########################### Plot ###########################
# make annotation to add n DEGs
n_degs <- long_toptable %>% 
  filter(abs(logFC) > log_FC, 
         adj.P.Val < p_val) %>% 
  group_by(DE) %>%
  summarize(n = n())

annot_df <- long_toptable %>% dplyr::select(DE, Outcome, Timepoint) %>% distinct() %>%
  left_join(n_degs, by = "DE") %>%
  mutate(n = ifelse(is.na(n), 0, n),
         xpos = ifelse(Outcome == "GDM", -1.888353,
                       ifelse(Outcome == "PE", -1.845137,
                              ifelse(Outcome == "SGA", -1.448821, -2.980605))),
         ypos = ifelse(Timepoint == "12wk", 2.23005, 
                       ifelse(Timepoint == "20wk", 1.369143,
                              ifelse(Timepoint == "28wk", 4.890148, 
                                     ifelse(Timepoint == "36wk", 2.790472, 1.910223)))),
         Timepoint = gsub("wk", " weeks", Timepoint),
         Timepoint = ifelse(Timepoint == "all_tps", "All time-points", Timepoint),
         Outcome = ifelse(Outcome == "PTD", "sPTB", Outcome),
         Outcome = factor(Outcome, levels = c("GDM", "PE", "SGA", "sPTB")),
         n = paste0("n DEGs = ", n)) 

# plot
long_toptable %>% 
  mutate(reg = ifelse(((logFC > log_FC) & (adj.P.Val < p_val)), "UP", 
                      ifelse(((logFC < -log_FC) & (adj.P.Val < p_val)), "DOWN", "no_change")), 
         adj.P.Val = ifelse(adj.P.Val < .Machine$double.xmin, .Machine$double.xmin, adj.P.Val),
         Timepoint = gsub("wk", " weeks", Timepoint),
         Timepoint = ifelse(Timepoint == "all_tps", "All time-points", Timepoint),
         Outcome = ifelse(Outcome == "PTD", "sPTB", Outcome),
         Outcome = factor(Outcome, levels = c("GDM", "PE", "SGA", "sPTB"))
  ) %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val), 
             color = reg)) +
  scale_color_manual(name = "reg",
                     values =c("no_change"= "#999999", 
                               "UP" = "red", 
                               "DOWN" = "red")) +
  geom_point(shape = 19, show.legend = FALSE, size = 0.2) +
  geom_hline(yintercept = -log10(p_val), linetype="longdash", colour="grey", size=0.5) +
  geom_vline(xintercept = log_FC, linetype="longdash", colour="#999999", size=0.5) +
  geom_vline(xintercept = -log_FC, linetype="longdash", colour="#999999", size=0.5) + 
  labs(title="Volcano plots of DEGs across outcome groups at different time-points",
       subtitle = paste("Fold-change > 1.5, adj.P.Val <", p_val)) + 
  theme(legend.position="none", panel.grid = element_blank(), strip.background = element_blank()) + 
  coord_cartesian(clip = "off") + 
  xlab("Log 2 fold change") + 
  ylab("-log10(adjusted p-value)") + 
  facet_grid(Timepoint~Outcome, scales = "free") + 
  geom_text(
    data = annot_df,
    aes(
      x = xpos,
      y = ypos,
      label = n
    ),
    size = 3,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1
  )





# 05_3_varying_peer_factors_single_tp.R

################################################################################

# 5.3. Compare eQTL results after including increasing number of PEER factors in the model

################################################################################

# Aim: Decide n peer factors to include in main effects eQTL model for pops2 TP

########################### Output paths ##########################

########################### Input paths ###########################
# path to input files
inpath <- "genotyping/analysis/eQTL/output_data/single_tp_analysis/peertests/"

########################### Parameters ############################

########################### Load packages ###########################
library(data.table)
library(ggplot2)
options(stringsAsFactors = FALSE)

########################### Load data ###########################

########################### Analysis ###########################

# Function to read in set of eQTL results
cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, full.names = T)
  print(paste("Number of files:", length(temp_eqtl)))
  eqtl <- lapply(temp_eqtl, readRDS)
  eqtl <- rbindlist(eqtl)
  print(c("eQTL results:", dim(eqtl)))
  saveRDS(eqtl, output_filename)
}

n.peer <- seq(0, 50, 5)

# Consolidate eQTL files from one per gene to one per mapping 
################################################################################
# 12 weeks only 
################################################################################
tp <- "12_weeks"
for(i in n.peer){
  cat.rds(paste0(inpath, tp, "/",  "peertest", i, "/chr1/"), "*.rds",
          paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
}

peer.res_12wk <- data.frame("NPeerFactors"=n.peer,
                            "NeGenes"=NA,
                            "NSigPairs"=NA)
sing_fit_df <- c()

for(i in n.peer){
  eqtl <- readRDS(paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
  print(length(unique(eqtl$Gene)))
  eqtl.o <- eqtl[order(eqtl$eQTL_pval), ]
  eqtl.o$padj <- p.adjust(eqtl.o$eQTL_pval, method="bonferroni")
  print(table(eqtl.o$padj < 0.05))
  p = which(i == n.peer)
  peer.res_12wk[p, 3] <- length(which(eqtl.o$padj < 0.05))
  eqtl.u <- subset(eqtl.o, !duplicated(Gene))
  print(table(eqtl.u$padj < 0.05))
  peer.res_12wk[p, 2] <- length(which(eqtl.u$padj < 0.05))
  prop_sing_fit_null <- eqtl %>% filter(singular_fit_null != 0) %>% nrow()/ eqtl %>% nrow()
  prop_sing_fit_test <- eqtl %>% filter(singular_fit_test != 0) %>% nrow()/ eqtl %>% nrow()
  sing_fit_df <- rbind(sing_fit_df, c(prop_sing_fit_null, prop_sing_fit_test))
}

peer.res_12wk$DiffeGenes <- c(peer.res_12wk$NeGenes[1], diff(peer.res_12wk$NeGenes, lag=1))
peer.res_12wk$DiffSigPairs <- c(peer.res_12wk$NSigPairs[1], diff(peer.res_12wk$NSigPairs, lag=1))
write.table(peer.res_12wk, paste0(inpath, "12_week_chr1_eqtl_number_by_peer.txt"), sep="\t")

p1_12wk <- ggplot(peer.res_12wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(aes(y=NeGenes + 10, label=NeGenes)) +
  geom_line() +
  theme_bw() +
  xlab("Number of PEER factors included") +
  ylab("Number of chr1 eGenes identified")

p2_12wk <- ggplot(peer.res_12wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(aes(y=DiffeGenes + 5, label=DiffeGenes)) +
  geom_line() +
  theme_bw() +
  ylim(c(0, 330)) +
  xlab("Number of PEER factors included") +
  ylab("Increase in number of chr1 eGenes identified")

# visualize
gridExtra::grid.arrange(p1_12wk, p2_12wk, nrow=1)

################################################################################
# 20 weeks only 
################################################################################
tp <- "20_weeks"
for(i in n.peer){
  print(i)
  cat.rds(paste0(inpath, tp, "/",  "peertest", i, "/chr1/"), "*.rds",
          paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
}

peer.res_20wk <- data.frame("NPeerFactors"=n.peer,
                            "NeGenes"=NA,
                            "NSigPairs"=NA)
sing_fit_df_20wk <- c()

for(i in n.peer){
  print(i)
  eqtl <- readRDS(paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
  print(length(unique(eqtl$Gene)))
  eqtl.o <- eqtl[order(eqtl$eQTL_pval), ]
  eqtl.o$padj <- p.adjust(eqtl.o$eQTL_pval, method="bonferroni")
  print(table(eqtl.o$padj < 0.05))
  p = which(i == n.peer)
  peer.res_20wk[p, 3] <- length(which(eqtl.o$padj < 0.05))
  eqtl.u <- subset(eqtl.o, !duplicated(Gene))
  print(table(eqtl.u$padj < 0.05))
  peer.res_20wk[p, 2] <- length(which(eqtl.u$padj < 0.05))
  prop_sing_fit_null <- eqtl %>% filter(singular_fit_null != 0) %>% nrow()/ eqtl %>% nrow()
  prop_sing_fit_test <- eqtl %>% filter(singular_fit_test != 0) %>% nrow()/ eqtl %>% nrow()
  sing_fit_df_20wk <- rbind(sing_fit_df_20wk, c(prop_sing_fit_null, prop_sing_fit_test))
}

peer.res_20wk$DiffeGenes <- c(peer.res_20wk$NeGenes[1], diff(peer.res_20wk$NeGenes, lag=1))
peer.res_20wk$DiffSigPairs <- c(peer.res_20wk$NSigPairs[1], diff(peer.res_20wk$NSigPairs, lag=1))
write.table(peer.res_20wk, paste0(inpath, "20_week_chr1_eqtl_number_by_peer.txt"), sep="\t")

p1_20wk <- ggplot(peer.res_20wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(aes(y=NeGenes + 10, label=NeGenes)) +
  geom_line() +
  theme_bw() +
  xlab("Number of PEER factors included") +
  ylab("Number of chr1 eGenes identified")

p2_20wk <- ggplot(peer.res_20wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(aes(y=DiffeGenes + 5, label=DiffeGenes)) +
  geom_line() +
  theme_bw() +
  ylim(c(0, 330)) +
  xlab("Number of PEER factors included") +
  ylab("Increase in number of chr1 eGenes identified")

# visualize
gridExtra::grid.arrange(p1_20wk, p2_20wk, nrow=1)

################################################################################
# 28 weeks only 
################################################################################
tp <- "28_weeks"
for(i in n.peer){
  print(i)
  cat.rds(paste0(inpath, tp, "/",  "peertest", i, "/chr1/"), "*.rds",
          paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
}

peer.res_28wk <- data.frame("NPeerFactors"=n.peer,
                            "NeGenes"=NA,
                            "NSigPairs"=NA)
sing_fit_df_28wk <- c()

for(i in n.peer){
  print(i)
  eqtl <- readRDS(paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
  print(length(unique(eqtl$Gene)))
  eqtl.o <- eqtl[order(eqtl$eQTL_pval), ]
  eqtl.o$padj <- p.adjust(eqtl.o$eQTL_pval, method="bonferroni")
  print(table(eqtl.o$padj < 0.05))
  p = which(i == n.peer)
  peer.res_28wk[p, 3] <- length(which(eqtl.o$padj < 0.05))
  eqtl.u <- subset(eqtl.o, !duplicated(Gene))
  print(table(eqtl.u$padj < 0.05))
  peer.res_28wk[p, 2] <- length(which(eqtl.u$padj < 0.05))
  prop_sing_fit_null <- eqtl %>% filter(singular_fit_null != 0) %>% nrow()/ eqtl %>% nrow()
  prop_sing_fit_test <- eqtl %>% filter(singular_fit_test != 0) %>% nrow()/ eqtl %>% nrow()
  sing_fit_df_28wk <- rbind(sing_fit_df_28wk, c(prop_sing_fit_null, prop_sing_fit_test))
}

peer.res_28wk$DiffeGenes <- c(peer.res_28wk$NeGenes[1], diff(peer.res_28wk$NeGenes, lag=1))
peer.res_28wk$DiffSigPairs <- c(peer.res_28wk$NSigPairs[1], diff(peer.res_28wk$NSigPairs, lag=1))
write.table(peer.res_28wk, paste0(inpath, "28_week_chr1_eqtl_number_by_peer.txt"), sep="\t")

p1_28wk <- ggplot(peer.res_28wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(aes(y=NeGenes + 10, label=NeGenes)) +
  geom_line() +
  theme_bw() +
  xlab("Number of PEER factors included") +
  ylab("Number of chr1 eGenes identified")

p2_28wk <- ggplot(peer.res_28wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(aes(y=DiffeGenes + 5, label=DiffeGenes)) +
  geom_line() +
  theme_bw() +
  ylim(c(0, 330)) +
  xlab("Number of PEER factors included") +
  ylab("Increase in number of chr1 eGenes identified")

# visualize
gridExtra::grid.arrange(p1_28wk, p2_28wk, nrow=1)

################################################################################
# 36 weeks only 
################################################################################
tp <- "36_weeks"
for(i in n.peer){
  print(i)
  cat.rds(paste0(inpath, tp, "/",  "peertest", i, "/chr1/"), "*.rds",
          paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
}

peer.res_36wk <- data.frame("NPeerFactors"=n.peer,
                            "NeGenes"=NA,
                            "NSigPairs"=NA)
sing_fit_df_36wk <- c()

for(i in n.peer){
  print(i)
  eqtl <- readRDS(paste0(inpath, tp, "_chr1_results_", i, "peer_factors.rds"))
  print(length(unique(eqtl$Gene)))
  eqtl.o <- eqtl[order(eqtl$eQTL_pval), ]
  eqtl.o$padj <- p.adjust(eqtl.o$eQTL_pval, method="bonferroni")
  print(table(eqtl.o$padj < 0.05))
  p = which(i == n.peer)
  peer.res_36wk[p, 3] <- length(which(eqtl.o$padj < 0.05))
  eqtl.u <- subset(eqtl.o, !duplicated(Gene))
  print(table(eqtl.u$padj < 0.05))
  peer.res_36wk[p, 2] <- length(which(eqtl.u$padj < 0.05))
  prop_sing_fit_null <- eqtl %>% filter(singular_fit_null != 0) %>% nrow()/ eqtl %>% nrow()
  prop_sing_fit_test <- eqtl %>% filter(singular_fit_test != 0) %>% nrow()/ eqtl %>% nrow()
  sing_fit_df_36wk <- rbind(sing_fit_df_36wk, c(prop_sing_fit_null, prop_sing_fit_test))
}

peer.res_36wk$DiffeGenes <- c(peer.res_36wk$NeGenes[1], diff(peer.res_36wk$NeGenes, lag=1))
peer.res_36wk$DiffSigPairs <- c(peer.res_36wk$NSigPairs[1], diff(peer.res_36wk$NSigPairs, lag=1))
write.table(peer.res_36wk, paste0(inpath, "36_week_chr1_eqtl_number_by_peer.txt"), sep="\t")

p1_36wk <- ggplot(peer.res_36wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(aes(y=NeGenes + 10, label=NeGenes)) +
  geom_line() +
  theme_bw() +
  xlab("Number of PEER factors included") +
  ylab("Number of chr1 eGenes identified")

p2_36wk <- ggplot(peer.res_36wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(aes(y=DiffeGenes + 5, label=DiffeGenes)) +
  geom_line() +
  theme_bw() +
  ylim(c(0, 330)) +
  xlab("Number of PEER factors included") +
  ylab("Increase in number of chr1 eGenes identified")

# visualize
gridExtra::grid.arrange(p1_36wk, p2_36wk, nrow=1)

########################### Plot ###########################
# peer factors single tp
peer.res_12wk <- read.table(paste0(inpath, "12_week_chr1_eqtl_number_by_peer.txt"))
peer.res_20wk <- read.table(paste0(inpath, "20_week_chr1_eqtl_number_by_peer.txt"))
peer.res_28wk <- read.table(paste0(inpath, "28_week_chr1_eqtl_number_by_peer.txt"))
peer.res_36wk <- read.table(paste0(inpath, "36_week_chr1_eqtl_number_by_peer.txt"))

p1_12wk <- ggplot(peer.res_12wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(
    aes(label = NeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_text(
    aes(label = NeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  xlab("N PEER factors included") +
  ylab("N chr1 eGenes identified") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_12wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("12 weeks")

p2_12wk <- ggplot(peer.res_12wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(
    aes(label = DiffeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  ylim(c(0, 330)) +
  xlab("N PEER factors included") +
  ylab("Increase in n eGenes") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_12wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("12 weeks")

p1_20wk <- ggplot(peer.res_20wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(
    aes(label = NeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  xlab("N PEER factors included") +
  ylab("N chr1 eGenes identified") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_20wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("20 weeks")

p2_20wk <- ggplot(peer.res_20wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(
    aes(label = DiffeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  ylim(c(0, 330)) +
  xlab("N PEER factors included") +
  ylab("Increase in n eGenes") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_20wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("20 weeks")

p1_28wk <- ggplot(peer.res_28wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(
    aes(label = NeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  xlab("N PEER factors included") +
  ylab("N chr1 eGenes identified") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_28wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("28 weeks")

p2_28wk <- ggplot(peer.res_28wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(
    aes(label = DiffeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  ylim(c(0, 330)) +
  xlab("N PEER factors included") +
  ylab("Increase in n eGenes") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_28wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("28 weeks")

p1_36wk <- ggplot(peer.res_36wk, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(
    aes(label = NeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  xlab("N PEER factors included") +
  ylab("N chr1 eGenes identified") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_36wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("36 weeks")

p2_36wk <- ggplot(peer.res_36wk, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(
    aes(label = DiffeGenes),
    vjust = -0.5,             # position labels slightly above the point
    size = 3                  # adjust size as needed
  ) + 
  geom_line() +
  ylim(c(0, 330)) +
  xlab("N PEER factors included") +
  ylab("Increase in n eGenes") + 
  scale_x_continuous(
    breaks = seq(0, max(peer.res_36wk$NPeerFactors), by = 5)
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggtitle("36 weeks")


gridExtra::grid.arrange(p1_12wk,
                                              p2_12wk,
                                              p1_20wk,
                                              p2_20wk,
                                              p1_28wk,
                                              p2_28wk,
                                              p1_36wk,
                                              p2_36wk,
                                              ncol = 2, nrow = 4, 
                                              widths = c(1, 1), 
                                              heights = c(1, 1, 1, 1),
                                              layout_matrix = rbind(c(1,2),
                                                                    c(3,4), 
                                                                    c(5,6),
                                                                    c(7,8)))

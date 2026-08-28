# 05_1_heteroskedasticity_btwn_tps.R 

################################################################################

# 5.1. Heteroskedasticity between time-points 

################################################################################

# Aim: Identify genes with heteroskedasticity (difference in variance) between time-points
# Note: script written by Celeste Cohe. Takes ~4hrs to run. 

########################### Output paths ##########################
outfile="rna-seq/analysis/variance/results/heteroskedasticity_betweentimepoints.tsv"
outfile_corrected="rna-seq/analysis/variance/results/adj_heteroskedasticity_betweentimepoints.tsv"

########################### Input paths ###########################
# sample info 
meta_fname = "rna-seq/data/sample_covariates_clin_tech.csv"
# count matrix 
cpm_fname="rna-seq/analysis/data_preprocessing/star-fc-genecounts_log2cpm_filt_samples_HLApm_2.txt"

########################### Parameters ############################
read.fast=function(fname,sep="\t",col.names=TRUE,row.names=TRUE){
  #message("Reading lines...")
  lines=readLines(fname)
  if(col.names){
    header=lines[1]
    firstrow=2
  } else {
    firstrow=1
  }
  #message("Splitting lines...")
  rows=strsplit(lines[firstrow:length(lines)],sep)
  #message("Making a dataframe...")
  df=as.data.frame(do.call(rbind,rows))
  if(row.names){
    rownames(df) <- df[[1]]
    df <- df[, -1]
  }
  if(col.names){
    colnames(df)=strsplit(header,sep)[[1]]}
  return(df)   
}

# Set cutoffs for plotting 
log_FC <- log2(1.5)
p_val <- 0.05

########################### Load packages ###########################
library(dplyr)
library(ggplot2)
library(stringr)
library(onewaytests)
library(car)

# set ggplot theme 
theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.ticks = element_blank(),
                             strip.background = element_blank()))

########################### Load data ###########################
message("Reading data...")
meta=read.csv(meta_fname)
cpm=read.fast(cpm_fname)

########################### Analysis ###########################



message("Preparing data...")
#Removing miscarriages
miscarriages=c("pops212970284", "pops212971619")
meta=meta[!meta$RNA_sanger_sample_id %in% miscarriages,]

#formatting 
meta$Sampletakenat=gsub(" ","",meta$Sampletakenat)
meta$Sampletakenat=paste0("timepoint_",meta$Sampletakenat)
timepoints=unique(meta$Sampletakenat)
tcpm=as.data.frame(t(cpm))
colnames(tcpm)=gsub("-","_",colnames(tcpm)) #so HLA-DRB stays complete in the bf.test formula
genes=unique(colnames(tcpm))
tcpm$RNA_sanger_sample_id=rownames(tcpm)

tp_pairs=combn(timepoints, 2, simplify = FALSE)

message("Selecting patients with all four timepoints...") #added after
counts=as.data.frame(table(meta$ANON_ID))
alltps_ids=counts$Var1[which(counts$Freq==4)]
meta=meta[which(meta$ANON_ID %in% alltps_ids),]

message("Testing data...")
res_genes=c()
res_timepoints=c()
res_pvals=c()
res_variance_ratios=c()
res_Fstats=c()

for(pair in tp_pairs){
  tp1=pair[1]
  tp2=pair[2]
  message(paste("Processing",tp1,"versus",tp2,"..."))
  for(gene in genes){
    meta_temp=meta %>%  select(ANON_ID,RNA_sanger_sample_id,id_run_position,Sampletakenat) %>% #this used to be subsetting to inds with paired samples in both timepoints but can be kept anyway
      filter(Sampletakenat %in% c(tp1,tp2)) %>%
      group_by(ANON_ID) %>%
      filter(n_distinct(Sampletakenat)==2) %>%
      ungroup() %>%
      left_join(tcpm[,c("RNA_sanger_sample_id",gene)],by="RNA_sanger_sample_id") %>% 
      na.omit()
    
    #batch correction
    meta_temp[[gene]]=as.numeric(meta_temp[[gene]])  
    meta_temp$resid=lm(paste(gene,"~ id_run_position"),data=meta_temp)$res
    
    #modelling
    lt=leveneTest(resid ~ as.factor(Sampletakenat),data=meta_temp)
    pval=lt[["Pr(>F)"]][1]
    Fstat=lt[["F value"]][1]
    
    var_earlytp=var(as.numeric(meta_temp$resid[meta_temp$Sampletakenat==tp1]))
    var_latetp=var(as.numeric(meta_temp$resid[meta_temp$Sampletakenat==tp2]))
    variance_ratio=var_latetp/var_earlytp
    
    res_genes=c(res_genes,gene)
    res_timepoints=c(res_timepoints,paste0(tp1,"_vs_",tp2)) #tp1 is the reference timepoint for the variance ratio
    res_variance_ratios=c(res_variance_ratios,variance_ratio)
    res_Fstats=c(res_Fstats,Fstat)
    res_pvals=c(res_pvals,pval)}}

res_genes=gsub("_","-",res_genes) #turn back into the original gene name
res_df=data.frame(gene=res_genes,
                  timepoints=res_timepoints,
                  variance_ratio=res_variance_ratios,
                  Fstat=res_Fstats,
                  pval=res_pvals)

message("Correcting p-values...")
res_df$qval=p.adjust(res_df$pval, method = "BH")
res_df$BYqval=p.adjust(res_df$pval, method = "BY")
res_df$log2FC=log2(res_df$variance_ratio)

message("Writing out data...")
write.table(res_df,outfile,quote=F,row.names=F,col.names=T,sep="\t")

# Perform multiple testing correction the same way it was done
# in the dream differntial expression analysis: 
#       Within each time-point comparison, perform BH correction 
# in topTable default adjust.method = "BH"
var_tests <- read.delim(outfile)

var_tests_update <- var_tests %>% group_by(timepoints) %>%
  mutate(BH_adj_pval_within_tp = p.adjust(pval, method = "BH")) %>%
  ungroup() 


# write out df with new p-values
var_tests_update %>% write.table(file = outfile_corrected, sep = "\t", row.names = FALSE)

message("Done!")

########################### Plot ###########################

heterosk_res <- read.delim(outfile_corrected) %>%
  mutate(comp = gsub("weeks", "", gsub("weeks_vs_", "v", gsub("timepoint_", "", timepoints))))

# Volcano plots
# make annotation to add n DEGs
n_degs <- heterosk_res %>% 
  dplyr::filter(abs(log2FC) > log_FC, 
                BH_adj_pval_within_tp < p_val) %>% 
  group_by(comp) %>%
  summarize(n = n())

heterosk_res %>% 
  left_join(n_degs, by = "comp") %>%
  mutate(reg = ifelse(((log2FC > log_FC) & (BH_adj_pval_within_tp < p_val)), "UP", 
                      ifelse(((log2FC < -log_FC) & (BH_adj_pval_within_tp < p_val)), "DOWN", "no_change")), 
         BH_adj_pval_within_tp = ifelse(BH_adj_pval_within_tp < .Machine$double.xmin, .Machine$double.xmin, BH_adj_pval_within_tp),
         comp = paste0(comp, ", n significant = ", n)
  ) %>%
  ggplot(aes(x=log2FC, y=-log10(BH_adj_pval_within_tp),# text = paste("Symbol:", gtf.gene_name), label = gtf.gene_name,
             color = reg)) +
  scale_color_manual(name = "reg",
                     values =c("no_change"= "#999999", 
                               "UP" = "black", 
                               "DOWN" = "black")) +
  geom_point(shape = 19, show.legend = FALSE, size = 0.5) +
  geom_hline(yintercept = -log10(p_val), linetype="longdash", colour="grey", size=0.5) +
  geom_vline(xintercept = log_FC, linetype="longdash", colour="#999999", size=0.5) +
  geom_vline(xintercept = -log_FC, linetype="longdash", colour="#999999", size=0.5) + 
  # labs(title="Volcano plots of differentially variable genes for all time-point contrasts") + #,
  #      #subtitle = paste("Fold-change > 1.5, adj.P.Val < 0.05")) + 
  theme(legend.position="none", 
        plot.title = element_text(size = 15),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)) + 
  coord_cartesian(clip = "off") + 
  xlab("Log 2 fold change") + 
  ylab("-log10(adjusted p-value)") + 
  facet_wrap(~comp, scales = "free", nrow = 2)

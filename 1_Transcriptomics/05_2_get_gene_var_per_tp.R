# 05_2_get_gene_var_per_tp.R 

################################################################################

# 5.2. Get gene variance per time-point 

################################################################################

# Aim: Find the variance of each gene at each time-point 
# Note: Code written by Celeste Cohen. Takes ~20 mins to run. 

########################### Output paths ##########################
outfile="rna-seq/analysis/variance/results/gene_variance_mean_pertimepoint.tsv"

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
    firstrow=1}
  #message("Splitting lines...")
  rows=strsplit(lines[firstrow:length(lines)],sep)
  #message("Making a dataframe...")
  df=as.data.frame(do.call(rbind,rows))
  if(row.names){
    rownames(df) <- df[[1]]
    df <- df[, -1]}
  if(col.names){
    colnames(df)=strsplit(header,sep)[[1]]}
  return(df)}

########################### Load packages ###########################
library(dplyr)

########################### Load data ###########################
message("Reading data...")
meta=read.csv(meta_fname)
cpm=read.fast(cpm_fname)

########################### Analysis ###########################

#Removing miscarriages
miscarriages=c("pops212970284", "pops212971619")
meta=meta[!meta$RNA_sanger_sample_id %in% miscarriages,]

#added below
tcpm=as.data.frame(t(cpm))
colnames(tcpm)=gsub("-","_",colnames(tcpm)) #so HLA-DRB stays complete in the bf.test formula
allgenes=unique(colnames(tcpm))
tcpm$RNA_sanger_sample_id=rownames(tcpm)

message("Selecting patients with all four timepoints...")
counts=as.data.frame(table(meta$ANON_ID))
alltps_ids=counts$Var1[which(counts$Freq==4)]
meta=meta[which(meta$ANON_ID %in% alltps_ids),]

message("Calculating variances and means...")
tps=c()
genes=c()
vars=c()
means=c()
uncor_vars=c()
uncor_means=c()
for(tp in unique(meta$Sampletakenat)){
  message(paste0("Calculating variances for ",tp,"..."))
  for(gene in allgenes){ #used to be rownames(cpm)
    meta_temp=meta %>%  select(ANON_ID,RNA_sanger_sample_id,id_run_position,Sampletakenat) %>% #this used to be subsetting to inds with paired samples in both timepoints but can be kept anyway
      filter(Sampletakenat==tp) %>%
      left_join(tcpm[,c("RNA_sanger_sample_id",gene)],by="RNA_sanger_sample_id") %>% 
      na.omit()
    
    uncor_genevar=var(as.numeric(meta_temp[[gene]]))
    uncor_genemean=mean(as.numeric(meta_temp[[gene]]))
    
    #batch correction
    meta_temp[[gene]]=as.numeric(meta_temp[[gene]])  
    meta_temp$resid=lm(paste(gene,"~ id_run_position"),data=meta_temp)$res
    
    genevar=var(as.numeric(meta_temp$resid))
    genemean=mean(as.numeric(meta_temp$resid))
    
    tps=c(tps,tp)
    genes=c(genes,gene)
    vars=c(vars,genevar)
    means=c(means,genemean)
    uncor_vars=c(uncor_vars,uncor_genevar)
    uncor_means=c(uncor_means,uncor_genemean)}}

message("Writing out data...")
genes=gsub("_","-",genes)
res=data.frame(timepoint=tps,gene=genes,uncorrected_variance=uncor_vars,variance=vars,uncorrected_mean=uncor_means,mean=means)
#saveRDS(res,outfile)
write.table(res,outfile,quote=F,row.names=F,col.names=T,sep="\t")

message("Done!")
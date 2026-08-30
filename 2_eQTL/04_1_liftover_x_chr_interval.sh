# 04_1_liftover_x_chr_interval.sh

################################################################################

# 4.1. Lift over interval X chr

################################################################################

# Aim: Lift over X chr in interval dataset (all other chrs had previously been lifted over)
# Note: following the steps post-imputation from 
# https://github.com/JonMarten/RNAseq/blob/master/covid-19/1_2_filter_chrx.sh

# start interactive job
fash 20
# load module with plink 
module load HGI/softpack/groups/trynka/popgen/1.0
# cd to correct dir
cd INTERVAL_RNAseq/genotype/chrX

# From the list of individuals in the other genotyping files
# make a list of samples to extract from the VCF files (with paired RNA-seq)
# IID IID
awk '{print $2, $2}' INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_AllAutosomes.fam > chrX/affy_ids_all_intervalRNA.txt

# turn the vcf into a plink file 
# 1. vcf input
# 2. restrict analysis to chr 23 (X)
# 3. keep only samples listed in file (2 columns FID and IID)
# 4. converts VCF into plink file 
# 5. set prefix for output files
plink\
 --vcf INTERVAL_X_imp_ann_filt_v2.vcf.gz \
 --chr 23 \
 --keep affy_ids_all_intervalRNA.txt \
 --make-bed \
 --out INTERVAL_chrX_merged_cleaned_RNAseq 
 
 # update positions and rename SNPs to lift over to b38
 # 1. plink input
 # 2. update variant positions using provided file to b38
 # 3. makes new plink dataset after updating the map 
 # 4. sets output file prefix 
 plink\
  --bfile INTERVAL_chrX_merged_cleaned_RNAseq \
  --update-map INTERVAL_chrX_update_b38pos.txt \
  --make-bed \
  --out INTERVAL_chrX_merged_cleaned_RNAseq_b38
  
# update positions to RSIDs
 # 1. plink input
 # 2. update variant names to rsid using file provided
 # 3. makes new plink dataset after updating the map 
 # 4. sets output file prefix 
 plink\
  --bfile INTERVAL_chrX_merged_cleaned_RNAseq_b38\
  --update-name INTERVAL_chrX_update_rsid.txt\
  --make-bed\
  --out INTERVAL_chrX_merged_cleaned_RNAseq_b38_rsids
 
 # # remove one duplicated SNP with different b37 names
 # made in https://github.com/JonMarten/RNAseq/blob/master/covid-19/1_2a_make_SNP_update_maps.R
  plink\
  --bfile INTERVAL_chrX_merged_cleaned_RNAseq_b38_rsids\
  --extract INTERVAL_chrX_b38snpfilter.txt\
  --make-bed\
  --out INTERVAL_chrX_merged_cleaned_RNAseq_b38_rsids

# filter for MAF 
plink \
  --bfile INTERVAL_chrX_merged_cleaned_RNAseq_b38_rsids \
  --maf 0.005 \
  --geno 0.05 \
  --make-bed \
  --out ../INTERVAL_RNAseq_imputed_b38_biallelic_MAF0.005_chrX

# run sex check 
plink \
--bfile INTERVAL_chrX_merged_cleaned_RNAseq_b38_rsids \
--check-sex \
--out INTERVAL_chrX_sexcheck

 
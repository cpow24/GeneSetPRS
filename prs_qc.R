#NOTE: code adapted from https://choishingwan.github.io/PRS-Tutorial/

#Setting working directory
setwd("ENTER DIRECTORY PATH")

#Loading required packages
library(data.table)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(scales)
library(magrittr)

#---------------------------------------------------------------------------------------
# GWAS Quality Control
#---------------------------------------------------------------------------------------

#Read GWAS file and rename columns
dat <- fread("Height_gwas_2018.txt.gz")
colnames(dat) <- c('CHR', 'BP', 'SNP', 'A1', 'A2', 'MAF', 'BETA', 'SE', 'P', 'N')

#Filter out SNPs with MAF < 0.05 (as N < 1000 in target data)
result <- dat[MAF > 0.05]

#Write the updated file as a .gz file
fwrite(result, "Height.gz", sep="\t")

#Run the following command line code to remove duplicate SNPs
#gunzip -c Height.gz |
#awk '{seen[$3]++; if(seen[$3]==1){ print}}' |
#gzip - > Height.nodup.gz

#Remove ambiguous SNPs in the command line
#gunzip -c Height.nodup.gz |
#awk '!( ($4=="A" && $5=="T") || 
#        ($4=="T" && $5=="A") || 
#        ($4=="G" && $5=="C") || 
#        ($4=="C" && $5=="G")) {print}' |
#gzip > Height.QC.gz

#---------------------------------------------------------------------------------------
# Target Data QC
#---------------------------------------------------------------------------------------

#Create initial PLINK files (BIM, BED, and FAM)
system(paste ("plink --vcf genotyping_data_fullset_train.vcf.gz --biallelic-only strict --make-bed --out OpenSNP_train",sep=" "))
system(paste ("plink --vcf genotyping_data_fullset_test.vcf.gz --biallelic-only strict --make-bed --out OpenSNP_test",sep=" "))

#Reading train and test bim files
train_bim <- fread("OpenSNP_train.bim")

#Creating column of Chromosome-rsid variable
train_bim <- train_bim %>%
  mutate(V7 = paste(V1, V2, sep="_"))

#Remove duplicated SNP IDs
train_bim2 <- subset(train_bim, !duplicated(V2))

#Remove '.' from VCF file
train_bim3 <- subset(train_bim2,V2!= ".")

write.table(x=train_bim3$V2,file=c("SNPS_keep.txt"),sep="",quote = F,col.names = F,row.names = F)
rm(train_bim, train_bim2, train_bim3)

#Create files with non-duplicated SNPs
system(paste ("plink --bfile OpenSNP_train --extract SNPS_keep.txt --make-bed --out target_train_nosex",sep=" "))
system(paste ("plink --bfile OpenSNP_test --extract SNPS_keep.txt --make-bed --out target_test_nosex",sep=" "))

#Impute sex
system(paste ("plink --bfile target_train_nosex --impute-sex ycount --make-bed --out  target_train_sex",sep=" "))
system(paste ("plink --bfile target_test_nosex --impute-sex ycount --make-bed --out  target_test_sex",sep=" "))

#Reading sexheck files
sexcheck_train <-fread("target_train_sex.sexcheck")
sexcheck_test <-fread("target_test_sex.sexcheck")

#Imputing sex based on F-statistic thresholds of 0.2 and 0.8
sexcheck_train <- sexcheck_train %>% 
  mutate(Sex= if_else(F <= 0.2, 2, if_else(F >= 0.8, 1, 0)))
sexcheck_test <- sexcheck_test %>% 
  mutate(Sex = if_else(F <= 0.2, 2, if_else(F >= 0.8, 1, 0)))

#Adding Sex to .FAM files
fam_train <- fread('target_train_sex.FAM')
fam_train$V5 <- sexcheck_train$Sex

fam_test <- fread('target_test_sex.FAM')
fam_test$V5 <- sexcheck_test$Sex

#Get summary of imputed sex values
table(fam_train$V5)
table(fam_test$V5)

#Get IID and FID of ambiguous sex individuals
no_sex_df <- fam_train[fam_train$V5 == 0]
no_sex_ids <- data.table(no_sex_df[, c(V1)])
no_sex_ids$V2 <- no_sex_df$V2
write.table(x=no_sex_ids,file=c("SNPS_nosex.txt"),sep=" ",quote = F,col.names = F,row.names = F)

write.table(fam_train, "target_train_sex.FAM", col.names=F, row.names=F, quote=F)
write.table(fam_test, "target_test_sex.FAM", col.names=F, row.names=F, quote=F)

#Remove individuals with ambiguous sex from PLINK files
system(paste ("plink --bfile target_train_sex --remove SNPS_nosex.txt --make-bed --out  target_train",sep=" "))
system(paste ("plink --bfile target_test_sex --make-bed --out  target_test",sep=" "))

#Ensuring train and test phenotype files are in the correct format
train_pheno <- fread('train_height.txt')
test_pheno <- fread('test_height.txt')

#Creating required datatable with FID and IID columns
cols <- c('FID', 'IID', 'height')
train_pheno$FID <- train_pheno$id
train_pheno$IID <- train_pheno$id
train_pheno <- subset(train_pheno, ,cols)

#Remove ambiguous sex ids
train_pheno <- subset(train_pheno, !(FID %in% c(no_sex_ids$V2)))

test_pheno$FID <- test_pheno$id
test_pheno$IID <- test_pheno$id
test_pheno <- subset(test_pheno, ,cols)

#Writing files
write.table(train_pheno, "height_train.HEIGHT", col.names=T, row.names=F, quote=F )
write.table(test_pheno, "height_test.HEIGHT", col.names=T, row.names=F, quote=F )

#Merge train and test data for QC
system(paste ("plink --bfile target_train --bmerge target_test --make-bed --out target_combined",sep=" "))

#Removing variants with multiple positions from both datasets
system(paste ("plink --bfile target_train --exclude target_combined-merge.missnp --make-bed --out target_train1",sep=" "))
system(paste ("plink --bfile target_test --exclude target_combined-merge.missnp --make-bed --out target_test1",sep=" "))

#Creating error-free merged file
system(paste ("plink --bfile target_train1 --bmerge target_test1 --make-bed --out target_combined",sep=" "))

#Remove SNPS with low genotyping rate, low MAF, and out of Hardy-Weinberg Equilibrium
system(paste ("plink --bfile target_combined --maf 0.05 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out target_combined.QC",sep=" "))

#Prune highly correlated SNPs
system(paste ("plink --bfile target_combined --keep target_combined.QC.fam --extract target_combined.QC.snplist --indep-pairwise 200 50 0.25 --out target_combined.QC",sep=" "))

#Compute heterozygosity rates
system(paste ("plink --bfile target_combined --extract target_combined.QC.prune.in --keep target_combined.QC.fam --het --out target_combined.QC",sep=" "))

#Remove individuals with F coefficients (from heterozygosity calcs) > 3 SD from mean
# Read in file
dat_target <- fread("target_combined.QC.het")

#Get samples with F coefficient within 3 SD of the population mean
valid_target <- dat_target[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 

#print FID and IID for valid samples
fwrite(valid_target[,c("FID","IID")], "target_combined.valid.sample", sep="\t") 

#Resolve mismatching SNPs

#Reading in .BIM file
bim_target <- fread("target_combined.bim") %>%
  setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
  .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]


#Reading in summary statistics
height <- fread("Height.QC.gz") %>%
  .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]

#Reading in qulaity controlled SNPs
qc <- fread("target_combined.QC.snplist", header=F)

#Merging summary statistic data and individual-level data
info <- merge(bim_target, height, by=c("SNP", "CHR", "BP")) %>%
  .[SNP %in% qc[,V1]]

#Function to identify complementary alleles
comp <- function(x){
  switch (x,"A" = "T",
          "C" = "G",
          "T" = "A",
          "G" = "C",
          return(NA)
  )
} 
#Identify SNPs with same alleles in base and target
info_match <- info[A1 == B.A1 & A2 == B.A2, SNP]

#Identify complementary SNPs
comp_snps <- info[sapply(B.A1, comp) == A1 &
                   sapply(B.A2, comp) == A2, SNP]

#Updating BIM file
bim_target[SNP %in% comp_snps, c("B.A1", "B.A2") :=
      list(sapply(B.A1, comp),
           sapply(B.A2, comp))]

#Identifying SNPs that need to be recoded
snp_rec <- info[B.A1==A2 & B.A2==A1, SNP]

#Update the bim file
bim_target[SNP %in% snp_rec, c("B.A1", "B.A2") :=
      list(B.A2, B.A1)]

#Identify SNPs that need to be recoded
recode_comp <- info[sapply(B.A1, comp) == A2 &
                     sapply(B.A2, comp) == A1, SNP]

#Updating BIM file
bim_target[SNP %in% recode_comp, c("B.A1", "B.A2") :=
      list(sapply(B.A2, comp),
           sapply(B.A1, comp))]

#Writing updated BIM file
fwrite(bim_target[,c("SNP", "B.A1")], "target_combined.a1", col.names=F, sep="\t")

#Identify SNPs that have different alleles in base and target data
mismatch <- bim_target[!(SNP %in% info_match |
                    SNP %in% comp_snps |
                    SNP %in% snp_rec |
                    SNP %in% recode_comp), SNP]
write.table(mismatch, "target_combined.mismatch", quote=F, row.names=F, col.names=F)


#Removing closely related individuals
system(paste ("plink --bfile target_combined --extract target_combined.QC.prune.in --keep target_combined.valid.sample --rel-cutoff 0.125 --out target_combined.QC",sep=" "))

#Run the following in the command line to remove duplicates
#cut -f 2 target_combined.bim | sort | uniq -d > target_combined.dups

system(paste ("plink --bfile target_combined  --exclude target_combined.dups --make-bed --out target_combined.QC_nodup",sep=" "))

#Generating quality controlled files
system(paste ("plink --bfile target_combined.QC_nodup --make-bed --keep target_combined.QC.rel.id --out target_combined.QC --extract target_combined.QC.snplist --exclude target_combined.mismatch --a1-allele target_combined.a1",sep=" "))

#Align GWAS and individual-level data
df <- fread("target_combined.QC.bim")
sumstats <- fread("Height.QC.gz")

sumstats <- sumstats %>% rename(rsid=SNP, chr=CHR, pos=BP, a1=A1, a0=A2, beta=BETA)
df <- df %>% rename(rsid=V2, chr=V1, unk = V3, pos=V4, a1=V5, a0=V6)

out <- snp_match(
  sumstats,
  df,
  strand_flip = TRUE,
  join_by_pos = FALSE,
  remove_dups = TRUE,
  match.min.prop = 0.2,
  return_flip_and_rev = TRUE
)

#Writing aligned GWAS data
write.table(out, "Height.gwas.aligned.txt", col.names=T, row.names=F, quote=F)

#Creating covariate files with sex included
sexcheck_train <-fread("target_train_sex.sexcheck")
sexcheck_test <-fread("target_test_sex.sexcheck")

sexcheck_train <- sexcheck_train %>% 
  mutate(Sex= if_else(F <= 0.2, 2, if_else(F >= 0.8, 1, 0)))
sexcheck_test <- sexcheck_test %>% 
  mutate(Sex = if_else(F <= 0.2, 2, if_else(F >= 0.8, 1, 0)))

write.table(sexcheck_train[,c('FID', 'IID', 'Sex')], 'target_train.COV',col.names=T, row.names=F, quote=F)
write.table(sexcheck_test[,c('FID', 'IID', 'Sex')], 'target_test.COV',col.names=T, row.names=F, quote=F)

#Getting list of all valid individual IDs after QC
fam <- fread('target_combined.QC.fam')
train_data <- fread('height_train.Height')
test_data <- fread('height_test.Height')

train_ids <- intersect(fam$V1, train_data$FID)
test_ids <- intersect(fam$V1, test_data$FID)

train_data <- subset(train_data, FID %in% train_ids)
test_data <- subset(test_data, FID %in% test_ids)

#Writing files
write.table(train_data, "height_train.HEIGHT", col.names=T, row.names=F, quote=F )
write.table(test_data, "height_test.HEIGHT", col.names=T, row.names=F, quote=F )

#Create separate training and test QC files
train_ids_df <- data.table(train_ids)
train_ids_df$V2 <- train_ids
colnames(train_ids_df) <- c('V1','V2')

write.table(x=train_ids_df,file=c("train_ids.txt"),sep=" ",quote = F,col.names = F,row.names = F)

test_ids_df <- data.table(test_ids)
test_ids_df$V2 <- test_ids
colnames(test_ids_df) <- c('V1','V2')

write.table(x=test_ids_df,file=c("test_ids.txt"),sep=" ",quote = F,col.names = F,row.names = F)

#Final quality controlled files
system(paste ("plink --bfile target_combined --keep train_ids.txt --make-bed --out  target_train.QC",sep=" "))
system(paste ("plink --bfile target_combined --keep test_ids.txt --make-bed --out  target_test.QC",sep=" "))

#Create files containing just Chromosomes 1-22
system(paste ("plink --bfile target_train.QC --chr 1-22 --make-bed --out target_train.QC_final",sep=" "))
system(paste ("plink --bfile target_test.QC --chr 1-22 --make-bed --out target_test.QC_final",sep=" "))

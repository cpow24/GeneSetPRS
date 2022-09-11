#Seting working directory
setwd('ENTER DIRECTORY PATH')

#Loading required packages
library(data.table)
library(bigsnpr)
library(dplyr)
library(scales)
library(magrittr)

#Reading in GWAS and saving as .Tranformed file
dat <- read.table('Height.gwas.aligned.txt', header=T)
dat <- dat %>% rename(SNP=rsid, CHR=chr, BP=pos, A1=a1, A2=a0, BETA=beta, POS=pos.ss)
dat <- dat[ c('CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'MAF')]
setcolorder(dat, c('CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'MAF'))
write.table(dat, "Height.QC.Transformed", col.names=T, quote=F, row.names=F)

#Perform Clumping (i.e retain weakly correlated SNPs with preference for keeping SNPs most associated with height)
system(paste ("plink --bfile target_train.QC --clump-p1 1 --clump-r2 0.1 --clump-kb 250
              --clump Height.QC.Transformed --clump-snp-field SNP --clump-field P --out target_train",sep=" "))

#Extract clumped SNP IDs using the following command line code:
#awk 'NR!=1{print $3}' target_train.clumped >  target_train.valid.snp

#Create file containing SNP IDs and their corresponding P-values
#awk '{print $3,$8}' Height.QC.Transformed > SNP.pvalue

#Create a file of P-value thresholds
#echo "0.00001 0 0.00001" > range_list
#echo "0.0001 0 0.0001" >> range_list
#echo "0.001 0 0.001" >> range_list 
#echo "0.05 0 0.05" >> range_list
#echo "0.1 0 0.1" >> range_list
#echo "0.2 0 0.2" >> range_list
#echo "0.3 0 0.3" >> range_list
#echo "0.4 0 0.4" >> range_list
#echo "0.5 0 0.5" >> range_list

#Generate polygenic risk scores for each of the P-value thresholds
#Training data
system(paste ("plink --bfile target_train.QC --score Height.QC.Transformed 3 4 9 header --q-score-range range_list SNP.pvalue --extract target_train.valid.snp --out target_train",sep=" "))

#Test data
system(paste ("plink --bfile target_test.QC --score Height.QC.Transformed 3 4 9 header --q-score-range range_list SNP.pvalue --extract target_train.valid.snp --out target_test",sep=" "))
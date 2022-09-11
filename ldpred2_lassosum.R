#Based on code from https://choishingwan.github.io/PRS-Tutorial/

#Seting working directory
setwd('ENTER DIRECTORY PATH')

#Loading required packages
library(data.table)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(scales)
library(magrittr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#Read in phenotype and covariate files
phenotype <- fread("target_train.height")
covariate <- fread("target_train.cov")

#Renaming columns
colnames(pcs) <- c("FID","IID", paste0("PC",1:10))

# generate required table
pheno <- merge(phenotype, covariate)

#Remove height_class and id variables
pheno <- pheno[,c('FID','IID','height','Sex')]

#Get HapMap3 SNPs
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"))

#Read summary statistic file
sumstats <- bigreadr::fread2('Height.QC.gz')

sumstats <- sumstats[ c('CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'MAF')]

#Renaming summary statistic columns for bignsnpr
names(sumstats) <-
  c("chr","pos","rsid","a1","a0","n_eff","beta_se","p","beta","MAF")

#Filtering for hapmap SNPs
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

#Max amount of cores
NCORES <- nb_cores()

#Creating temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

#Initialising variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL

#Order of samples in the bed file 
order_fam <- NULL

# preprocess the bed file (only need to do once for each data set)
snp_readBed("target_train.QC_final.bed")

# now attach the genotype object
obj.bigSNP <- snp_attach("target_train.QC_final.rds")

# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")

# perform SNP matching
info_snp <- snp_match(
  sumstats,
  map,
  strand_flip = TRUE,
  join_by_pos = FALSE,
  remove_dups = TRUE,
  match.min.prop = 0.2,
  return_flip_and_rev = TRUE
)

# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes

# Rename the data structures
CHR <- map$chr
POS <- map$pos

#Get genetic positions of SNPs using 1000 Genomes reference
genetic_pos <- snp_asGeneticPos(CHR, POS, dir = ".")

#Generating LD matrix
for (chr in 1:22) {
  #Extract chromosome SNPs
  ind_chr <- which(info_snp$chr == chr)
  ind_chr2 <- info_snp$`_NUM_ID_`[ind_chr]
  #Calculate LD score
  corr0 <- snp_cor(
    genotype,
    ind.col = ind_chr2,
    ncores = NCORES,
    infos.pos = genetic_pos[ind_chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

#Order of fam file
order_fam <- as.data.table(obj.bigSNP$fam)

#Rename fam order columns
setnames(order_fam,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

#Perform LD score regression
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)

h2_est <- ldsc[["h2"]]

write.table(ldsc, 'ldsc.txt', col.names=F, quote=F, row.names=F)
write.table(ld, 'ld_scores.txt', col.names=F, quote=F, row.names=F)

#Ensure ordering of phenotype and genotype files are the same
y <- pheno[order_fam, on = c("FID", "IID")]

#Prepare data for grid model - Round 1
p_seq <- c(signif(seq_log(1e-4, 1, length.out = 8), 2))
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid_params <-
  expand.grid(p = p_seq,
              h2 = h2_seq,
              sparse = c(FALSE, TRUE))

#Get updated beta from grid model - Round 2
beta_grid <-
  snp_ldpred2_grid(corr, df_beta, grid_params, ncores = NCORES)
write.table(beta_grid, 'beta_grid.txt', col.names=F, quote=F, row.names=F)

#Prepare data for grid model - Round 2
p_seq_2 <- c(signif(seq_log(0.3, 1, length.out = 8)))
h2_seq_2 <- round(h2_est * c(0.8, 0.9, 1.1), 4)
grid_params_2 <-
  expand.grid(p = p_seq_2,
              h2 = h2_seq_2,
              sparse = c(FALSE, TRUE))

#Get updated beta from grid model - Round 2
beta_grid_2 <-
  snp_ldpred2_grid(corr, df_beta, grid_params_2, ncores = NCORES)
write.table(beta_grid, 'beta_grid_2.txt', col.names=F, quote=F, row.names=F)

#Generate beta values for LDpred2-inf model
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
write.table(beta_inf, 'beta_inf.txt', col.names=F, quote=F, row.names=F)

#Get updated beta from the auto model
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
  ncores = NCORES
)

beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)

write.table(beta_auto, 'beta_auto.txt', col.names=F, quote=F, row.names=F)

#Get updated beta from lassosum
beta_lassosum <- snp_lassosum2(
  corr,
  df_beta,
  delta = signif(seq_log(0.001, 3, 6), 1),
  nlambda = 20,
  lambda.min.ratio = 0.01,
  dfmax = 2e+05,
  maxiter = 500,
  tol = 1e-05,
  ncores = NCORES
)

lasso_params <- attr(beta_lassosum, "grid_param")

write.csv(lasso_params, 'lassosum_params.csv')
write.table(beta_lassosum, 'beta_lassosum.txt', col.names=F, quote=F, row.names=F)

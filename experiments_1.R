#Seting working directory
setwd('ENTER DIRECTORY PATH')

#Loading required packages
library(data.table)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(scales)
library(magrittr)
library(MLmetrics)
set.seed(1234)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#Getting list of valid sample IDs after QC
train_ids <- fread("train_ids.txt")
test_ids <- fread("test_ids.txt")

#Read in phenotype and covariate files
phenotype <- fread("target_train.height")
covariate <- fread("target_train.cov")
test_phenotype <- fread("target_test.height")
test_covariate <- fread("target_test.cov")

#Filtering for required columns only
phenotype <- phenotype[,c('FID', 'IID', 'height')]
test_phenotype <- test_phenotype[,c('FID', 'IID', 'height')]

#Filtering for valid IDs
phenotype <- subset(phenotype, IID %in% train_ids$V1)
covariate <- subset(covariate, IID %in% train_ids$V1)
test_phenotype <- subset(test_phenotype, IID %in% test_ids$V1)
test_covariate <- subset(test_covariate, IID %in% test_ids$V1)

#Merging phenotype with covariates
pheno <- merge(phenotype, covariate)
test_pheno <- merge(test_phenotype, test_covariate)

#Get HapMap3 SNPs
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"))

#Read summary statistic file
sumstats <- bigreadr::fread2('Height.QC.gz')

sumstats <- sumstats[ c('CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'MAF')]

names(sumstats) <-
  c("chr","pos","rsid","a1","a0","n_eff","beta_se","p","beta","MAF")

#Filter for hapmap SNPs only
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

#Ordering of samples in the bed file 
fam_order <- NULL
test_fam_order <- NULL

#Preprocess bed files
snp_readBed("target_train.QC_final.bed")
snp_readBed("target_test.QC_final.bed")

#Attach genotype object
obj_bigSNP <- snp_attach("target_train.QC_final.rds")
test_obj_bigSNP <- snp_attach("target_test.QC_final.rds")

#Extract the SNP information from the genotype
map <- obj_bigSNP$map[-3]
test_map <- test_obj_bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
names(test_map) <- c("chr", "rsid", "pos", "a1", "a0")

#Rename chromosome and genetic position data structures
CHR <- map$chr
POS <- map$pos

#Match SNPs between target and GWAS data
info_snp <- snp_match(
  sumstats,
  map,
  strand_flip = TRUE,
  join_by_pos = FALSE,
  remove_dups = TRUE,
  match.min.prop = 0.2,
  return_flip_and_rev = TRUE
)

#Create train and test genotype variables
genotype <- obj_bigSNP$genotypes
test_genotype <- test_obj_bigSNP$genotypes

#Generating correct fam order
fam_order <- as.data.table(obj_bigSNP$fam)
test_fam_order <- as.data.table(test_obj_bigSNP$fam)

#Rename fam order columns
setnames(fam_order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

setnames(test_fam_order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

#Number of cores available
NCORES <- nb_cores()

#Reading in pre-saved values for h2 and LD
ld <- fread('ld_scores.txt')
ld <- ld$V1
ldsc <- fread('ldsc.txt')
h2_est <- as.double(ldsc[2])

#Ensuring the genotype and phenotype files have the same ordering
y <- pheno[fam_order, on = c("FID", "IID")]
test_y <- test_pheno[test_fam_order, on = c("FID", "IID")]

#LDPRED-2 Grid - Round 1 Parameters
p_seq <- signif(seq_log(1e-4, 1, length.out = 8), 2)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <-
  expand.grid(p = p_seq,
              h2 = h2_seq,
              sparse = c(FALSE, TRUE))

param_search_vals <- grid.param %>%
  mutate(col_title = paste0(as.character(p),'_',as.character(h2),'_',sparse))
grid_titles <- param_search_vals$col_title

#LDPRED-2 Grid - Round 2 Parameters
p_seq_2 <- c(signif(seq_log(0.3, 1, length.out = 8)))
h2_seq_2 <- round(h2_est * c(0.8, 0.9, 1.1), 4)
grid.param_2 <-
  expand.grid(p = p_seq_2,
              h2 = h2_seq_2,
              sparse = c(FALSE, TRUE))

param_search_vals_2 <- grid.param_2 %>%
  mutate(col_title = paste0(as.character(p),'_',as.character(h2),'_',sparse))
grid_titles_2 <- param_search_vals_2$col_title

#Reading in LDPRED-2 grid beta values - Gibbs Sampling Round 1
beta_grid <- fread('beta_grid.txt')
colnames(beta_grid) <- grid_titles

#Reading in LDPRED-2 grid beta values - Gibbs Sampling Round 2
beta_grid_2 <- fread('beta_grid_2.txt')
colnames(beta_grid_2) <- grid_titles_2

#Remove NA columns (i.e. where Gibbs sampler failed to converge)
beta_grid <- beta_grid %>%
  select_if(~ !any(is.na(.)))

beta_grid_2 <- beta_grid_2 %>%
  select_if(~ !any(is.na(.)))

#Reading in LDPred-2-inf beta values
beta_inf <- fread('beta_inf.txt')

#LDPred-2 Auto Beta values
beta_auto <- fread('beta_auto.txt')

pred_auto <-
  big_prodMat(genotype,
              as.matrix(beta_auto),
              ind.col = info_snp$`_NUM_ID_`)

#Scale the PRS generated from auto model
final_beta_auto <-
  as.data.table(rowMeans(beta_auto))

#Read in Lassosum beta values
beta_lasso <- fread('beta_lassosum.txt')
lasso_params <- read.csv('lassosum_params.csv')

#Name lassosum values based on values of lambda and s
lasso_titles <- list()
for (i in 1:120){
  col_str <- paste0('lasso_',round(lasso_params$lambda[i],5),'_',lasso_params$delta[i])
  lasso_titles[[i]] <- col_str
}

colnames(beta_lasso) <- c(unlist(lasso_titles, recursive=F))

#Merging grid, auto and inf beta values
beta_all <- beta_grid
beta_all <- cbind(beta_all, beta_grid_2)
beta_all$inf <- beta_inf$V1
beta_all$auto <- final_beta_auto$V1

#Subsetting valid Lassosum results
beta_lasso_valid <- beta_lasso[,c(17,18,19,20,37,38,39,40,57,58,59,60,77,78,79,80,97,98,99)]
beta_all <- cbind(beta_all, beta_lasso_valid)

#Creating IDs of folds for 10-fold CV
k <- 5
ids <- c(1:nrow(pheno))
ids_rand <- sample(ids)
set_size <- round(nrow(pheno)/k)
folds <- split(ids_rand, ceiling(seq_along(ids_rand)/set_size))

#Exporting 5 data folds for use in Python scripts
for (i in 1:5){
  write.table(x=folds[i],file=c(paste0("fold_",i,".txt")),sep=" ",quote = F,col.names = F,row.names = F)
}

#Generating principal Components for training-validation pairs in cross-validation
for (i in 1:5){
  
  #Train and validation IDs
  val_fold <- unlist(folds[i], recursive=FALSE)
  train_fold <- unlist(folds[-i], recursive=FALSE)
  
  pca_res <- bed_projectPCA(
    bed('target_train.QC_final.bed'),
    bed('target_train.QC_final.bed'),
    k = 20,
    ind.row.new = val_fold,
    ind.row.ref = train_fold,
    ind.col.ref = info_snp$'_NUM_ID_',
    strand_flip = TRUE,
    join_by_pos = FALSE,
    match.min.prop = 0.2,
    verbose = TRUE,
    ncores = NCORES
  )
  
  assign(paste0("pca_res_", i),pca_res)
  
}

#Scree plots of 5 training folds
all_pca_res <- data.frame(pca_res_1$obj.svd.ref$d)
colnames(all_pca_res) <- c('fold_1')
all_pca_res$fold_2 <- pca_res_2$obj.svd.ref$d
all_pca_res$fold_3 <- pca_res_3$obj.svd.ref$d
all_pca_res$fold_4 <- pca_res_4$obj.svd.ref$d
all_pca_res$fold_5 <- pca_res_5$obj.svd.ref$d
all_pca_res$pc_id <- c(1:20)

ggplot(all_pca_res, aes(x=pc_id))+
  geom_line(aes(y=fold_1, color='fold_1'), size=1)+
  geom_line(aes(y=fold_2, color='fold_2'), size=1)+
  geom_line(aes(y=fold_3, color='fold_3'), size=1)+
  geom_line(aes(y=fold_4, color='fold_4'), size=1)+
  geom_line(aes(y=fold_5, color='fold_5'), size=1)+
  labs(title="Scree Plot",
       x ="PC Index", y = "Singular Value", fill="Fold", colour='Fold') +
  theme(
    plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5, vjust=2),
    axis.title.x = element_text(color="Black", size=14, face="bold"),
    axis.title.y = element_text(color="Black", size=14, angle=90, face="bold"),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"),  complete = TRUE,
    legend.position = c(0.9, 0.85),
    legend.title = element_text(color = "black", size = 14, face="bold"),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = "grey92"),
    legend.key = element_rect(fill = "grey92", color = NA)
  )


#Projecting training PCs onto test data
test_pca <- bed_projectPCA(
  bed('target_train.QC_final.bed'),
  bed('target_test.QC_final.bed'),
  k = 20,
  ind.col.ref = info_snp$'_NUM_ID_',
  strand_flip = TRUE,
  join_by_pos = FALSE,
  match.min.prop = 0.2,
  verbose = TRUE,
  ncores = NCORES
)

#Generating PRS predictions using beta weights
pred_grid <- big_prodMat(genotype,
                         as.matrix(beta_all),
                         ind.col = info_snp$'_NUM_ID_')

#Appending clumping and thresholding PRS predictions
threshold_list <- c('0.00001','0.0001','0.001','0.05','0.1','0.2','0.3','0.4','0.5')
ct_results <- data.frame(matrix(ncol = 0, nrow = 738))
ct_results_test <- data.frame(matrix(ncol = 0, nrow = 129))

#Iterating through C+T score files and adding them to prediction grid - Training Set
for (thresh in threshold_list){
  score_file <- fread(paste0('target_train.',thresh,'.profile'))
  scores <- score_file$SCORE
  ct_results <- cbind(ct_results,scores)
  
}

colnames(ct_results) <- threshold_list

pred_grid <- cbind(pred_grid, ct_results)

#Iterating through C+T score files and adding them to prediction grid - Test Set
for (thresh in threshold_list){
  score_file <- fread(paste0('target_test.',thresh,'.profile'))
  scores <- score_file$SCORE
  ct_results_test <- cbind(ct_results_test,scores)
  
}

colnames(ct_results_test) <- threshold_list

test_pred_grid <- big_prodMat(test_genotype,
                              as.matrix(beta_all),
                              ind.col = info_snp$'_NUM_ID_')

test_pred_grid <- cbind(test_pred_grid, ct_results_test)

#Naming columns by model name
colnames(pred_grid) <- c('LDPRED2_1_0.2858_FALSE', 'LDPRED2_1_0.4083_FALSE','LDPRED2_1_0.2858_TRUE','LDPRED2_1_0.4083_TRUE','LDPRED2_0.708934_0.3267_FALSE',
                         'LDPRED2_0.841982_0.3267_FALSE', 'LDPRED2_1_0.3267_FALSE', 'LDPRED2_0.841982_0.3675_FALSE', 'LDPRED2_1_0.3675_FALSE', 'LDPRED2_1_0.4492_FALSE', 
                         'LDPRED2_0.708934_0.3267_TRUE', 'LDPRED2_0.841982_0.3267_TRUE', 'LDPRED2_1_0.3267_TRUE', 'LDPRED2_0.841982_0.3675_TRUE',
                         'LDPRED2_1_0.3675_TRUE','LDPRED2_1_0.4492_TRUE', 'LDPRED2_inf','LDPRED2_auto',
                         'ct_0.00001','ct_0.0001','ct_0.001','ct_0.05','ct_0.1','ct_0.2','ct_0.3','ct_0.4','ct_0.5', 'lasso_0.0013_0.001', 'lasso_0.00104_0.001', 
                         'lasso_0.00082_0.001', 'lasso_0.00065_0.001', 'lasso_0.0013_0.005', 'lasso_0.00104_0.005', 'lasso_0.00082_0.005', 'lasso_0.00065_0.005',
                         'lasso_0.0013_0.02', 'lasso_0.00104_0.02', 'lasso_0.00082_0.02', 'lasso_0.00065_0.02', 'lasso_0.0013_0.1', 'lasso_0.00104_0.1',
                         'lasso_0.00082_0.1', 'lasso_0.00065_0.1', 'lasso_0.0013_0.6', 'lasso_0.00104_0.6', 'lasso_0.00082_0.6')

colnames(test_pred_grid) <- c('LDPRED2_1_0.2858_FALSE', 'LDPRED2_1_0.4083_FALSE','LDPRED2_1_0.2858_TRUE','LDPRED2_1_0.4083_TRUE','LDPRED2_0.708934_0.3267_FALSE',
                              'LDPRED2_0.841982_0.3267_FALSE', 'LDPRED2_1_0.3267_FALSE', 'LDPRED2_0.841982_0.3675_FALSE', 'LDPRED2_1_0.3675_FALSE', 'LDPRED2_1_0.4492_FALSE', 
                              'LDPRED2_0.708934_0.3267_TRUE', 'LDPRED2_0.841982_0.3267_TRUE', 'LDPRED2_1_0.3267_TRUE', 'LDPRED2_0.841982_0.3675_TRUE',
                              'LDPRED2_1_0.3675_TRUE','LDPRED2_1_0.4492_TRUE', 'LDPRED2_inf','LDPRED2_auto',
                              'ct_0.00001','ct_0.0001','ct_0.001','ct_0.05','ct_0.1','ct_0.2','ct_0.3','ct_0.4','ct_0.5', 'lasso_0.0013_0.001', 'lasso_0.00104_0.001', 
                              'lasso_0.00082_0.001', 'lasso_0.00065_0.001', 'lasso_0.0013_0.005', 'lasso_0.00104_0.005', 'lasso_0.00082_0.005', 'lasso_0.00065_0.005',
                              'lasso_0.0013_0.02', 'lasso_0.00104_0.02', 'lasso_0.00082_0.02', 'lasso_0.00065_0.02', 'lasso_0.0013_0.1', 'lasso_0.00104_0.1',
                              'lasso_0.00082_0.1', 'lasso_0.00065_0.1', 'lasso_0.0013_0.6', 'lasso_0.00104_0.6', 'lasso_0.00082_0.6')

#Merging PRS scores with data
y <- cbind(y, pred_grid)
test_y <- cbind(test_y, test_pred_grid)

#Creating full data matrix for each cross-validation fold
#Validation values are given their projection onto the PCs computed with each training fold

create_folds <- function(fold, pca_res, data, fold_ids){
  
  #Getting validation and training IDs
  val_ids <- unlist(fold_ids[fold], recursive=FALSE)
  train_ids <- unlist(fold_ids[-fold], recursive=FALSE)
  
  #Creating dataframe
  y_new <- data
  
  #Assigning appropriate PCA projections
  #PC1
  y_new$PC1 <- 1:nrow(y)
  y_new$PC1[val_ids] <- pca_res$OADP_proj[,1]
  y_new$PC1[train_ids] <- predict(pca_res$obj.svd.ref)[,1]
  
  #PC2
  y_new$PC2 <- 1:nrow(y)
  y_new$PC2[val_ids] <- pca_res$OADP_proj[,2]
  y_new$PC2[train_ids] <- predict(pca_res$obj.svd.ref)[,2]
  
  #PC3
  y_new$PC3 <- 1:nrow(y)
  y_new$PC3[val_ids] <- pca_res$OADP_proj[,3]
  y_new$PC3[train_ids] <- predict(pca_res$obj.svd.ref)[,3]
  
  #PC4
  y_new$PC4 <- 1:nrow(y)
  y_new$PC4[val_ids] <- pca_res$OADP_proj[,4]
  y_new$PC4[train_ids] <- predict(pca_res$obj.svd.ref)[,4]
  
  #PC5
  y_new$PC5 <- 1:nrow(y)
  y_new$PC5[val_ids] <- pca_res$OADP_proj[,5]
  y_new$PC5[train_ids] <- predict(pca_res$obj.svd.ref)[,5]
  
  return(y_new)
  
}

#Creating and exporting full datasets for each round of CV
y1 <- create_folds(1, pca_res_1, y, folds)
y2 <- create_folds(2, pca_res_2, y, folds)
y3 <- create_folds(3, pca_res_3, y, folds)
y4 <- create_folds(4, pca_res_4, y, folds)
y5 <- create_folds(5, pca_res_5, y, folds)


write.csv(y1, 'folds_1.csv', row.names=FALSE)
write.csv(y2, 'folds_2.csv', row.names=FALSE)
write.csv(y3, 'folds_3.csv', row.names=FALSE)
write.csv(y4, 'folds_4.csv', row.names=FALSE)
write.csv(y5, 'folds_5.csv', row.names=FALSE)


#Test Set Principal Components
#PC1
test_y$PC1 <- 1:nrow(test_y)
test_y$PC1 <- test_pca$OADP_proj[,1]

#PC2
test_y$PC2 <- 1:nrow(test_y)
test_y$PC2 <-test_pca$OADP_proj[,2]

#PC3
test_y$PC3 <- 1:nrow(test_y)
test_y$PC3 <- test_pca$OADP_proj[,3]

#PC4
test_y$PC4 <- 1:nrow(test_y)
test_y$PC4 <- test_pca$OADP_proj[,4]

#PC5
test_y$PC5 <- 1:nrow(test_y)
test_y$PC5 <- test_pca$OADP_proj[,5]

#Dataframe of 5-fold CV results
cv_res <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(cv_res) <- c('fold', 'model', 'train_r2','null_r2', 'val_rmse', 'train_rmse', 'null_train_rmse', 'null_val_mse', 'sex_val_mse', 'val_mse', 'null_val_r2','sex_val_r2', 'val_r2')

fold_list <- list(y1,y2,y3,y4,y5)

#5-fold cross validation for all genome-wide models
for (i in 1:k){
  val_fold <- unlist(folds[i], recursive=FALSE)
  train_fold <- unlist(folds[-i], recursive=FALSE)
  y_curr <- fold_list[[i]]
  
  null_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5")), data = y_curr[train_fold,])
  null_summary <- summary(null_model)
  
  null_r2 <- null_summary$r.squared
  null_train_rmse <- sqrt(mean(null_summary$residuals^2))
  null_preds <- predict(null_model, y_curr[val_fold,]) 
  null_val_mse <- mean((y_curr[val_fold,]$height - null_preds)^2)
  null_val_r2 <- R2_Score(null_preds, y_curr[val_fold,]$height)
  
  sex_model <- lm(as.formula(paste0("height~Sex")), data = y_curr[train_fold,])
  sex_summary <- summary(null_model)
  
  sex_r2 <- sex_summary$r.squared
  sex_train_rmse <- sqrt(mean(sex_summary$residuals^2))
  sex_preds <- predict(sex_model, y_curr[val_fold,]) 
  sex_val_mse <- mean((y_curr[val_fold,]$height - sex_preds)^2)
  sex_val_r2 <- R2_Score(sex_preds, y_curr[val_fold,]$height)
  
  #Fit each model and get test + train R2
  for (model in colnames(pred_grid)){
    prs_model <- lm(as.formula(paste0('height ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + ',model)), data=y_curr[train_fold,]) 
    prs_summary <- summary(prs_model)
    train_r2 <- prs_summary$r.squared
    train_rmse <- sqrt(mean(prs_summary$residuals^2))
    
    preds <- predict(prs_model, y_curr[val_fold,])
    
    val_rmse <- sqrt(mean((y_curr[val_fold,]$height - preds)^2))
    val_mse <- mean((y_curr[val_fold,]$height - preds)^2)
    val_r2 <- R2_Score(preds, y_curr[val_fold,]$height)
    
    cv_res <- rbindlist(list(cv_res,list(i, model, train_r2, null_r2, val_rmse, train_rmse, null_train_rmse, null_val_mse, sex_val_mse, val_mse, null_val_r2, sex_val_r2, val_r2)))
  }
  
}

#Dataframe to store summary of CV results
model_summary <- data.frame(matrix(ncol = 23, nrow = 0))
colnames(model_summary) <- c('model', 'mean_train_r2', 'sd_train_r2', 'mean_null_train_r2', 'sd_null_train_r2', 'mean_null_val_r2','sd_null_val_r2','mean_val_rmse', 'sd_val_rmse', 
                             'mean_val_mse', 'sd_val_mse', 'mean_val_r2', 'sd_val_r2', 'mean_train_rmse', 'sd_train_rmse', 'mean_null_train_rmse', 'sd_null_train_rmse',
                             'mean_null_val_mse', 'sd_null_val_mse', 'mean_sex_val_mse', 'sd_sex_val_mse', 'mean_sex_val_r2', 'sd_sex_val_r2')

#Calculating model performance statistics
for (model_name in colnames(pred_grid)){
  temp_data <- subset(cv_res, model == model_name)
  
  #Train R2
  mean_train_r2 <- mean(temp_data$train_r2)
  sd_train_r2 <- sd(temp_data$train_r2)
  
  #Null Train R2
  mean_null_train_r2 <- mean(temp_data$null_r2)
  sd_null_train_r2 <- sd(temp_data$null_r2)
  
  #Null validation R2
  mean_null_val_r2 <- mean(temp_data$null_val_r2)
  sd_null_val_r2 <- sd(temp_data$null_val_r2)
  
  #Validation RMSE
  mean_val_rmse <- mean(temp_data$val_rmse)
  sd_val_rmse <- sd(temp_data$val_rmse)
  
  #Validation MSE
  mean_val_mse <- mean(temp_data$val_mse)
  sd_val_mse <- sd(temp_data$val_mse)
  
  #Validation R2
  mean_val_r2 <- mean(temp_data$val_r2)
  sd_val_r2 <- sd(temp_data$val_r2)
  
  #Train RMSE
  mean_train_rmse <- mean(temp_data$train_rmse)
  sd_train_rmse <- sd(temp_data$train_rmse)
  
  #Null model validation MSE
  mean_null_val_mse <- mean(temp_data$null_val_mse)
  sd_null_val_mse <- sd(temp_data$null_val_mse)
  
  #Null model train RMSE
  mean_null_train_rmse <- mean(temp_data$null_train_rmse)
  sd_null_train_rmse <- sd(temp_data$null_train_rmse)
  
  #Sex-only model validation RMSE
  mean_sex_val_mse <- mean(temp_data$sex_val_mse)
  sd_sex_val_mse <- sd(temp_data$sex_val_mse)
  
  #Sex-only model validation R2
  mean_sex_val_r2 <- mean(temp_data$sex_val_r2)
  sd_sex_val_r2 <- sd(temp_data$sex_val_r2)
  
  model_summary <- rbindlist(list(model_summary,list(model_name, mean_train_r2, sd_train_r2, mean_null_train_r2, sd_null_train_r2,
                                                     mean_null_val_r2, sd_null_val_r2, mean_val_rmse, sd_val_rmse, mean_val_mse, sd_val_mse,
                                                     mean_val_r2, sd_val_r2,
                                                     mean_train_rmse, sd_train_rmse, mean_null_train_rmse, sd_null_train_rmse,
                                                     mean_null_val_mse, sd_null_val_mse, mean_sex_val_mse, sd_sex_val_mse, mean_sex_val_r2, sd_sex_val_r2)))
  
}

write.csv(model_summary, 'model_summary.csv', row.names=FALSE)

#Best model based on Validation RMSE is LDpred-2 Grid model -> LDPRED2_0.841982_0.3267_TRUE

#Subsetting relevant test set columns inclding best PRs model scores
test_y <- test_y[,c('height','Sex','LDPRED2_0.841982_0.3267_TRUE','PC1','PC2','PC3','PC4','PC5', 'ct_0.5', 'LDPRED2_auto', 'LDPRED2_inf', 'lasso_0.00065_0.005')]
write.csv(test_y, 'final_test_data.csv', row.names=FALSE)

#--------------------------------------------------------------------------------
#GENE SET ANALYSIS - Training Set
#--------------------------------------------------------------------------------
best_preds <- pred_grid$LDPRED2_0.841982_0.3267_TRUE
best_betas <- beta_all$`0.841982_0.3267_TRUE`
gene_set_scores <- data.frame(matrix(ncol = 0, nrow = 738))

#Creating gene set scores for training data
for (i in 1:64){
  #Reading in set of rsids
  snp_set <- fread(paste0('snplist_set_',i,'.set'))
  colnames(snp_set) <- c('rsid')
  
  #Row indices of subset SNPs
  subset_ids <- which(info_snp$rsid %in% snp_set$rsid)
  info_snp_subset <- subset(info_snp, rsid %in% snp_set$rsid)
  
  #Get scores
  set_scores <- big_prodVec(genotype,
                            best_betas[subset_ids],
                            ind.col = info_snp_subset$`_NUM_ID_`)
  
  gene_set_scores <- cbind(gene_set_scores, set_scores)
}

colnames(gene_set_scores) <- paste0('set_', 1:64)
gene_set_scores$full_prs <- best_preds
gene_set_scores$sex <- y$Sex
gene_set_scores$height <- y$height
write.csv(gene_set_scores, 'set_scores.csv', row.names=FALSE)

#--------------------------------------------------------------------------------
#GENE SET ANALYSIS - Test Set
#--------------------------------------------------------------------------------
best_preds_test <- test_pred_grid$LDPRED2_0.841982_0.3267_TRUE
best_betas_test <- beta_all$`0.841982_0.3267_TRUE`
test_gene_set_scores <- data.frame(matrix(ncol = 0, nrow = 129))

#Creating gene set scores for test data
for (i in 1:64){
  #Reading in set of rsids
  snp_set <- fread(paste0('snplist_set_',i,'.set'))
  colnames(snp_set) <- c('rsid')
  
  #Row indices of subset SNPs
  subset_ids <- which(info_snp$rsid %in% snp_set$rsid)
  info_snp_subset <- subset(info_snp, rsid %in% snp_set$rsid)
  
  #Get scores
  test_set_scores <- big_prodVec(test_genotype,
                            best_betas_test[subset_ids],
                            ind.col = info_snp_subset$`_NUM_ID_`)
  
  test_gene_set_scores <- cbind(test_gene_set_scores, test_set_scores)
}

colnames(test_gene_set_scores) <- paste0('set_', 1:64)
test_gene_set_scores$full_prs <- best_preds_test
test_gene_set_scores$sex <- test_y$Sex
test_gene_set_scores$height <- test_y$height
write.csv(test_gene_set_scores, 'test_data_set_scores.csv', row.names=FALSE)


#--------------------------------------------------------------------------------
#TEST SET EVALUATION
#--------------------------------------------------------------------------------
#Full training set
train_y <- y[,c('height','Sex','LDPRED2_0.841982_0.3267_TRUE', 'ct_0.5', 'LDPRED2_auto', 'LDPRED2_inf', 'lasso_0.00065_0.005')] 

#Generating principal components for full training set
train_y$PC1 <- predict(test_pca$obj.svd.ref)[,1] 
train_y$PC2 <- predict(test_pca$obj.svd.ref)[,2] 
train_y$PC3 <- predict(test_pca$obj.svd.ref)[,3] 
train_y$PC4 <- predict(test_pca$obj.svd.ref)[,4] 
train_y$PC5 <- predict(test_pca$obj.svd.ref)[,5] 

write.csv(train_y, 'full_train_set.csv', row.names=FALSE)

#Null Model
full_null_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5")), data = train_y)
test_null_preds <- predict(full_null_model, test_y)
test_null_mse <- (mean((test_y$height - test_null_preds)^2))
test_null_r2 <- R2_Score(test_null_preds, test_y$height)

#LDpred2-grid model
full_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5 + LDPRED2_0.841982_0.3267_TRUE")), data = train_y)
test_preds <- predict(full_model, test_y)
test_mse <- (mean((test_y$height - test_preds)^2))
test_r2 <- R2_Score(test_preds, test_y$height)

#LDpred2-auto model
auto_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5 + LDPRED2_auto")), data = train_y)
auto_test_preds <- predict(auto_model, test_y)
auto_test_mse <- (mean((test_y$height - auto_test_preds)^2))
auto_test_r2 <- R2_Score(auto_test_preds, test_y$height)

#LDpred2-inf model
inf_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5 + LDPRED2_inf")), data = train_y)
inf_test_preds <- predict(inf_model, test_y)
inf_test_mse <- (mean((test_y$height - inf_test_preds)^2))
inf_test_r2 <- R2_Score(inf_test_preds, test_y$height)

#Best C+T model
ct_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5 + ct_0.5")), data = train_y)
ct_test_preds <- predict(ct_model, test_y)
ct_test_mse <- (mean((test_y$height - ct_test_preds)^2))
ct_test_r2 <- R2_Score(ct_test_preds, test_y$height)

#Best Lassosum model
lasso_model <- lm(as.formula(paste0("height~Sex + PC1 + PC2 + PC3 + PC4 + PC5 + lasso_0.00065_0.005")), data = train_y)
lasso_test_preds <- predict(lasso_model, test_y)
lasso_test_mse <- (mean((test_y$height - lasso_test_preds)^2))
lasso_test_r2 <- R2_Score(lasso_test_preds, test_y$height)

#Creating and exporting test set results
test_results <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(test_results) <- c('model', 'test_r2', 'test_mse')
test_results <- rbindlist(list(test_results,list('null_model', test_null_r2, test_null_mse)))
test_results <- rbindlist(list(test_results,list('LDPRED2_grid', test_r2, test_mse)))
test_results <- rbindlist(list(test_results,list('LDPRED2_auto', auto_test_r2, auto_test_mse)))
test_results <- rbindlist(list(test_results,list('LDPRED2_inf', inf_test_r2, inf_test_mse)))
test_results <- rbindlist(list(test_results,list('ct_0.5', ct_test_r2, ct_test_mse)))
test_results <- rbindlist(list(test_results,list('lasso_0.00065_0.005', lasso_test_r2, lasso_test_mse)))

write.csv(test_results, 'base_test_results.csv', row.names=FALSE)

#Exporting test set predictions for statistical significance test
test_pred_df <- data.frame(matrix(ncol=0, nrow=129))
test_pred_df$full_prs <- test_preds
test_pred_df$actual <- test_y$height
write.csv(test_pred_df, 'base_preds.csv', row.names=FALSE)
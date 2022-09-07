#Set working directory
setwd('ENTER DIRECTORY PATH')

#Import libraries
library(data.table)
library(dplyr)
library(psych)
set.seed(1234)

#Reading in model test set predictions
preds_py <- fread('all_preds.csv')
preds_r <- fread('base_preds.csv')

#Williams Test of Best Results
corr_enet_true <- cor(preds_py$pred_enet, preds_py$actual)
corr_prs_true <- cor(preds_r$full_prs, preds_py$actual)
corr_prs_enet <- cor(preds_r$full_prs, preds_py$pred_enet)
test <- paired.r(corr_enet_true, corr_prs_true, corr_prs_enet, n=129, twotailed=FALSE)

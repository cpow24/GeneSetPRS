#Setting working directory
setwd('ENTER DIRECTORY PATH')

#Loading required packages
library(data.table)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(ramify)
library(grid)
library(cowplot)

#Gene set titles
title_list = list()

#Gene set SNP count
snp_count = list()

#Gene set SNPs
snp_list = list()

#Set ID List
set_id_list <- list()

for (i in 1:64){
  snp_set <- fread(paste0('snplist_set_',i,'.set'))
  snp_set <- head(snp_set, -2)
  title_list <- append(title_list, colnames(snp_set)[1])
  snp_count <- append(snp_count, nrow(snp_set))
  snp_list <- append(snp_list, snp_set)
  set_id_list <- append(set_id_list, paste0('set_',i))
}

#Calculating number of overlapping SNPs in each gene set
common <- outer(snp_list, snp_list, Vectorize(function(x, y) sum(x %in% y)))

#Calculating Mean % Overlap for each gene set
mean_overlap_list <- list()

for (i in 1:64){

  prop_overlap <- (sum(common[i,-i] / as.numeric(snp_count[i])) / 63)
  mean_overlap_list <- append(mean_overlap_list, prop_overlap)
  
}

gene_set_df <- data.frame(matrix(ncol=0, nrow=64))
gene_set_df$set_name <- as.character(title_list)
gene_set_df$snp_count <- snp_count
gene_set_df$mean_overlap <- mean_overlap_list
gene_set_df$set <- set_id_list
gene_set_df <- apply(gene_set_df,2,as.character)
gene_set_df <- data.frame(gene_set_df)
gene_set_df$snp_count <- as.numeric(gene_set_df$snp_count)
gene_set_df$mean_overlap <- as.numeric(gene_set_df$mean_overlap)

options(scipen=999)

ggplot(gene_set_df, aes(x=set, y=mean_overlap)) +
  geom_bar (stat="identity", color="darkblue", fill="dodgerblue") +
  labs(title="Mean Gene Set Overlap (%)",
       x ="Gene Set", y = "Mean % Overlap") +
  theme(
    plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5, vjust=2),
    axis.title.x = element_text(color="Black", size=13, face="bold"),
    axis.title.y = element_text(color="Black", size=14, angle=90, face="bold"),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE
  ) +
  coord_flip() +
  scale_x_discrete(position = "bottom", limits  = unique(gene_set_df[["set"]])) +
  scale_y_continuous(expand=c(0,0),labels= scales::percent, limits=c(0,0.4))

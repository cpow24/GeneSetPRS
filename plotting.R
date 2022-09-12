#Setting working directory
setwd('ENTER DIRECTORY PATH')

#Loading required packages
library(data.table)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(cowplot)


#Plotting distribution of height
phenotype <- fread("target_train.height")
ggplot(phenotype, aes(x=height))+
  geom_histogram(color="darkblue", fill="dodgerblue", bins=25)+
  labs(title="Distribution of Height",
       x ="Height", y = "Count") +
  theme(
    plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5, vjust=2),
    axis.title.x = element_text(color="Black", size=14, face="bold"),
    axis.title.y = element_text(color="Black", size=14, angle=90, face="bold"),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                    colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
                                    panel.grid.minor = element_line(size = rel(0.5)), 
                                    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                     colour = NA), complete = TRUE
  )


#Plotting F coefficients from Inbreeding Coefficient (females = < 0.2, males = > 0.8)
sexcheck_train <-fread("target_train_sex.sexcheck")
sexcheck_test <-fread("target_test_sex.sexcheck")
full_sex <- rbind(sexcheck_train, sexcheck_test)
full_sex <- full_sex %>% 
  mutate(Sex= if_else(F <= 0.2, 2, if_else(F >= 0.8, 1, 0)))
full_sex <- full_sex %>% 
  mutate(sex_label= if_else(F <= 0.2, 'Female', if_else(F >= 0.8, 'Male', 'Unknown')))

ggplot(full_sex, aes(x=F, fill=sex_label))+
  geom_histogram(bins=50, alpha=0.8)+
  labs(title="Imputed Sex",
       x ="F Coefficient", y = "Count", fill="Sex") +
  theme(
    plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5, vjust=2),
    axis.title.x = element_text(color="Black", size=14, face="bold"),
    axis.title.y = element_text(color="Black", size=14, angle=90, face="bold"),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"),
    complete = TRUE,
    legend.position = c(0.1, 0.85),
    legend.title = element_text(color = "black", size = 14, face="bold"),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = "grey92"),
    legend.key = element_rect(fill = "grey92", color = NA)
  )
table(full_sex$sex_label)

#Plotting QQ-Plot of Polygenic Scores
risk_scores <- read.csv('full_train_set.csv')
qqnorm(risk_scores$LDPRED2_0.841982_0.3267_TRUE, pch = 1, frame = FALSE)
qqline(risk_scores$LDPRED2_0.841982_0.3267_TRUE, col = "dodgerblue", lwd = 2)
box(col='black')

#Plotting results of feature selection
feat_res <- read.csv('results_feat_num.csv')

ggplot(data=feat_res, aes(x=num_features, y=r2, group=model)) +
  geom_line(aes(color=model),size=1.3)+
  geom_point(size=1.8)+
  labs(title="MRMR Feature Selection",
       x ="Number of Features", y = "R2") +
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
    legend.position = c(0.1, 0.2),
    legend.title = element_text(color = "black", size = 14, face="bold"),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = "grey92"),
    legend.key = element_rect(fill = "grey92", color = NA)
  )


#Plotting Gene Set Polygenic Scores
set_scores <- read.csv('set_scores.csv')
list <-lapply(1:64,
              function(col) ggplot2::ggplot(set_scores, aes(x=set_scores[[paste0('set_',col)]])) + geom_histogram(color="darkblue", fill="dodgerblue", bins=25)
              +labs(title=' ',
                     x ="Score", y = "Count") 
              +theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(), 
                      axis.text.y=element_blank(),  
                      axis.ticks.y=element_blank(),
                     plot.title = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                ))


title <- ggdraw() + 
  draw_label(
    "'Distribution of Gene Set Polygenic Risk Scores",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

cowplot::plot_grid(plotlist=list, labels=c(paste0(1:64)), scale=1.1)

#Plotting Features and Rules from RuleFit Model
rf_df <- read.csv('rfit_rules_labelled.csv')
rules <- subset(rf_df, coef != 0)

#Disabling scientific notation
options(scipen=999)

ggplot(rules[1:10,], aes(x=reorder(label, importance), y=importance)) +
  geom_bar (stat="identity", color="darkblue", fill="dodgerblue") +
  labs(title="RuleFit Feature Importance",
       x ="Feature", y = "Importance") +
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
                                                                                  colour = NA), complete = TRUE) +
  coord_flip()


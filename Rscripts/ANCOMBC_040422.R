setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/phyloseq/")
library(rRDPData)
library(debar)
library(tidyverse)
library(Rcpp)
library(plyr)
library("phyloseq")
library(ggplot2)
library(vegan)
library(ANCOMBC)
library(gridExtra)
library(microbiome)
library(stringr)
library(ggpubr)
library(cowplot)


#Load in the phyloseq object 
load('phylo_main_no_cont_NoOutliers.RData')



#Aggregate to genus lvl
phylo_main_no_cont_NoOutliers_Genus <- aggregate_taxa(phylo_main_no_cont_NoOutliers, "Genus")

#Filtering all taxa wiuth read sum < 100
phylo_main_no_cont_NoOutliers_Genus_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_Genus, function (x) {sum(x > 100) > 0}, prune=TRUE)

#summary(colSums(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter)))

#Remove Phaeobacter OTU before making differential analyses
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_Genus_filter, Genus != "Phaeobacter_inhibens")


#Planktonic suspension Control set as REF  ----
phylo_main_no_cont_NoOutliers_SW <- subset_samples(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_SW_day1 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="1")
phylo_main_no_cont_NoOutliers_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="4")
phylo_main_no_cont_NoOutliers_SW_day10 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="10")


#Order Treatments so WT is the comparing group
sample_data(phylo_main_no_cont_NoOutliers_SW_day1)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_day1)$Treatment, levels = c("Control", "WT", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_SW_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_day4)$Treatment, levels = c("Control", "WT", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_SW_day10)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_day10)$Treatment, levels = c("Control", "WT", "dTDA"))


out_phylo_main_no_cont_NoOutliers_SW_Genus_day1 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_day1, formula = "Treatment",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_SW_Genus_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_day4, formula = "Treatment",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_SW_Genus_day10 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_day10, formula = "Treatment",
                                                            p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                            group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                            tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                            alpha = 0.01, global = TRUE)

#Results SW Day1
res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1 <- out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$res
#Results SW Day4
res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4 <- out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$res
#Results SW Day10
res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10 <- out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$res

#### Dataframe Day 1
df_fig1_SW_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$se * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day1)[-1] = paste0(colnames(df_fig2_SW_day1)[-1], "SD")
colnames(df_fig1_SW_day1) <- sub("Treatment", "", colnames(df_fig1_SW_day1))
colnames(df_fig2_SW_day1) <- sub("Treatment", "", colnames(df_fig2_SW_day1))

df_fig_SW_day1_dat = df_fig1_SW_day1 %>% left_join(df_fig2_SW_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day1 = df_fig1_SW_day1 %>% left_join(df_fig2_SW_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day1 <- as.data.frame(df_fig_SW_day1)
df_fig_SW_day1 <- df_fig_SW_day1[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day1 <- reshape(df_fig_SW_day1,  direction='long', 
                          varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                          timevar='var',
                          times=c('WT', 'dTDA'),
                          v.names=c('LogFoldChange', 'SD'))

df_fig_SW_day1 <- df_fig_SW_day1 %>% mutate(group = ifelse(LogFoldChange > 0, "WT & dTDA", "Control"))
str(df_fig_SW_day1)
#df_fig_SW_day1$taxon_id = factor(df_fig_SW_day1$taxon_id, levels = df_fig_SW_day1$taxon_id)
df_fig_SW_day1$Day <- rep("Day 1")

#Count negative and positive values
df_fig_SW_day1_sum <- df_fig_SW_day1 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 4
df_fig1_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$se * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day4)[-1] = paste0(colnames(df_fig2_SW_day4)[-1], "SD")
colnames(df_fig1_SW_day4) <- sub("Treatment", "", colnames(df_fig1_SW_day4))
colnames(df_fig2_SW_day4) <- sub("Treatment", "", colnames(df_fig2_SW_day4))

df_fig_SW_day4_dat = df_fig1_SW_day4 %>% left_join(df_fig2_SW_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day4 = df_fig1_SW_day4 %>% left_join(df_fig2_SW_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_SW_day4 <- as.data.frame(df_fig_SW_day4)
df_fig_SW_day4 <- df_fig_SW_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day4 <- reshape(df_fig_SW_day4, direction='long', 
                          varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                          timevar='var',
                          times=c('WT', 'dTDA'),
                          v.names=c('LogFoldChange', 'SD'))

df_fig_SW_day4 <- df_fig_SW_day4 %>% mutate(group = ifelse(LogFoldChange > 0,  "WT & dTDA", "Control"))
str(df_fig_SW_day4)
#df_fig_SW_day4$taxon_id = factor(df_fig_SW_day4$taxon_id, levels = df_fig_SW_day4$taxon_id)
df_fig_SW_day4$Day <- rep("Day 4")

#Count negative and positive values
df_fig_SW_day4_sum <- df_fig_SW_day4 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 10
df_fig1_SW_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$se * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day10)[-1] = paste0(colnames(df_fig2_SW_day10)[-1], "SD")
colnames(df_fig1_SW_day10) <- sub("Treatment", "", colnames(df_fig1_SW_day10))
colnames(df_fig2_SW_day10) <- sub("Treatment", "", colnames(df_fig2_SW_day10))

df_fig_SW_day10_dat = df_fig1_SW_day10 %>% left_join(df_fig2_SW_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day10 = df_fig1_SW_day10 %>% left_join(df_fig2_SW_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_SW_day10 <- as.data.frame(df_fig_SW_day10)
df_fig_SW_day10 <- df_fig_SW_day10[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day10 <- reshape(df_fig_SW_day10,  direction='long', 
                           varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                           timevar='var',
                           times=c('WT', 'dTDA'),
                           v.names=c('LogFoldChange', 'SD'))

df_fig_SW_day10 <- df_fig_SW_day10 %>% mutate(group = ifelse(LogFoldChange > 0,  "WT & dTDA", "Control"))
str(df_fig_SW_day10)
#df_fig_SW_day10$taxon_id = factor(df_fig_SW_day10$taxon_id, levels = df_fig_SW_day10$taxon_id)
df_fig_SW_day10$Day <- rep("Day 10")

#Count negative and positive values
df_fig_SW_day10_sum <- df_fig_SW_day10 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))


df_fig_SW_day1_sum_dat <- as.data.frame(df_fig_SW_day1_sum)
#rownames(df_fig_SW_day1_sum_dat) <- c("Control", "dTDA")
df_fig_SW_day1_sum_dat$Day <- rep("Day 1")
df_fig_SW_day4_sum_dat <- as.data.frame(df_fig_SW_day4_sum)
#rownames(df_fig_SW_day4_sum_dat) <- c("Control", "dTDA")
df_fig_SW_day4_sum_dat$Day <- rep("Day 4")
df_fig_SW_day10_sum_dat <- as.data.frame(df_fig_SW_day10_sum)
#rownames(df_fig_SW_day10_sum_dat) <- c("Control", "dTDA")
df_fig_SW_day10_sum_dat$Day <- rep("Day 10")



df_fig_SW_sum_dat_all <- rbind(df_fig_SW_day1_sum_dat,df_fig_SW_day4_sum_dat,df_fig_SW_day10_sum_dat)
colnames(df_fig_SW_sum_dat_all) <- c("Treatment", "LogFC > 0", "LogFC < 0", "Day")
library(reshape)
df_fig_SW_sum_dat_all_long <- reshape(df_fig_SW_sum_dat_all, direction='long',
                                      varying=c('LogFC > 0', 'LogFC < 0'),
                                      timevar='var',
                                      times=c('LogFC > 0', 'LogFC < 0'),
                                      v.names=c('Count'))

df_fig_SW_sum_dat_all_long$Day <- factor(df_fig_SW_sum_dat_all_long$Day, levels = c("Day 1", "Day 4", "Day 10"))


#Biofilm Control set as REF ----
phylo_main_no_cont_NoOutliers_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_Biofilm_day1 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="1")
phylo_main_no_cont_NoOutliers_Biofilm_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="4")
phylo_main_no_cont_NoOutliers_Biofilm_day10 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="10")


#Order Treatments so WT is the comparing group
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day1)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day1)$Treatment, levels = c("Control", "WT", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day4)$Treatment, levels = c("Control", "WT", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day10)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day10)$Treatment, levels = c("Control", "WT", "dTDA"))


out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_day1, formula = "Treatment",
                                                                p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_day4, formula = "Treatment",
                                                                p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_day10, formula = "Treatment",
                                                                 p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                 group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                 tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                 alpha = 0.01, global = TRUE)

#Results Biofilm Day1
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$res
#Results Biofilm Day4
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$res
#Results Biofilm Day10
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$res

#### Dataframe Day 1
df_fig1_Biofilm_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day1)[-1] = paste0(colnames(df_fig2_Biofilm_day1)[-1], "SD")
colnames(df_fig1_Biofilm_day1) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day1))
colnames(df_fig2_Biofilm_day1) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day1))

df_fig_Biofilm_day1_dat = df_fig1_Biofilm_day1 %>% left_join(df_fig2_Biofilm_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day1 = df_fig1_Biofilm_day1 %>% left_join(df_fig2_Biofilm_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day1 <- as.data.frame(df_fig_Biofilm_day1)
df_fig_Biofilm_day1 <- df_fig_Biofilm_day1[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day1 <- reshape(df_fig_Biofilm_day1,  direction='long', 
                               varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                               timevar='var',
                               times=c('WT', 'dTDA'),
                               v.names=c('LogFoldChange', 'SD'))

df_fig_Biofilm_day1 <- df_fig_Biofilm_day1 %>% mutate(group = ifelse(LogFoldChange > 0, "WT & dTDA", "Control"))
str(df_fig_Biofilm_day1)
#df_fig_Biofilm_day1$taxon_id = factor(df_fig_Biofilm_day1$taxon_id, levels = df_fig_Biofilm_day1$taxon_id)
df_fig_Biofilm_day1$Day <- rep("Day 1")

#Count negative and positive values
df_fig_Biofilm_day1_sum <- df_fig_Biofilm_day1 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 4
df_fig1_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day4)[-1] = paste0(colnames(df_fig2_Biofilm_day4)[-1], "SD")
colnames(df_fig1_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day4))
colnames(df_fig2_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day4))

df_fig_Biofilm_day4_dat = df_fig1_Biofilm_day4 %>% left_join(df_fig2_Biofilm_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day4 = df_fig1_Biofilm_day4 %>% left_join(df_fig2_Biofilm_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_Biofilm_day4 <- as.data.frame(df_fig_Biofilm_day4)
df_fig_Biofilm_day4 <- df_fig_Biofilm_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day4 <- reshape(df_fig_Biofilm_day4, direction='long', 
                               varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                               timevar='var',
                               times=c('WT', 'dTDA'),
                               v.names=c('LogFoldChange', 'SD'))

df_fig_Biofilm_day4 <- df_fig_Biofilm_day4 %>% mutate(group = ifelse(LogFoldChange > 0,  "WT & dTDA", "Control"))
str(df_fig_Biofilm_day4)
#df_fig_Biofilm_day4$taxon_id = factor(df_fig_Biofilm_day4$taxon_id, levels = df_fig_Biofilm_day4$taxon_id)
df_fig_Biofilm_day4$Day <- rep("Day 4")

#Count negative and positive values
df_fig_Biofilm_day4_sum <- df_fig_Biofilm_day4 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 10
df_fig1_Biofilm_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day10)[-1] = paste0(colnames(df_fig2_Biofilm_day10)[-1], "SD")
colnames(df_fig1_Biofilm_day10) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day10))
colnames(df_fig2_Biofilm_day10) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day10))

df_fig_Biofilm_day10_dat = df_fig1_Biofilm_day10 %>% left_join(df_fig2_Biofilm_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day10 = df_fig1_Biofilm_day10 %>% left_join(df_fig2_Biofilm_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_Biofilm_day10 <- as.data.frame(df_fig_Biofilm_day10)
df_fig_Biofilm_day10 <- df_fig_Biofilm_day10[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day10 <- reshape(df_fig_Biofilm_day10,  direction='long', 
                                varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                                timevar='var',
                                times=c('WT', 'dTDA'),
                                v.names=c('LogFoldChange', 'SD'))

df_fig_Biofilm_day10 <- df_fig_Biofilm_day10 %>% mutate(group = ifelse(LogFoldChange > 0,  "WT & dTDA", "Control"))
str(df_fig_Biofilm_day10)
#df_fig_Biofilm_day10$taxon_id = factor(df_fig_Biofilm_day10$taxon_id, levels = df_fig_Biofilm_day10$taxon_id)
df_fig_Biofilm_day10$Day <- rep("Day 10")

#Count negative and positive values
df_fig_Biofilm_day10_sum <- df_fig_Biofilm_day10 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))


df_fig_Biofilm_day1_sum_dat <- as.data.frame(df_fig_Biofilm_day1_sum)
#rownames(df_fig_Biofilm_day1_sum_dat) <- c("Control", "dTDA")
df_fig_Biofilm_day1_sum_dat$Day <- rep("Day 1")
df_fig_Biofilm_day4_sum_dat <- as.data.frame(df_fig_Biofilm_day4_sum)
#rownames(df_fig_Biofilm_day4_sum_dat) <- c("Control", "dTDA")
df_fig_Biofilm_day4_sum_dat$Day <- rep("Day 4")
df_fig_Biofilm_day10_sum_dat <- as.data.frame(df_fig_Biofilm_day10_sum)
#rownames(df_fig_Biofilm_day10_sum_dat) <- c("Control", "dTDA")
df_fig_Biofilm_day10_sum_dat$Day <- rep("Day 10")



df_fig_Biofilm_sum_dat_all <- rbind(df_fig_Biofilm_day1_sum_dat,df_fig_Biofilm_day4_sum_dat,df_fig_Biofilm_day10_sum_dat)
colnames(df_fig_Biofilm_sum_dat_all) <- c("Treatment", "LogFC > 0", "LogFC < 0", "Day")
library(reshape)
df_fig_Biofilm_sum_dat_all_long <- reshape(df_fig_Biofilm_sum_dat_all, direction='long',
                                           varying=c('LogFC > 0', 'LogFC < 0'),
                                           timevar='var',
                                           times=c('LogFC > 0', 'LogFC < 0'),
                                           v.names=c('Count'))

df_fig_Biofilm_sum_dat_all_long$Day <- factor(df_fig_Biofilm_sum_dat_all_long$Day, levels = c("Day 1", "Day 4", "Day 10"))


###### ANCOM-BC Log FC counts Control as REFERENCE GROUP------
df_fig_SW_sum_dat_all_long$Env <- rep("Planktonic suspension") 
df_fig_Biofilm_sum_dat_all_long$Env <- rep("Biofilm")

LogFoldCounts_sum_dat_all_long_dTDA <- rbind(df_fig_Biofilm_sum_dat_all_long,df_fig_SW_sum_dat_all_long)
sum(LogFoldCounts_sum_dat_all_long_dTDA$Count)
LogFoldCounts_sum_dat_all_long_dTDA %>%
  group_by(var) %>% dplyr::summarise(Count = sum(Count)) %>% mutate(percentage = Count / sum(LogFoldCounts_sum_dat_all_long_dTDA$Count)*100)

LogFoldCounts_sum_dat_all_long_dTDA %>%
  group_by(var, Treatment) %>% dplyr::summarise(Count = sum(Count)) 



#Planktonic suspension WT set as REF ----
phylo_main_no_cont_NoOutliers_SW <- subset_samples(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_SW_day1 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="1")
phylo_main_no_cont_NoOutliers_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="4")
phylo_main_no_cont_NoOutliers_SW_day10 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="10")


#Order Treatments so WT is the comparing group
sample_data(phylo_main_no_cont_NoOutliers_SW_day1)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_day1)$Treatment, levels = c("WT", "Control", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_SW_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_day4)$Treatment, levels = c("WT", "Control", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_SW_day10)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_day10)$Treatment, levels = c("WT", "Control", "dTDA"))


out_phylo_main_no_cont_NoOutliers_SW_Genus_day1 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_day1, formula = "Treatment",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_SW_Genus_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_day4, formula = "Treatment",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_SW_Genus_day10 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_day10, formula = "Treatment",
                                                            p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                            group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                            tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                            alpha = 0.01, global = TRUE)

#Results SW Day1
res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1 <- out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$res
#Results SW Day4
res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4 <- out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$res
#Results SW Day10
res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10 <- out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$res

#### Dataframe Day 1
df_fig1_SW_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$se * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day1)[-1] = paste0(colnames(df_fig2_SW_day1)[-1], "SD")
colnames(df_fig1_SW_day1) <- sub("Treatment", "", colnames(df_fig1_SW_day1))
colnames(df_fig2_SW_day1) <- sub("Treatment", "", colnames(df_fig2_SW_day1))

df_fig_SW_day1_dat = df_fig1_SW_day1 %>% left_join(df_fig2_SW_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day1 = df_fig1_SW_day1 %>% left_join(df_fig2_SW_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day1 <- as.data.frame(df_fig_SW_day1)
df_fig_SW_day1 <- df_fig_SW_day1[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day1 <- reshape(df_fig_SW_day1, direction='long', 
                          varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                          timevar='var',
                          times=c('Control', 'dTDA'),
                          v.names=c('LogFoldChange', 'SD'))

df_fig_SW_day1 <- df_fig_SW_day1 %>% mutate(group = ifelse(LogFoldChange > 0, "Control & dTDA", "WT"))
str(df_fig_SW_day1)
#df_fig_SW_day1$taxon_id = factor(df_fig_SW_day1$taxon_id, levels = df_fig_SW_day1$taxon_id)
df_fig_SW_day1$Day <- rep("Day 1")

#Count negative and positive values
df_fig_SW_day1_sum <- df_fig_SW_day1 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 4
df_fig1_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$se * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day4)[-1] = paste0(colnames(df_fig2_SW_day4)[-1], "SD")
colnames(df_fig1_SW_day4) <- sub("Treatment", "", colnames(df_fig1_SW_day4))
colnames(df_fig2_SW_day4) <- sub("Treatment", "", colnames(df_fig2_SW_day4))

df_fig_SW_day4_dat = df_fig1_SW_day4 %>% left_join(df_fig2_SW_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day4 = df_fig1_SW_day4 %>% left_join(df_fig2_SW_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_SW_day4 <- as.data.frame(df_fig_SW_day4)
df_fig_SW_day4 <- df_fig_SW_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day4 <- reshape(df_fig_SW_day4, direction='long', 
                          varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                          timevar='var',
                          times=c('Control', 'dTDA'),
                          v.names=c('LogFoldChange', 'SD'))

df_fig_SW_day4 <- df_fig_SW_day4 %>% mutate(group = ifelse(LogFoldChange > 0, "Control & dTDA", "WT"))
str(df_fig_SW_day4)
#df_fig_SW_day4$taxon_id = factor(df_fig_SW_day4$taxon_id, levels = df_fig_SW_day4$taxon_id)
df_fig_SW_day4$Day <- rep("Day 4")

#Count negative and positive values
df_fig_SW_day4_sum <- df_fig_SW_day4 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 10
df_fig1_SW_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$se * res_out_phylo_main_no_cont_NoOutliers_SW_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day10)[-1] = paste0(colnames(df_fig2_SW_day10)[-1], "SD")
colnames(df_fig1_SW_day10) <- sub("Treatment", "", colnames(df_fig1_SW_day10))
colnames(df_fig2_SW_day10) <- sub("Treatment", "", colnames(df_fig2_SW_day10))

df_fig_SW_day10_dat = df_fig1_SW_day10 %>% left_join(df_fig2_SW_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_SW_day10 = df_fig1_SW_day10 %>% left_join(df_fig2_SW_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_SW_day10 <- as.data.frame(df_fig_SW_day10)
df_fig_SW_day10 <- df_fig_SW_day10[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day10 <- reshape(df_fig_SW_day10, direction='long', 
                           varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                           timevar='var',
                           times=c('Control', 'dTDA'),
                           v.names=c('LogFoldChange', 'SD'))

df_fig_SW_day10 <- df_fig_SW_day10 %>% mutate(group = ifelse(LogFoldChange > 0, "Control & dTDA", "WT"))
str(df_fig_SW_day10)
#df_fig_SW_day10$taxon_id = factor(df_fig_SW_day10$taxon_id, levels = df_fig_SW_day10$taxon_id)
df_fig_SW_day10$Day <- rep("Day 10")

#Count negative and positive values
df_fig_SW_day10_sum <- df_fig_SW_day10 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))


df_fig_SW_day1_sum_dat <- as.data.frame(df_fig_SW_day1_sum)
#rownames(df_fig_SW_day1_sum_dat) <- c("Control", "dTDA")
df_fig_SW_day1_sum_dat$Day <- rep("Day 1")
df_fig_SW_day4_sum_dat <- as.data.frame(df_fig_SW_day4_sum)
#rownames(df_fig_SW_day4_sum_dat) <- c("Control", "dTDA")
df_fig_SW_day4_sum_dat$Day <- rep("Day 4")
df_fig_SW_day10_sum_dat <- as.data.frame(df_fig_SW_day10_sum)
#rownames(df_fig_SW_day10_sum_dat) <- c("Control", "dTDA")
df_fig_SW_day10_sum_dat$Day <- rep("Day 10")



df_fig_SW_sum_dat_all <- rbind(df_fig_SW_day1_sum_dat,df_fig_SW_day4_sum_dat,df_fig_SW_day10_sum_dat)
colnames(df_fig_SW_sum_dat_all) <- c("Treatment", "LogFC > 0", "LogFC < 0", "Day")
library(reshape)
df_fig_SW_sum_dat_all_long <- reshape(df_fig_SW_sum_dat_all, direction='long',
                                      varying=c('LogFC > 0', 'LogFC < 0'),
                                      timevar='var',
                                      times=c('LogFC > 0', 'LogFC < 0'),
                                      v.names=c('Count'))

df_fig_SW_sum_dat_all_long$Day <- factor(df_fig_SW_sum_dat_all_long$Day, levels = c("Day 1", "Day 4", "Day 10"))



# 
# 
# 
# 
# rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10, df_fig_SW_day1, df_fig_SW_day4, df_fig_SW_day10) %>% filter(var!="Control") %>% filter(LogFoldChange<0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
# rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10, df_fig_SW_day1, df_fig_SW_day4, df_fig_SW_day10) %>% filter(var!="Control") %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
# 
# 
# rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10, df_fig_SW_day1, df_fig_SW_day4, df_fig_SW_day10) %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
# 
# 
# rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10, df_fig_SW_day1, df_fig_SW_day4, df_fig_SW_day10) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
# rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10, df_fig_SW_day1, df_fig_SW_day4, df_fig_SW_day10) %>% filter(LogFoldChange<0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
# 

df_fig_Biofilm_day1 %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
df_fig_Biofilm_day4 %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
df_fig_Biofilm_day10 %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())

rbind(df_fig_SW_day1, df_fig_SW_day4, df_fig_SW_day10) %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
df_fig_SW_day1 %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
df_fig_SW_day4 %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())
df_fig_SW_day10 %>% group_by(var) %>% filter(LogFoldChange>0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())



rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10) %>% group_by(var) %>% filter(LogFoldChange<0) %>% dplyr::summarise(unique(taxon_id)) %>% dplyr::summarise(n())



rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10) %>% group_by(var) %>% filter(LogFoldChange>0) %>%  dplyr::summarise(n())
rbind(df_fig_Biofilm_day1, df_fig_Biofilm_day4, df_fig_Biofilm_day10) %>% group_by(var) %>% filter(group=="Control") %>%  dplyr::summarise(n())




#Biofilm WT set as REF  ----

phylo_main_no_cont_NoOutliers_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_Biofilm_day1 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="1")
phylo_main_no_cont_NoOutliers_Biofilm_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="4")
phylo_main_no_cont_NoOutliers_Biofilm_day10 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="10")

#Order Treatments so WT is the comparing group
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day1)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day1)$Treatment, levels = c("WT", "Control", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day4)$Treatment, levels = c("WT", "Control", "dTDA"))
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day10)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_day10)$Treatment, levels = c("WT", "Control", "dTDA"))


out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_day1, formula = "Treatment",
                                                                p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_day4, formula = "Treatment",
                                                                p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                alpha = 0.01, global = TRUE)

out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_day10, formula = "Treatment",
                                                                 p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                 group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                 tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                 alpha = 0.01, global = TRUE)


#Results Biofilm Day1
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$res
#Results Biofilm Day4
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$res


#Results Biofilm Day10
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$res

#### Dataframe Day 1
df_fig1_Biofilm_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day1 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day1$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day1)[-1] = paste0(colnames(df_fig2_Biofilm_day1)[-1], "SD")
colnames(df_fig1_Biofilm_day1) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day1))
colnames(df_fig2_Biofilm_day1) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day1))

df_fig_Biofilm_day1_dat = df_fig1_Biofilm_day1 %>% left_join(df_fig2_Biofilm_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day1 = df_fig1_Biofilm_day1 %>% left_join(df_fig2_Biofilm_day1, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day1 <- as.data.frame(df_fig_Biofilm_day1)
df_fig_Biofilm_day1 <- df_fig_Biofilm_day1[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day1 <- reshape(df_fig_Biofilm_day1, direction='long', 
                               varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                               timevar='var',
                               times=c('Control', 'dTDA'),
                               v.names=c('LogFoldChange', 'SD'))

df_fig_Biofilm_day1 <- df_fig_Biofilm_day1 %>% mutate(group = ifelse(LogFoldChange > 0, "Control & dTDA", "WT"))
str(df_fig_Biofilm_day1)
#df_fig_Biofilm_day1$taxon_id = factor(df_fig_Biofilm_day1$taxon_id, levels = df_fig_Biofilm_day1$taxon_id)
df_fig_Biofilm_day1$Day <- rep("Day 1")

#Count negative and positive values
df_fig_Biofilm_day1_sum <- df_fig_Biofilm_day1 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 4
df_fig1_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day4)[-1] = paste0(colnames(df_fig2_Biofilm_day4)[-1], "SD")
colnames(df_fig1_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day4))
colnames(df_fig2_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day4))


#Save df_fig_Biofilm_day4 for input for later network combined analysis in NetCoMi_030122.R
# save.image("ANCOMBC_030122")  # save workspace to disk
# rm(list=ls()) # remove everything from workspace
# tmp.env <- new.env() # create a temporary environment
# load("ANCOMBC_030122", envir=tmp.env) # load workspace into temporary environment
# df_fig1_Biofilm_day4 <- get("df_fig1_Biofilm_day4", pos=tmp.env) # get the objects you need into your globalenv()
# rm(tmp.env) # remove the temporary environment to free up memory
# save.image("df_fig1_Biofilm_day4")



df_fig_Biofilm_day4_dat = df_fig1_Biofilm_day4 %>% left_join(df_fig2_Biofilm_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day4 = df_fig1_Biofilm_day4 %>% left_join(df_fig2_Biofilm_day4, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day4 <- as.data.frame(df_fig_Biofilm_day4)
df_fig_Biofilm_day4 <- df_fig_Biofilm_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day4 <- reshape(df_fig_Biofilm_day4, direction='long', 
                               varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                               timevar='var',
                               times=c('Control', 'dTDA'),
                               v.names=c('LogFoldChange', 'SD'))

df_fig_Biofilm_day4 <- df_fig_Biofilm_day4 %>% mutate(group = ifelse(LogFoldChange > 0, "Control & dTDA", "WT"))
str(df_fig_Biofilm_day4)
#df_fig_Biofilm_day4$taxon_id = factor(df_fig_Biofilm_day4$taxon_id, levels = df_fig_Biofilm_day4$taxon_id)
df_fig_Biofilm_day4$Day <- rep("Day 4")



#Count negative and positive values
df_fig_Biofilm_day4_sum <- df_fig_Biofilm_day4 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))

#### Dataframe Day 10
df_fig1_Biofilm_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day10 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Genus_day10$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day10)[-1] = paste0(colnames(df_fig2_Biofilm_day10)[-1], "SD")
colnames(df_fig1_Biofilm_day10) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day10))
colnames(df_fig2_Biofilm_day10) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day10))

df_fig_Biofilm_day10_dat = df_fig1_Biofilm_day10 %>% left_join(df_fig2_Biofilm_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day10 = df_fig1_Biofilm_day10 %>% left_join(df_fig2_Biofilm_day10, by = "taxon_id") %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 


df_fig_Biofilm_day10 <- as.data.frame(df_fig_Biofilm_day10)
df_fig_Biofilm_day10 <- df_fig_Biofilm_day10[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day10 <- reshape(df_fig_Biofilm_day10, direction='long', 
                                varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                                timevar='var',
                                times=c('Control', 'dTDA'),
                                v.names=c('LogFoldChange', 'SD'))

df_fig_Biofilm_day10 <- df_fig_Biofilm_day10 %>% mutate(group = ifelse(LogFoldChange > 0, "Control & dTDA", "WT"))
str(df_fig_Biofilm_day10)
#df_fig_Biofilm_day10$taxon_id = factor(df_fig_Biofilm_day10$taxon_id, levels = df_fig_Biofilm_day10$taxon_id)
df_fig_Biofilm_day10$Day <- rep("Day 10")

#Count negative and positive values
df_fig_Biofilm_day10_sum <- df_fig_Biofilm_day10 %>%
  group_by(var) %>%
  dplyr::summarise(pos = sum(LogFoldChange>0),
                   neg = sum(LogFoldChange<0))


df_fig_Biofilm_day1_sum_dat <- as.data.frame(df_fig_Biofilm_day1_sum)
#rownames(df_fig_Biofilm_day1_sum_dat) <- c("Control", "dTDA")
df_fig_Biofilm_day1_sum_dat$Day <- rep("Day 1")
df_fig_Biofilm_day4_sum_dat <- as.data.frame(df_fig_Biofilm_day4_sum)
#rownames(df_fig_Biofilm_day4_sum_dat) <- c("Control", "dTDA")
df_fig_Biofilm_day4_sum_dat$Day <- rep("Day 4")
df_fig_Biofilm_day10_sum_dat <- as.data.frame(df_fig_Biofilm_day10_sum)
#rownames(df_fig_Biofilm_day10_sum_dat) <- c("Control", "dTDA")
df_fig_Biofilm_day10_sum_dat$Day <- rep("Day 10")



df_fig_Biofilm_sum_dat_all <- rbind(df_fig_Biofilm_day1_sum_dat,df_fig_Biofilm_day4_sum_dat,df_fig_Biofilm_day10_sum_dat)
colnames(df_fig_Biofilm_sum_dat_all) <- c("Treatment", "LogFC > 0", "LogFC < 0", "Day")
library(reshape)
df_fig_Biofilm_sum_dat_all_long <- reshape(df_fig_Biofilm_sum_dat_all, direction='long',
                                           varying=c('LogFC > 0', 'LogFC < 0'),
                                           timevar='var',
                                           times=c('LogFC > 0', 'LogFC < 0'),
                                           v.names=c('Count'))

df_fig_Biofilm_sum_dat_all_long$Day <- factor(df_fig_Biofilm_sum_dat_all_long$Day, levels = c("Day 1", "Day 4", "Day 10"))

###### ANCOM-BC Log FC counts WT as REFERENCE GROUP ------
df_fig_SW_sum_dat_all_long$Env <- rep("Planktonic suspension") 
df_fig_Biofilm_sum_dat_all_long$Env <- rep("Biofilm")

LogFoldCounts_sum_dat_all_long_WT <- rbind(df_fig_Biofilm_sum_dat_all_long,df_fig_SW_sum_dat_all_long)
sum(LogFoldCounts_sum_dat_all_long_WT$Count)
LogFoldCounts_sum_dat_all_long_WT %>%
  group_by(var) %>% dplyr::summarise(Count = sum(Count)) %>% mutate(percentage = Count / sum(LogFoldCounts_sum_dat_all_long_WT$Count)*100)

LogFoldCounts_sum_dat_all_long_WT %>%
  group_by(var, Treatment) %>% dplyr::summarise(Count = sum(Count)) 


LogFoldCounts_sum_dat_all_long_WT <- LogFoldCounts_sum_dat_all_long_WT[,c(1,2,4,3,6)]




##### COMBINING ALL COMBINATIONS -----

LogFoldCounts_sum_dat_all_long_dTDA <- LogFoldCounts_sum_dat_all_long_dTDA[,c(1,2,4,3,6)]
LogFoldCounts_sum_dat_all_long_dTDA_Control <- LogFoldCounts_sum_dat_all_long_dTDA %>% filter(Treatment == "dTDA") 
LogFoldCounts_sum_dat_all_long_dTDA_Control$Comp <- ifelse(LogFoldCounts_sum_dat_all_long_dTDA_Control$Treatment=="dTDA", "Control vs. dTDA", LogFoldCounts_sum_dat_all_long$Treatment)
colnames(LogFoldCounts_sum_dat_all_long_dTDA_Control)

LogFoldCounts_sum_dat_all_long_WT$Comp <- ifelse(LogFoldCounts_sum_dat_all_long_WT$Treatment=="Control", "WT vs. Control", LogFoldCounts_sum_dat_all_long_WT$Treatment)
LogFoldCounts_sum_dat_all_long_WT$Comp <- ifelse(LogFoldCounts_sum_dat_all_long_WT$Comp=="dTDA", "WT vs. dTDA", LogFoldCounts_sum_dat_all_long_WT$Comp)
colnames(LogFoldCounts_sum_dat_all_long_WT)


LogFoldCounts_sum_dat_all_long_ALLComp <- rbind(LogFoldCounts_sum_dat_all_long_WT, LogFoldCounts_sum_dat_all_long_dTDA_Control)
LogFoldCounts_sum_dat_all_long_ALLComp$Comp <- factor(LogFoldCounts_sum_dat_all_long_ALLComp$Comp, levels = c("WT vs. Control", "WT vs. dTDA", "Control vs. dTDA"))

#Change Days to factoral numbers
LogFoldCounts_sum_dat_all_long_ALLComp$Day <-  gsub('Day ','', LogFoldCounts_sum_dat_all_long_ALLComp$Day)

LogFoldCounts_sum_dat_all_long_ALLComp$Day <- factor(LogFoldCounts_sum_dat_all_long_ALLComp$Day, levels = c("1","4","10"))
###Figure 2A -----
library(ggpubr)


### SELECT ONLY WT AND DTDA COMPARISONS


Figure_2A <- LogFoldCounts_sum_dat_all_long_ALLComp %>% filter(Comp == "WT vs. dTDA") %>% ggballoonplot(
  x = "var", y="Day",
  size = "Count", fill = "Count", show.label = TRUE, 
  facet.by = c("Env", "Comp"), size.range = c(1, 17),
  ggtheme = guides(size = FALSE)
) +
  scale_fill_gradientn(colors = c("white","#fdb863","#e66101","#5e3c99"), values=c(0,0.4,1)) +
  scale_y_discrete(limits=rev) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) 


#Count total sum of unique genera increased by the WT and the dTDA respectively

WT_increased <- rbind(df_fig_SW_day1[,1:5], df_fig_SW_day4, df_fig_SW_day10[,1:5], 
      df_fig_Biofilm_day1[,1:5], df_fig_Biofilm_day4, df_fig_Biofilm_day10[,1:5]) %>% filter(var != "Control") %>% filter(LogFoldChange < 0)
length(unique(WT_increased$taxon_id))

dTDA_increased <- rbind(df_fig_SW_day1[,1:5], df_fig_SW_day4, df_fig_SW_day10[,1:5], 
                      df_fig_Biofilm_day1[,1:5], df_fig_Biofilm_day4, df_fig_Biofilm_day10[,1:5]) %>% filter(var != "Control") %>% filter(LogFoldChange > 0)

length(unique(dTDA_increased$taxon_id))
  

# dev.off()


#Taxonomic profile of Day 4 at order level (Figure 2B)

#Planktonic suspension Day 4 ANCOM-BC Order lvl. Comparing ref group set to WT ---- 
phylo_main_no_cont_NoOutliers_SW <- subset_samples(phylo_main_no_cont_NoOutliers, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="4")
phylo_main_no_cont_NoOutliers_SW_day4_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_SW_day4, Genus!="Phaeobacter_inhibens")



#SW
#Gloom/aggregate 


phylo_main_no_cont_NoOutliers_SW_Order_day4 <- aggregate_taxa(phylo_main_no_cont_NoOutliers_SW_day4_noPhaeo, "Order")
sample_data(phylo_main_no_cont_NoOutliers_SW_Order_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_Order_day4)$Treatment, levels = c("WT", "Control", "dTDA"))


#Run ANCOMBC differential analysis
out_phylo_main_no_cont_NoOutliers_SW_Order_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_Order_day4, formula = "Treatment",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.01, global = TRUE)



out_phylo_main_no_cont_NoOutliers_SW_Order_day4$res$diff_abn
#Results SW Day4
res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4 <- out_phylo_main_no_cont_NoOutliers_SW_Order_day4$res

#Add Relative abundances for color gradient in later plotting
#Transform to TSS in percentage
#Relative abundances 
phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
phylo_main_no_cont_NoOutliers_scalled_SW <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_SW, Day=="4")
phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, Genus!="Phaeobacter_inhibens")

phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, "Order")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order, function (x) {sum(x > 0) > 0}, prune=TRUE)
sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))

RelAbs_SW_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order))
RelAbs_SW_Order_day4$taxon_id <- rownames(RelAbs_SW_Order_day4)



colnames(RelAbs_SW_Order_day4)
RelAbs_SW_Order_day4_long <- reshape(RelAbs_SW_Order_day4, direction = "long",
                                     v.names = "Rel_Abs",
                                     varying = 1:ncol(RelAbs_SW_Order_day4)-1,
                                     times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
)

#Mean Rel_Abs per treatment per tax
RelAbs_SW_Order_day4_long_mean <- RelAbs_SW_Order_day4_long %>%
  group_by(taxon_id, time) %>%
  dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 



#dplyr::filter(RelAbs_SW_Order_day4_long_mean, grepl("Chloroplast",taxon_id))

#### Dataframe Day 4 for logFC
df_fig1_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$se * res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day4)[-1] = paste0(colnames(df_fig2_SW_day4)[-1], "SD")
colnames(df_fig1_SW_day4) <- sub("Treatment", "", colnames(df_fig1_SW_day4))
colnames(df_fig2_SW_day4) <- sub("Treatment", "", colnames(df_fig2_SW_day4))


df_fig_SW_day4 = df_fig1_SW_day4 %>% left_join(df_fig2_SW_day4, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

#df_fig_SW_day4 %>% filter(Control > 0) %>% mutate()

#reshape to lon format
df_fig_SW_day4 <- as.data.frame(df_fig_SW_day4)
df_fig_SW_day4 <- df_fig_SW_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day4 <- reshape(df_fig_SW_day4, direction='long', 
                          varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                          timevar='var',
                          times=c('Control', 'dTDA'),
                          v.names=c('LogFoldChange', 'SD'))

#Mutate rel abundances from mutant and Control treatment, then Relative abudances from log_obs_abn_adj_SW_Order_day4_mean
#Dataframe for LogFC for control enrichments (LogFC > 0) and rel abs
df_fig_SW_day4_Control_RelAbsMeans <- df_fig_SW_day4  %>% filter(var=='Control') %>% left_join(RelAbs_SW_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='Control')
#Dataframe for LogFC for mutant enrichments (LogFC > 0) and rel abs
df_fig_SW_day4_dTDA_RelAbsMeans <- df_fig_SW_day4  %>% filter(var=='dTDA') %>% left_join(RelAbs_SW_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='dTDA')
#Bind the two
df_fig_SW_day4_ControldTDA_RelAbsMeans <- rbind(df_fig_SW_day4_Control_RelAbsMeans, df_fig_SW_day4_dTDA_RelAbsMeans)
#Removing Rel Abs where LogFC is zero
Rel_Abs1_1 <-df_fig_SW_day4_ControldTDA_RelAbsMeans  %>% filter(LogFoldChange==0) %>% mutate(Rel_Abs_1 = LogFoldChange)#Add new column with zeros for RelAbs mutate()
Rel_Abs1_2 <- df_fig_SW_day4_ControldTDA_RelAbsMeans  %>% filter(LogFoldChange!=0) %>% mutate(Rel_Abs_1 = Rel_Abs) #Add new column with RelAbs mutate()
#Bind the two
df_fig_SW_day4_ControldTDA_RelAbsMeans1 <- rbind(Rel_Abs1_1, Rel_Abs1_2)


#Biofilm Day 4 ANCOM-BC Order lvl. Comparing ref group set to WT ---- 
phylo_main_no_cont_NoOutliers_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_Biofilm_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="4")
phylo_main_no_cont_NoOutliers_Biofilm_day4_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_Biofilm_day4, Genus!="Phaeobacter_inhibens")

phylo_main_no_cont_NoOutliers_Biofilm_Order_day4 <- aggregate_taxa(phylo_main_no_cont_NoOutliers_Biofilm_day4_noPhaeo, "Order")
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_Order_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_Order_day4)$Treatment, levels = c("WT", "Control", "dTDA"))


out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_Order_day4, formula = "Treatment",
                                                                p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                alpha = 0.01, global = TRUE)


#Results Biofilm Day4
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$res


#Relative abundances 
phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
sample_data(phylo_main_no_cont_NoOutliers_scalled)
phylo_main_no_cont_NoOutliers_scalled_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_Biofilm, Day=="4")
#Remove Phaeobacter
phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4, Genus!="Phaeobacter_inhibens")

phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4, "Order")
sample_data(phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))

RelAbs_Biofilm_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_Biofilm_day4_Order))
RelAbs_Biofilm_Order_day4$taxon_id <- rownames(RelAbs_Biofilm_Order_day4)

RelAbs_Biofilm_Order_day4_long <- reshape(RelAbs_Biofilm_Order_day4, direction = "long",
                                          v.names = "Rel_Abs",
                                          varying = 1:ncol(RelAbs_Biofilm_Order_day4)-1,
                                          times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
)

#Mean Rel_Abs per treatment per tax
RelAbs_Biofilm_Order_day4_long_mean <- RelAbs_Biofilm_Order_day4_long %>%
  group_by(taxon_id, time) %>%
  dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 


#### Dataframe Day 4 for logFC
df_fig1_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day4)[-1] = paste0(colnames(df_fig2_Biofilm_day4)[-1], "SD")
colnames(df_fig1_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day4))
colnames(df_fig2_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day4))


df_fig_Biofilm_day4 = df_fig1_Biofilm_day4 %>% left_join(df_fig2_Biofilm_day4, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_Biofilm_day4 %>% filter(Control > 0) %>% mutate()

#reshape to lon format
df_fig_Biofilm_day4 <- as.data.frame(df_fig_Biofilm_day4)
df_fig_Biofilm_day4 <- df_fig_Biofilm_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day4 <- reshape(df_fig_Biofilm_day4, direction='long', 
                               varying=c('Control','ControlSD','dTDA', 'dTDASD'), 
                               timevar='var',
                               times=c('Control', 'dTDA'),
                               v.names=c('LogFoldChange', 'SD'))

#Mutate rel abundances from mutant and Control treatment, then Relative abudances from log_obs_abn_adj_Biofilm_Order_day4_mean

#Dataframe for LogFC for control enrichments (LogFC > 0) and rel abs
df_fig_Biofilm_day4_Control_RelAbsMeans <- df_fig_Biofilm_day4  %>% filter(var=='Control') %>% left_join(RelAbs_Biofilm_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='Control')
#Dataframe for LogFC for mutant enrichments (LogFC > 0) and rel abs
df_fig_Biofilm_day4_dTDA_RelAbsMeans <- df_fig_Biofilm_day4  %>% filter(var=='dTDA') %>% left_join(RelAbs_SW_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='dTDA')
#Bind the two
df_fig_Biofilm_day4_ControldTDA_RelAbsMeans <- rbind(df_fig_Biofilm_day4_Control_RelAbsMeans, df_fig_Biofilm_day4_dTDA_RelAbsMeans)
#Removing Rel Abs where LogFC is zero
Rel_Abs1_1 <-df_fig_Biofilm_day4_ControldTDA_RelAbsMeans  %>% filter(LogFoldChange==0) %>% mutate(Rel_Abs_1 = LogFoldChange)#Replace rel_abundances for LogFC = 0 with zero in Order to have a better color scale in plot downstream
Rel_Abs1_2 <- df_fig_Biofilm_day4_ControldTDA_RelAbsMeans  %>% filter(LogFoldChange!=0) %>% mutate(Rel_Abs_1 = Rel_Abs) #Add new column with RelAbs
#Bind the two
df_fig_Biofilm_day4_ControldTDA_RelAbsMeans1 <- rbind(Rel_Abs1_1, Rel_Abs1_2)


##### Figure 2B Part 1----
df_fig_SW_day4_ControldTDA_RelAbsMeans1$Env <- rep("Planktonic suspension")
df_fig_Biofilm_day4_ControldTDA_RelAbsMeans1$Env <- rep("Biofilm")

df_fig_day4_ControldTDA_RelAbsMeans1_all <- rbind(df_fig_SW_day4_ControldTDA_RelAbsMeans1, df_fig_Biofilm_day4_ControldTDA_RelAbsMeans1)
#remove unOrderified
#df_fig_day4_ControldTDA_RelAbsMeans1_all<-df_fig_day4_ControldTDA_RelAbsMeans1_all[df_fig_day4_ControldTDA_RelAbsMeans1_all$taxon_id!="UnOrderified_UnOrderified_UnOrderified_UnOrderified" & df_fig_day4_ControldTDA_RelAbsMeans1_all$taxon_id!="Bacteria_Proteobacteria_UnOrderified",]

df_fig_day4_ControldTDA_RelAbsMeans1_all$time <- factor(df_fig_day4_ControldTDA_RelAbsMeans1_all$time)
levels(df_fig_day4_ControldTDA_RelAbsMeans1_all$time) <- c("WT vs. Control", "WT vs. dTDA") # Change names of treatment labels


#List according to rel abundance
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr <- df_fig_day4_ControldTDA_RelAbsMeans1_all %>% group_by(Env) %>% arrange(desc(Rel_Abs_1))
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr
#Dataframe rel abs for WT
WT_RelAbs_Biofilm_Order_day4_long_mean  <- RelAbs_Biofilm_Order_day4_long_mean %>% filter(time=='WT') 
WT_RelAbs_SW_Order_day4_long_mean  <- RelAbs_SW_Order_day4_long_mean %>% filter(time=='WT') 
WT_RelAbs_Biofilm_Order_day4_long_mean$Env <- rep("Biofilm")
WT_RelAbs_SW_Order_day4_long_mean$Env <- rep("Planktonic suspension")
#Bind the two envs
WT_RelAbs_Biofilm_SW_Order_day4_long_mean <- rbind(WT_RelAbs_Biofilm_Order_day4_long_mean, WT_RelAbs_SW_Order_day4_long_mean)
#Sort by drecreasing Order
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr %>% left_join(WT_RelAbs_Biofilm_Order_day4_long_mean, by = "taxon_id") 
#Join WT rel to dataframe
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr %>% left_join(WT_RelAbs_Biofilm_SW_Order_day4_long_mean, by = c("taxon_id", "Env"))  
#Make new colnames
colnames(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT) <- c("taxon_id","Treatment","LogFoldChange","SD","id","Comparison","Rel_Abs","Rel_Abs_print","Env","WT","Rel_Abs_WT")

#Clean up 
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT[,c(1,2,3,4,6,7,8,9,11)]
#Add new column to dataframe with relative abundances if LogFold < 0 add values from Rel_Abs_WT and else Rel_Abs
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel <- mutate(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT, Rel_Abs_ControldTDAWT = ifelse(LogFoldChange < 0, Rel_Abs_WT, Rel_Abs))
#Remove low abundant Orders
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel %>% filter(Rel_Abs_ControldTDAWT>=0.1)

#Round to one decimal
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt[,c(6,9)] <- round(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt[,c(6,9)], 1)
#WT- dTDA:
head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt  %>% filter (Comparison == "WT vs. dTDA") %>% filter (Env == "Biofilm"), 10)[,c(1,3,10)]
head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt  %>% filter (Comparison == "WT vs. dTDA") %>% filter (Env == "Planktonic suspension"), 10)[,c(1,3,10)]
#Counts
unique((df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt %>% filter (Comparison == "WT vs. Control") %>% filter (LogFoldChange != 0 ) )$taxon_id)
unique((df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt %>% filter (Comparison == "WT vs. dTDA") %>% filter (LogFoldChange != 0 ) )$taxon_id)
#WT- Control:
head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt  %>% filter (Comparison == "WT vs. Control") %>% filter (Env == "Biofilm"), 20) [,c(1,3,10)]
head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt  %>% filter (Comparison == "WT vs. Control") %>% filter (Env == "Planktonic suspension"), 10)[,c(1,3,10)]

#df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt_1 <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt[df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt$Rel_Abs_print!=0,]
df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt_1 <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt[df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt$Rel_Abs_ControldTDAWT!=0,]


#Planktonic suspension Day 4 ANCOM-BC Order lvl. Comparing ref group set to Control ---- 
phylo_main_no_cont_NoOutliers_SW <- subset_samples(phylo_main_no_cont_NoOutliers, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_SW, Day=="4")
phylo_main_no_cont_NoOutliers_SW_day4_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_SW_day4, Genus!="Phaeobacter_inhibens")

#phylo_main_no_cont_NoOutliers_SW_day4_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_SW_day4, Genus!="Phaeobacter_inhibens")
#Remove low abundant taxa
#phylo_main_no_cont_NoOutliers_SW_day4_noPhaeo <- 

#SW
#Gloom/aggregate 
library(microbiome)
library(ANCOMBC)

phylo_main_no_cont_NoOutliers_SW_Order_day4 <- aggregate_taxa(phylo_main_no_cont_NoOutliers_SW_day4_noPhaeo, "Order")
sample_data(phylo_main_no_cont_NoOutliers_SW_Order_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_SW_Order_day4)$Treatment, levels = c("Control", "WT", "dTDA"))


#Run ANCOMBC differential analysis
out_phylo_main_no_cont_NoOutliers_SW_Order_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_SW_Order_day4, formula = "Treatment",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.005, global = TRUE)



out_phylo_main_no_cont_NoOutliers_SW_Order_day4$res$diff_abn
#Results SW Day4
res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4 <- out_phylo_main_no_cont_NoOutliers_SW_Order_day4$res


#### Dataframe Day 4 for logFC
df_fig1_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$beta * res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_SW_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$se * res_out_phylo_main_no_cont_NoOutliers_SW_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_SW_day4)[-1] = paste0(colnames(df_fig2_SW_day4)[-1], "SD")
colnames(df_fig1_SW_day4) <- sub("Treatment", "", colnames(df_fig1_SW_day4))
colnames(df_fig2_SW_day4) <- sub("Treatment", "", colnames(df_fig2_SW_day4))


df_fig_SW_day4 = df_fig1_SW_day4 %>% left_join(df_fig2_SW_day4, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

#df_fig_SW_day4 %>% filter(Control > 0) %>% mutate()

#reshape to lon format
df_fig_SW_day4 <- as.data.frame(df_fig_SW_day4)
df_fig_SW_day4 <- df_fig_SW_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_SW_day4 <- reshape(df_fig_SW_day4, direction='long', 
                          varying=c('WT','WTSD','dTDA', 'dTDASD'), 
                          timevar='var',
                          times=c('WT', 'dTDA'),
                          v.names=c('LogFoldChange', 'SD'))

#Mutate rel abundances from mutant and Control treatment, then Relative abudances from log_obs_abn_adj_SW_Order_day4_mean
#Dataframe for LogFC for control enrichments (LogFC > 0) and rel abs
df_fig_SW_day4_WT_RelAbsMeans <- df_fig_SW_day4  %>% filter(var=='WT') %>% left_join(RelAbs_SW_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='WT')
#Dataframe for LogFC for mutant enrichments (LogFC > 0) and rel abs
df_fig_SW_day4_dTDA_RelAbsMeans <- df_fig_SW_day4  %>% filter(var=='dTDA') %>% left_join(RelAbs_SW_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='dTDA')
#Bind the two
df_fig_SW_day4_WTdTDA_RelAbsMeans <- rbind(df_fig_SW_day4_WT_RelAbsMeans, df_fig_SW_day4_dTDA_RelAbsMeans)
#Removing Rel Abs where LogFC is zero
Rel_Abs1_1 <-df_fig_SW_day4_WTdTDA_RelAbsMeans  %>% filter(LogFoldChange==0) %>% mutate(Rel_Abs_1 = LogFoldChange)#Add new column with zeros for RelAbs mutate()
Rel_Abs1_2 <- df_fig_SW_day4_WTdTDA_RelAbsMeans  %>% filter(LogFoldChange!=0) %>% mutate(Rel_Abs_1 = Rel_Abs) #Add new column with RelAbs mutate()
#Bind the two
df_fig_SW_day4_WTdTDA_RelAbsMeans1 <- rbind(Rel_Abs1_1, Rel_Abs1_2)


#Biofilm Day 4 ANCOM-BC Order lvl. Comparing ref group set to Control ---- 
phylo_main_no_cont_NoOutliers_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_Biofilm_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_Biofilm, Day=="4")
phylo_main_no_cont_NoOutliers_Biofilm_day4_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_Biofilm_day4, Genus!="Phaeobacter_inhibens")

phylo_main_no_cont_NoOutliers_Biofilm_Order_day4 <- aggregate_taxa(phylo_main_no_cont_NoOutliers_Biofilm_day4_noPhaeo, "Order")
sample_data(phylo_main_no_cont_NoOutliers_Biofilm_Order_day4)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_Biofilm_Order_day4)$Treatment, levels = c("Control", "dTDA", "WT"))


out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4 <- ancombc(phyloseq = phylo_main_no_cont_NoOutliers_Biofilm_Order_day4, formula = "Treatment",
                                                                p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                                group = "Treatment", struc_zero = TRUE, neg_lb = FALSE,
                                                                tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                                alpha = 0.01, global = TRUE)


#Results Biofilm Day4
res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4 <- out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$res



#### Dataframe Day 4 for logFC
df_fig1_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$beta * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_Biofilm_day4 = data.frame(res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$se * res_out_phylo_main_no_cont_NoOutliers_Biofilm_Order_day4$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_Biofilm_day4)[-1] = paste0(colnames(df_fig2_Biofilm_day4)[-1], "SD")
colnames(df_fig1_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig1_Biofilm_day4))
colnames(df_fig2_Biofilm_day4) <- sub("Treatment", "", colnames(df_fig2_Biofilm_day4))


df_fig_Biofilm_day4 = df_fig1_Biofilm_day4 %>% left_join(df_fig2_Biofilm_day4, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

#df_fig_Biofilm_day4 %>% filter(Control > 0) %>% mutate()

#reshape to lon format
df_fig_Biofilm_day4 <- as.data.frame(df_fig_Biofilm_day4)
df_fig_Biofilm_day4 <- df_fig_Biofilm_day4[,c(1,2,4,3,5)]
library(reshape)
df_fig_Biofilm_day4 <- reshape(df_fig_Biofilm_day4, direction='long', 
                               varying=c('dTDA','dTDASD','WT', 'WTSD'), 
                               timevar='var',
                               times=c('dTDA', 'WT'),
                               v.names=c('LogFoldChange', 'SD'))

#Mutate rel abundances from mutant and Control treatment, then Relative abudances from log_obs_abn_adj_Biofilm_Order_day4_mean
#Dataframe for LogFC for control enrichments (LogFC > 0) and rel abs
#Dataframe for LogFC for control enrichments (LogFC > 0) and rel abs
df_fig_Biofilm_day4_WT_RelAbsMeans <- df_fig_Biofilm_day4  %>% filter(var=='WT') %>% left_join(RelAbs_Biofilm_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='WT')
#Dataframe for LogFC for mutant enrichments (LogFC > 0) and rel abs
df_fig_Biofilm_day4_dTDA_RelAbsMeans <- df_fig_Biofilm_day4  %>% filter(var=='dTDA') %>% left_join(RelAbs_SW_Order_day4_long_mean, by = "taxon_id")  %>% filter(time=='dTDA')
#Bind the two
df_fig_Biofilm_day4_WTdTDA_RelAbsMeans <- rbind(df_fig_Biofilm_day4_WT_RelAbsMeans, df_fig_Biofilm_day4_dTDA_RelAbsMeans)
#Removing Rel Abs where LogFC is zero
Rel_Abs1_1 <-df_fig_Biofilm_day4_WTdTDA_RelAbsMeans  %>% filter(LogFoldChange==0) %>% mutate(Rel_Abs_1 = LogFoldChange)#Replace rel_abundances for LogFC = 0 with zero in order to have a better color scale in plot downstream
Rel_Abs1_2 <- df_fig_Biofilm_day4_WTdTDA_RelAbsMeans  %>% filter(LogFoldChange!=0) %>% mutate(Rel_Abs_1 = Rel_Abs) #Add new column with RelAbs
#Bind the two
df_fig_Biofilm_day4_WTdTDA_RelAbsMeans1 <- rbind(Rel_Abs1_1, Rel_Abs1_2)


##### Figure 2B part 2 ----
df_fig_SW_day4_WTdTDA_RelAbsMeans1$Env <- rep("Planktonic suspension")
df_fig_Biofilm_day4_WTdTDA_RelAbsMeans1$Env <- rep("Biofilm")

df_fig_day4_WTdTDA_RelAbsMeans1_all <- rbind(df_fig_SW_day4_WTdTDA_RelAbsMeans1, df_fig_Biofilm_day4_WTdTDA_RelAbsMeans1)
#remove unclassified
#df_fig_day4_WTdTDA_RelAbsMeans1_all<-df_fig_day4_WTdTDA_RelAbsMeans1_all[df_fig_day4_WTdTDA_RelAbsMeans1_all$taxon_id!="UnOrderified_UnOrderified_UnOrderified_UnOrderified" & df_fig_day4_WTdTDA_RelAbsMeans1_all$taxon_id!="Bacteria_Proteobacteria_UnOrderified",]

df_fig_day4_WTdTDA_RelAbsMeans1_all$time <- factor(df_fig_day4_WTdTDA_RelAbsMeans1_all$time)
levels(df_fig_day4_WTdTDA_RelAbsMeans1_all$time) <- c("Control vs. dTDA", "Control vs. WT") # Change names of treatment labels

df_fig_day4_WTdTDA_RelAbsMeans1_all_decr <- df_fig_day4_WTdTDA_RelAbsMeans1_all %>% group_by(Env) %>% arrange(desc(Rel_Abs_1))

#Dataframe rel abs for WT
dTDA_RelAbs_Biofilm_Order_day4_long_mean  <- RelAbs_Biofilm_Order_day4_long_mean %>% filter(time=='dTDA') 
dTDA_RelAbs_SW_Order_day4_long_mean  <- RelAbs_SW_Order_day4_long_mean %>% filter(time=='dTDA') 
dTDA_RelAbs_Biofilm_Order_day4_long_mean$Env <- rep("Biofilm")
dTDA_RelAbs_SW_Order_day4_long_mean$Env <- rep("Planktonic suspension")
#Bind the two envs
dTDA_RelAbs_Biofilm_SW_Order_day4_long_mean <- rbind(dTDA_RelAbs_Biofilm_Order_day4_long_mean, dTDA_RelAbs_SW_Order_day4_long_mean)
#Sort by drecreasing order
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr %>% left_join(dTDA_RelAbs_Biofilm_Order_day4_long_mean, by = "taxon_id") 
#Join WT rel to dataframe
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA <- df_fig_day4_WTdTDA_RelAbsMeans1_all_decr %>% left_join(dTDA_RelAbs_Biofilm_SW_Order_day4_long_mean, by = c("taxon_id", "Env"))  
#Make new colnames
colnames(df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA) <- c("taxon_id","Treatment","LogFoldChange","SD","id","Comparison","Rel_Abs","Rel_Abs_print","Env","WT","Rel_Abs_dTDA")

#Clean up 
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA <- df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA[,c(1,2,3,4,6,7,8,9,11)]
#Add new column to dataframe with relative abundances if LogFold < 0 add values from Rel_Abs_WT and else Rel_Abs
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel <- mutate(df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA, Rel_Abs_ControldTDAWT = ifelse(LogFoldChange < 0, Rel_Abs_dTDA, Rel_Abs))
#Remove low abundant orders
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt <- df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel %>% filter(Rel_Abs_ControldTDAWT>=0.1)

#head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT , 20)
# #Remove tax for which the sum of WT and Control/dTDA rel abs is below 0.05 mean rel abs
# df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel %>% rowwise() %>% mutate(sum = sum(c(Rel_Abs, Rel_Abs_WT))) %>% filter(sum>=0.2)
# #manually remove Opitutales - don't know how else to do
# df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt  %>% filter (taxon_id != "Opitutales")
#Round to one decimal
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt[,c(6,9)] <- round(df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt[,c(6,9)], 1)
#Control vs. dTDA
head(df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt  %>% filter (Comparison == "dTDA - Control") %>% filter (Env == "Biofilm"), 10)[,c(1,3,10)]
head(df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt  %>% filter (Comparison == "dTDA - Control") %>% filter (Env == "Planktonic suspension"), 10)[,c(1,3,10)]
# #Counts
# unique((df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt %>% filter (Comparison == "WT vs. Control") %>% filter (LogFoldChange != 0 ) )$taxon_id)
# unique((df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt %>% filter (Comparison == "WT vs. dTDA") %>% filter (LogFoldChange != 0 ) )$taxon_id)
# #WT- Control:
# head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt  %>% filter (Comparison == "WT vs. Control") %>% filter (Env == "Biofilm"), 20) [,c(1,3,10)]
# head(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt  %>% filter (Comparison == "WT vs. Control") %>% filter (Env == "Planktonic suspension"), 10)[,c(1,3,10)]

#df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt_1 <- df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt[df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_filt$Rel_Abs_print!=0,]
df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt <- df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt[df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt$Rel_Abs_ControldTDAWT!=0,]


df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt_ControldTDA <- df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt %>% filter(Comparison =="Control vs. dTDA")

df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt_ControldTDA <- df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt %>% filter(Comparison =="Control vs. dTDA")


#Figure 2B part 1 + 2 -----
df_fig_day4_ANCOMBC_Final <- rbind(df_fig_day4_ControldTDA_RelAbsMeans1_all_decr_WT_rel_filt, df_fig_day4_WTdTDA_RelAbsMeans1_all_decr_dTDA_rel_filt_ControldTDA)

#Remove "Bacteria_", "Bacteria_Proteobacteria", "Unclassified_" from taxon_id

df_fig_day4_ANCOMBC_Final$taxon_id <- str_replace_all(df_fig_day4_ANCOMBC_Final$taxon_id, c("Bacteria_Proteobacteria_Unclassified_Unclassified"="Proteobacteria_Unclassified"))
df_fig_day4_ANCOMBC_Final$taxon_id <- str_replace_all(df_fig_day4_ANCOMBC_Final$taxon_id, c("Bacteria_Proteobacteria_" = "", "Bacteria_Planctomycetes_" = "", "Unclassified_"="", "Bacteria_"=""))
df_fig_day4_ANCOMBC_Final$taxon_id <- str_replace_all(df_fig_day4_ANCOMBC_Final$taxon_id, c("Unclassified"="Uncl."))

###ONLY SELECT FOR THE WT vs dTDA comparison

Figure_2B <- df_fig_day4_ANCOMBC_Final %>% filter(Comparison == "WT vs. dTDA")  %>% filter(LogFoldChange!=0) %>%
  ggplot(aes(x=LogFoldChange, y=reorder(taxon_id, LogFoldChange), label=Rel_Abs_ControldTDAWT)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (Rel_Abs_ControldTDAWT)),
           position = position_dodge(width = 0.7), color="black") +
  geom_errorbar(aes(xmin = LogFoldChange - SD, xmax = LogFoldChange + SD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  labs(y = NULL, x="Log fold change",
       title = "") + scale_fill_gradientn(name="Mean relative\nabundance (%)",colours=c("white","#fdb863","#e66101","#5e3c99") , values=c(0,0.1,1)) +#, breaks=c(0,0.2,3,6,9)) +
  #geom_text(aes(y = taxon_id, aes(label=ifelse(Rel_Abs==0, value,""))))
  #geom_text(position=position_dodge(0.9), aes(label=ifelse(Rel_Abs==0, "\u2666" ,"")), size = 5) +
  #geom_text(aes(label=ifelse(Rel_Abs==0, "*", "" ), ifelse(sign(LogFoldChange)==-1,1,-1)), size = 4) +
  facet_grid(vars(Env), vars(Comparison)) + 
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()) 

# dev.off()

#FINAL Figure 2: A and 2B merged in one plot ----


Figure_2B_1 <- plot_grid(NULL,Figure_2B + theme(legend.position = "bottom", text = element_text(size=10)), 
          rel_widths = c(1,1.6))



pdf(file = "Figures/Figure_3.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 7) # The height of the plot in inches

ggdraw(Figure_2B_1) + draw_label("B", x = 0.5, y = 0.97, fontface = "bold") + draw_label("A", x = 0.025, y = 0.97, fontface = "bold") +
  draw_plot(Figure_2A + theme(legend.position = "bottom", text = element_text(size=10)), 
            x = 0.45, y = 0.97, hjust = 1, vjust = 1, width = 0.35, height = 0.55)


dev.off()


tiff("Figures/Figure_3.tiff", units="in", width=6, height=6, res=300)

ggdraw(Figure_2B_1) + draw_label("B", x = 0.41, y = 0.97, fontface = "bold") + draw_label("A", x = 0.025, y = 0.97, fontface = "bold") +
  draw_plot(Figure_2A + theme(legend.position = "bottom", text = element_text(size=10)), 
            x = 0.35, y = 0.96, hjust = 1, vjust = 1, width = 0.35, height = 0.955)


dev.off()

 ####### venn ########
WT_biofilm_facilitated <- df_fig_Biofilm_day4 %>% filter (group == "WT" , var == "dTDA")  %>% select(taxon_id)
WT_biofilm_reduced <- df_fig_Biofilm_day4 %>% filter (var == "dTDA", group == "Control & dTDA") %>% select(taxon_id)


WT_biofilm_facilitated[[1]]

df_fig_SW_day4 %>% filter (group == "WT")
WT_SW_facilitated <- df_fig_SW_day4 %>% filter (group == "WT" , var == "dTDA")  %>% select(taxon_id)
WT_SW_reduced <- df_fig_SW_day4 %>% filter  (var == "dTDA", group == "Control & dTDA") %>% select(taxon_id)


library(venn)
venn_WT_Delta_Biofilm_vs_SW_increased <- venn(list(WT_biofilm_facilitated[[1]], WT_SW_facilitated[[1]]), snames = "Biofilm, Planktonic Suspension")
#63 + 55 + 54 = 172, 52% 32%  32% 
venn_WT_Delta_Biofilm_vs_SW_reduced <- venn(list(WT_biofilm_reduced[[1]], WT_SW_reduced[[1]]), snames = "Biofilm, Planktonic Suspension")
#39 + 20 + 29 = 88, 44% 22% 33% #Overvejende større andel som påvirkes i biofilm end planktoniske suspension 

# #Shared vectors of same Phyla 
# venn_WT_Delta_Biofilm <- venn(list(WT_tax_count$Phylum, Delta_tax_count$Phylum), snames = "WT_vectors, Delta_vectors")
# attr(x = venn_WT_Delta_Biofilm, "intersections")$WT_vectors
# #Shared vectors of same Class
# venn_WT_Delta_Biofilm <- venn(list(WT_tax_count$Class, Delta_tax_count$Class), snames = "WT_vectors, Delta_vectors")
# attr(x = venn_WT_Delta_Biofilm, "intersections")$WT_vectors


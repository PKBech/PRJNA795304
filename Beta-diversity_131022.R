setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/phyloseq/")
library(phyloseq)
library(ggplot2)
library(vegan)
library(microbiome)
library(gridExtra)

#Load in the phyloseq object 
load('phylo_main_no_cont_NoOutliers.RData')

phylo_main_no_cont_NoOutliers_Genus <- aggregate_taxa(phylo_main_no_cont_NoOutliers, "Genus")

#Filtering all taxa wiuth read sum < 100
phylo_main_no_cont_NoOutliers_Genus_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_Genus, function (x) {sum(x > 100) > 0}, prune=TRUE)

#Phaeobacter OTU removed
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_Genus_filter, Genus != "Phaeobacter_inhibens")

#Scalling to total sum TSS
phylo_main_no_cont_NoOutliers_Genus_filter_norm = transform_sample_counts(phylo_main_no_cont_NoOutliers_Genus_filter, function(x) 100000 * x/sum(x))
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo_norm = transform_sample_counts(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, function(x) 100000 * x/sum(x))
#Rarefying instead in order to make sure that scaling is neutral to the rarefying method
phylo_main_no_cont_NoOutliers_Genus_filter_norm = rarefy_even_depth(phylo_main_no_cont_NoOutliers_Genus_filter)
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo_norm = rarefy_even_depth(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo)


#Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions
genus_norm_dat <- as.data.frame((as(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter_norm), "matrix")))
genus_norm_dat <- round(genus_norm_dat)
genus_norm_dat_t <- t(genus_norm_dat)

noPhaeo_genus_norm_dat <- as.data.frame((as(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo_norm), "matrix")))
noPhaeo_genus_norm_dat <- round(noPhaeo_genus_norm_dat)
noPhaeo_genus_norm_dat_t <- t(noPhaeo_genus_norm_dat)

#Squar-root transform matrix
sqrt_genus_norm_dat_t= sqrt(genus_norm_dat_t)
sqrt_noPhaeo_genus_norm_dat_t = sqrt(noPhaeo_genus_norm_dat_t)

#Bray-curtis distance matrix
Genus_dist_bray = as.matrix((vegdist(sqrt_genus_norm_dat_t, "bray")))
noPhaeo_Genus_dist_bray = as.matrix((vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray")))

#Metadata dataframe
genus_meta_dat <- data.frame(sample_data(phylo_main_no_cont_NoOutliers_Genus_filter))

#Beta-disper ----
## Calculate multivariate dispersions from all 18 combinations
bdisp_treatment_day_env <- betadisper(vegdist(sqrt_genus_norm_dat_t, "bray"), 
                              factor(paste(genus_meta_dat$Treatment,genus_meta_dat$Day, genus_meta_dat$Environment, sep=" ")))


## Perform test
anova(bdisp_treatment_day_env)

## Permutation test for F
permutest(bdisp_treatment_day_env, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(bdisp.HSD <- TukeyHSD(bdisp_treatment_day_env))
plot(bdisp.HSD)
dim(bdisp.HSD$group)

# 
# #Remove "Seawater at Day 0" for down-stream analyses for betadisper and final adonis
phylo_main_no_cont_NoOutliers_Genus_filter <- subset_samples(phylo_main_no_cont_NoOutliers_Genus_filter, Environment != "Seawater at Day 0")
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo <- subset_samples(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, Environment != "Seawater at Day 0")


#Scalling to total sum TSS
phylo_main_no_cont_NoOutliers_Genus_filter_norm = transform_sample_counts(phylo_main_no_cont_NoOutliers_Genus_filter, function(x) 100000 * x/sum(x))
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo_norm = transform_sample_counts(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, function(x) 100000 * x/sum(x))


# #Rarefying instead
# phylo_main_no_cont_NoOutliers_Genus_filter_norm = rarefy_even_depth(phylo_main_no_cont_NoOutliers_Genus_filter, rngseed = 100)
# #7 ASVs are removed
# phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo_norm = rarefy_even_depth(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo, rngseed = 100)
# #23 ASVs are removed

#Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions
genus_norm_dat <- as.data.frame((as(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter_norm), "matrix")))
genus_norm_dat <- round(genus_norm_dat)
genus_norm_dat_t <- t(genus_norm_dat)


genus_norm_dat_noNORM <- as.data.frame((as(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo), "matrix")))

noPhaeo_genus_norm_dat <- as.data.frame((as(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo_norm), "matrix")))
noPhaeo_genus_norm_dat <- round(noPhaeo_genus_norm_dat)
noPhaeo_genus_norm_dat_t <- t(noPhaeo_genus_norm_dat)

genus_norm_dat_noNORM_t <- t(genus_norm_dat_noNORM)

#Squar-root transform matrix
sqrt_genus_norm_dat_t= sqrt(genus_norm_dat_t)
sqrt_noPhaeo_genus_norm_dat_t = sqrt(noPhaeo_genus_norm_dat_t)

#Bray-curtis distance matrix
Genus_dist_bray = as.matrix((vegdist(sqrt_genus_norm_dat_t, "bray")))
noPhaeo_Genus_dist_bray = as.matrix((vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray")))

genus_norm_dat_noNORM_t_dist_bray = as.matrix((vegdist(genus_norm_dat_noNORM_t, "bray")))
noPhaeo_genus_norm_dat_t_noSqrt_dist_bray = as.matrix((vegdist(noPhaeo_genus_norm_dat_t, "bray")))



#Metadata dataframe
genus_meta_dat <- data.frame(sample_data(phylo_main_no_cont_NoOutliers_Genus_filter_norm))

bdisp_time <- betadisper(vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray"), genus_meta_dat$Day)
bdisp_treatment <- betadisper(vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray"), genus_meta_dat$Treatment)
bdisp_environment <- betadisper(vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray"), genus_meta_dat$Environment)

genus_meta_dat_vec <- as.factor((paste(genus_meta_dat$Environment, genus_meta_dat$Treatment, genus_meta_dat$Day, sep = "_")))
bdisp_all <- betadisper(vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray"), genus_meta_dat_vec)


## Perform test
anova(bdisp_time)
anova(bdisp_treatment)
anova(bdisp_environment)

anova(bdisp_all)
## Permutation test for F
permutest(bdisp_all, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD <- TukeyHSD(bdisp_all))
plot(bdisp.HSD)




#PERMANOVA (Adonis) ----
adonis_Genus <- adonis(Genus_dist_bray ~ Treatment * Environment * Day, genus_meta_dat)
adonis_Genus_dat <- as.data.frame(adonis_Genus$aov.tab)
noPhao_adonis <- adonis(noPhaeo_Genus_dist_bray ~ Treatment * Environment * Day, genus_meta_dat)
noPhao_adonis_Genus_dat <- as.data.frame(noPhao_adonis$aov.tab)


#adonis(genus_norm_dat_noNORM_t_dist_bray ~ Treatment * Environment * Day, genus_meta_dat)
#adonis(noPhaeo_genus_norm_dat_t_noSqrt_dist_bray ~ Treatment * Environment * Day, genus_meta_dat)



write.table(adonis_Genus_dat, file = "adonis_Genus_dat.csv", sep = "\t", dec = ",")
write.table(noPhao_adonis_Genus_dat, file = "noPhao_adonis_Genus_dat.csv", sep = "\t", dec = ",")


#Pairwise PERMANOVA (Adonis2) ----
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairwise.adonis2_noPhaeo_Genus <- pairwise.adonis2(noPhaeo_Genus_dist_bray ~ Treatment * Day *  Environment, genus_meta_dat)

pairwise.adonis2(Genus_dist_bray ~ Treatment * Day *  Environment, genus_meta_dat)$dTDA_vs_WT
pairwise.adonis2_noPhaeo_Genus$dTDA_vs_WT

pairwise.adonis2(genus_norm_dat_noNORM_t_dist_bray ~ Treatment * Day *  Environment, genus_meta_dat)$dTDA_vs_WT


pairwise.adonis2_noPhaeo_Genus_WT_dTDA_dat <- pairwise.adonis2_noPhaeo_Genus
write.table(pairwise.adonis2_noPhaeo_Genus_WT_dTDA_dat, file = "pairwise.adonis2_noPhaeo_Genus_WT_dTDA_dat.csv", sep = "\t", dec = ",")



# #Subset for per day
# # day 1
# genus_meta_dat_day1 <- genus_meta_dat[genus_meta_dat$Day %in% c('1'),]
# Genus_dist_bray_day1 <- Genus_dist_bray[rownames(Genus_dist_bray) %in% rownames(genus_meta_dat_day1),rownames(Genus_dist_bray) %in% rownames(genus_meta_dat_day1)]
# pairwise.adonis2(Genus_dist_bray_day1 ~ Treatment * Environment, genus_meta_dat_day1)
# 
# # day 4
# genus_meta_dat_day4 <- genus_meta_dat[genus_meta_dat$Day %in% c('4'),]
# Genus_dist_bray_day4 <- Genus_dist_bray[rownames(Genus_dist_bray) %in% rownames(genus_meta_dat_day4),rownames(Genus_dist_bray) %in% rownames(genus_meta_dat_day4)]
# pairwise.adonis2(Genus_dist_bray_day4 ~ Treatment * Environment, genus_meta_dat_day4)
# 
# # day 10
# genus_meta_dat_day10 <- genus_meta_dat[genus_meta_dat$Day %in% c('10'),]
# Genus_dist_bray_day10 <- Genus_dist_bray[rownames(Genus_dist_bray) %in% rownames(genus_meta_dat_day10),rownames(Genus_dist_bray) %in% rownames(genus_meta_dat_day10)]
# pairwise.adonis2(Genus_dist_bray_day10 ~ Treatment * Environment, genus_meta_dat_day10)

#Include seawater at day 0 again for NMDS


# NDMS -----
library(goeveg)
dimcheckMDS(noPhaeo_Genus_dist_bray, k = 10) # Stress values < 0.1 with k=3
NMDS_ord_bray = metaMDS(noPhaeo_Genus_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(sqrt_genus_norm_dat_t)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species

#NMDS time------------------
NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS3 = NMDS3, Day = genus_meta_dat$Day)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))


NMDS_plot_time_NMDS1_2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, colour = Day, fill = Day)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_time,
    aesthetics = c("colour", "fill") 
  )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


NMDS_plot_time_NMDS1_2

NMDS_plot_time_NMDS2_3 <- ggplot(NMDS, aes(x=NMDS2, y=NMDS3, colour = Day, fill = Day)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_time,
    aesthetics = c("colour", "fill") 
  ) + xlim(-0.38,0.38) + ylim(-0.4,0.4) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


NMDS_plot_time_NMDS2_3

NMDS_plot_time_NMDS1_3 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS1, colour = Day, fill = Day)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_time,
    aesthetics = c("colour", "fill") 
  ) + ylim(-0.58,0.58) + xlim(-0.4,0.4) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


legend_time <- cowplot::get_legend(NMDS_plot_time_NMDS1_2 + theme(legend.position="bottom") )
str(legend_time)

#library(gridExtra)
#library(ggpubr)
library(gridExtra)

NMDS_plot_time <- grid.arrange(NMDS_plot_time_NMDS1_2, NMDS_plot_time_NMDS2_3, NMDS_plot_time_NMDS1_3, legend_time, ncol=3, nrow = 2, 
                               layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                               widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))

#NMDS environment-----------------

NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS3 = NMDS3, Environment = genus_meta_dat$Environment)
NMDS$Environment <- factor(NMDS$Environment, levels = c("Seawater at Day 0", "Biofilm", "Planktonic suspension"))
cols_env <- c("Seawater at Day 0" = "#f7f7f7", "Biofilm" = "#f1a340", "Planktonic suspension" = "#998ec3")

# NMDS$Environment <- factor(NMDS$Environment, levels = c( "Biofilm", "Planktonic suspension"))
# cols_env <- c("Biofilm" = "#f1a340", "Planktonic suspension" = "#998ec3")



NMDS_plot_env1_2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, colour = Environment, fill = Environment)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_env,
    aesthetics = c("colour", "fill") 
  ) + xlim(-0.58,0.65) + ylim(-0.4,0.4) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


NMDS_plot_env1_2

NMDS_plot_env2_3 <- ggplot(NMDS, aes(x=NMDS2, y=NMDS3, colour = Environment, fill = Environment)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_env,
    aesthetics = c("colour", "fill") 
  )  + xlim(-0.35,0.3) + ylim(-0.28,0.28) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

NMDS_plot_env2_3

NMDS_plot_env1_3 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS1, colour = Environment, fill = Environment)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_env,
    aesthetics = c("colour", "fill") 
  ) +  ylim(-0.58,0.65) + xlim(-0.22,0.22) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

NMDS_plot_env1_3


legend_env <- cowplot::get_legend(NMDS_plot_env1_2 + theme(legend.position="bottom") )


#library(gridExtra)
#library(ggpubr)

NMDS_plot_env <- grid.arrange(NMDS_plot_env1_2, NMDS_plot_env2_3, NMDS_plot_env1_3, legend_env, ncol=3, nrow = 2, 
                              layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                              widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))

#NMDS Treatmens----------------

NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS3 = NMDS3, Treatment = genus_meta_dat$Treatment)
NMDS$Treatment <- factor(NMDS$Treatment, levels = c("Seawater at Day 0","WT", "dTDA", "Control"))
cols_treat <- c("Seawater at Day 0" = "#f7f7f7" ,"WT" = "#a6611a", "dTDA" = "#dfc27d", "Control" = "#80cdc1")
# NMDS$Treatment <- factor(NMDS$Treatment, levels = c("WT", "dTDA", "Control"))
# cols_treat <- c("WT" = "#a6611a", "dTDA" = "#dfc27d", "Control" = "#80cdc1")


NMDS_plot_sys1_2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, colour = Treatment, fill = Treatment)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_treat,
    aesthetics = c("colour", "fill") 
  ) + xlim(-0.65,0.65) + ylim(-0.4,0.4) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

NMDS_plot_sys1_2


NMDS_plot_sys2_3 <- ggplot(NMDS, aes(x=NMDS2, y=NMDS3, colour = Treatment, fill = Treatment)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_treat,
    aesthetics = c("colour", "fill") 
  ) +  xlim(-0.32,0.32) + ylim(-0.25,0.25)+ theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

NMDS_plot_sys2_3

NMDS_plot_sys1_3 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS1, colour = Treatment, fill = Treatment)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) +
  scale_colour_manual(
    values = cols_treat,
    aesthetics = c("colour", "fill") 
  ) + ylim(-0.7,0.65) + xlim(-0.25,0.25) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

NMDS_plot_sys1_3


legend_sys <- cowplot::get_legend(NMDS_plot_sys1_2 + theme(legend.position="bottom") )


library(gridExtra)
library(ggpubr)

NMDS_plot_sys <- grid.arrange(NMDS_plot_sys1_2, NMDS_plot_sys2_3, NMDS_plot_sys1_3, legend_sys, ncol=3, nrow = 2, 
                              layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                              widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))


#Figure 1 ----
library(cowplot)
pdf(file = "Figures/Figure_1_NMDS_noPhaeo_rarefaction.pdf",   # The directory you want to save the file in
    width = 8.3, # The width of the plot in inches
    height = 9.7) # The height of the plot in inches

plot_grid(NMDS_plot_time, NMDS_plot_sys, NMDS_plot_env, ncol = 1, labels = c('A', 'B', 'C'))

dev.off()

tiff("Figures/Figure_1_NMDS_noPhaeo_rarefaction.tiff", units="in", width=8.5, height=9, res=300)

plot_grid(NMDS_plot_time, NMDS_plot_sys, NMDS_plot_env, ncol = 1, labels = c('A', 'B', 'C'))

dev.off()


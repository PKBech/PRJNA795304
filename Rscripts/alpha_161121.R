library(vegan)
library(dplyr)
library(phyloseq)
library(emmeans)
library(lme4)
library(DescTools)
library(cowplot)
library(ggh4x)
library(microbiome)
library(ggbeeswarm)

# library(vipor)

#Alpha diversity for Biofilm Genus level ------
setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/phyloseq/")

#Load in the phyloseq object 
load('phylo_main_no_cont_NoOutliers.RData')

phylo_main_no_cont_NoOutliers_Genus <- aggregate_taxa(phylo_main_no_cont_NoOutliers, "Genus")

#Filtering all taxa wiuth read sum < 100
phylo_main_no_cont_NoOutliers_Genus_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_Genus, function (x) {sum(x > 100) > 0}, prune=TRUE)

#Normalize to total sum
phylo_main_no_cont_NoOutliers_Genus_filter_norm = transform_sample_counts(phylo_main_no_cont_NoOutliers_Genus_filter, function(x) 100000 * x/sum(x))
#phylo_main_no_cont_NoOutliers_Genus_filter_rar <- rarefy_even_depth(phylo_main_no_cont_NoOutliers_Genus_filter, sample.size = min(sample_sums(phylo_main_no_cont_NoOutliers_Genus_filter)),
#                  rngseed = 100, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#Observed genera
Observed_Genera <- data.frame(t(estimateR(t(round(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter_norm)) ))))
Observed_Genera_tab <- data.frame(cbind(Obs = Observed_Genera$S.obs, Treatment=sample_data(phylo_main_no_cont_NoOutliers_Genus_filter_norm)$Treatment, Day=sample_data(phylo_main_no_cont_NoOutliers_Genus_filter_norm)$Day, Environment=sample_data(phylo_main_no_cont_NoOutliers_Genus_filter_norm)$Environment), Rep_tec=sample_data(phylo_main_no_cont_NoOutliers_Genus_filter_norm)$Group_ind_metadata)
Observed_Genera_tab$Obs <- as.numeric(Observed_Genera_tab$Obs)
#Observed_Genera_tab$Simpson_Genera <- as.numeric(Observed_Genera_tab$Simpson_Genera)

#LMer by emmeans ------
## Subset for Biofilm samples
Observed_Genera_tab_biofilm <- Observed_Genera_tab %>% filter(Environment=="Biofilm")

#Add replicate ID to dataframe
Observed_Genera_tab_biofilm$ID=gsub("(.*)(T[0-9])(.*)","\\1",Observed_Genera_tab_biofilm$Rep_tec)

#LMER
Observed_Genera_tab_biofilm_LMER=lmer(Obs~Treatment*Day + (1|ID), Observed_Genera_tab_biofilm)
anova(Observed_Genera_tab_biofilm_LMER)

Observed_Genera_tab_biofilm_LMER_EM=emmeans(Observed_Genera_tab_biofilm_LMER, ~Treatment | Day)
Observed_Genera_tab_biofilm_LMER_EM_stat <- contrast(Observed_Genera_tab_biofilm_LMER_EM, interaction = "pairwise", adjust = "Bonferroni")

## Subset for Planktonic suspension samples
Observed_Genera_tab_SW <- Observed_Genera_tab %>% filter(Environment=="Planktonic suspension")

#Add replicate ID to dataframe
Observed_Genera_tab_SW$ID=gsub("(.*)(T[0-9])(.*)","\\1",Observed_Genera_tab_SW$Rep_tec)

#LMER
Observed_Genera_tab_SW_LMER=lmer(Obs~Treatment*Day + (1|ID), Observed_Genera_tab_SW)
anova(Observed_Genera_tab_SW_LMER)

Observed_Genera_tab_SW_LMER_EM=emmeans(Observed_Genera_tab_SW_LMER, ~Treatment | Day)
Observed_Genera_tab_SW_LMER_EM_stat <- contrast(Observed_Genera_tab_SW_LMER_EM, interaction = "pairwise", adjust = "Bonferroni")

#Make  Dunnest test with seawater at day 0 as control 

#Make all combinations to individual groups
groups <- factor(paste(Observed_Genera_tab$Treatment,Observed_Genera_tab$Day, Observed_Genera_tab$Environment))

## Dunnett's test ----
DunnettTes_Day0 <- DunnettTest(Observed_Genera_tab$Obs, groups,control = "Seawater at Day 0 0 Seawater at Day 0")
DunnettTes_Day0_dat <- as.data.frame(DunnettTes_Day0$`Seawater at Day 0 0 Seawater at Day 0`)

write.table(DunnettTes_Day0_dat, file = "DunnettTes_Day0_dat.csv", dec=",", sep="_")


#Preparing plot ------
#Clean up dataframe for plotting
#Replace "Seawater at day 0" with "Planktonic suspension"
Observed_Genera_tab$Environment <-  ifelse(Observed_Genera_tab$Environment=="Seawater at Day 0", "Planktonic suspension", Observed_Genera_tab$Environment)
# and reorder days, treatment and environment
Observed_Genera_tab$Treatment <- factor(Observed_Genera_tab$Treatment, levels = c("Seawater at Day 0","Control","WT", "dTDA"))
Observed_Genera_tab$Environment <- factor(Observed_Genera_tab$Environment, levels = c("Biofilm", "Planktonic suspension"))
Observed_Genera_tab$Day <- factor(Observed_Genera_tab$Day, levels = c("0","1", "4", "10"))


#########Means and standard deviations#########
Observed_Genera_tab_stat <- Observed_Genera_tab %>% group_by(Environment, Treatment, Day) %>%
  dplyr::summarise(n=n(),'Observed richness' = mean(Obs), sd = sd(Obs), max=max(Obs))


#Define colors
colors <- c("Seawater at Day 0" ="grey", "Control" = "#80cdc1", "WT" = "#a6611a", "dTDA" = "#dfc27d")

## Beeswarm Plot -------
p <- ggplot(Observed_Genera_tab, aes(x = Day, y = Obs, col = Treatment)) + 
  geom_quasirandom(aes(col = Treatment), dodge.width=.8, cex=2)  +
  stat_summary(aes(group = Treatment), fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2, position=position_dodge(0.8)) +
  stat_summary(aes(group = Treatment),fun=mean, geom="point", color="black", position=position_dodge(0.8)) +
  facet_grid(cols = vars(Environment), scales="free_x") + 
  xlab("Day") + ylab("Observed richness") + #guides(fill = FALSE) +
  scale_colour_manual(
    values = colors,
    aesthetics = c("colour") 
  ) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


legend_p_plot <- cowplot::get_legend(p + theme(legend.position="bottom",  text = element_text(size=10)))



#Add emmeans significance (LETTERS) to bars
Observed_Genera_tab_biofilm_LMER_EM_stat_dat <-as.data.frame(summary(Observed_Genera_tab_biofilm_LMER_EM_stat)[c("Treatment_pairwise","p.value")])
Observed_Genera_tab_SW_LMER_EM_stat_dat <-as.data.frame(summary(Observed_Genera_tab_SW_LMER_EM_stat)[c("Treatment_pairwise","p.value")])
Observed_Genera_tab_biofilm_LMER_EM_stat_dat$Environment <- rep("Biofilm")
Observed_Genera_tab_SW_LMER_EM_stat_dat$Environment <- rep("Planktonic suspension")
emmeans_dat_all <- rbind(Observed_Genera_tab_biofilm_LMER_EM_stat_dat, Observed_Genera_tab_SW_LMER_EM_stat_dat)
#Manual add significans letters to dataframe
emmeans_dat_all$signif <- c("a", "b", "b",
                            "a", "a", "a",
                            "a", "a", "a",
                            "a", "b", "b",
                            "a", "a", "a",
                            "a", "a", "a")

emmeans_dat_all$Treatment <- c("Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA")

emmeans_dat_all$Day <- c( rep("1", 3), rep("4", 3), rep("10", 3), rep("1", 3), rep("4", 3), rep("10", 3))
#Clean up table
emmeans_dat_all <- emmeans_dat_all[,3:6]
#Add timepoint zero back to the dataframe
emmeans_dat_all[19,]<-c("Planktonic suspension", "", "Seawater at Day 0", "0")
emmeans_dat_all$Treatment <- factor(emmeans_dat_all$Treatment, levels = c("Seawater at Day 0","Control", "WT", "dTDA"))



#Add max values to define where labels should be plottet according to y axis
emmeans_dat_all <- merge(emmeans_dat_all, Observed_Genera_tab_stat, by=c("Environment", "Treatment", "Day"))

####Figure S1A#####


Observed_Genera_rich <- p + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+40, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Environment)) 
  



############################################################
             #Barplots of microbial composition
##########################################################

#Class: Bacterial community composition ----- 
#Change Class level to P. inhibens for only P. inhibens Genus in order to plot P.inhibens individually
main_no_cont_NoOutliers_PhaeoInClass <- data.frame(tax_table(phylo_main_no_cont_NoOutliers))

main_no_cont_NoOutliers_PhaeoInClass$Class <- ifelse(main_no_cont_NoOutliers_PhaeoInClass$Genus=="Phaeobacter_inhibens", "P. inhibens OTU", main_no_cont_NoOutliers_PhaeoInClass$Class)
#Make new phyloseq table with new tax_table
phylo_main_no_cont_NoOutliers_PhaeoInClass <- phyloseq(otu_table(phylo_main_no_cont_NoOutliers), 
                                                       tax_table(as.matrix(main_no_cont_NoOutliers_PhaeoInClass)), sample_data(phylo_main_no_cont_NoOutliers))

# agglomerate at Class level
phylo_main_no_cont_NoOutliers_Class <- aggregate_taxa(phylo_main_no_cont_NoOutliers_PhaeoInClass, "Class")
# Transform to rel. abundance
phylo_main_no_cont_NoOutliers_Class_norm <- transform_sample_counts(phylo_main_no_cont_NoOutliers_Class, function(x) 100 * x/sum(x))
# Melt to long format
phylo_main_no_cont_NoOutliers_Class_norm_melt <- psmelt(phylo_main_no_cont_NoOutliers_Class_norm)
#Transform to percentage
phylo_main_no_cont_NoOutliers_Class_norm_melt <- aggregate(Abundance ~ OTU + Environment + Treatment + Day + Sample, phylo_main_no_cont_NoOutliers_Class_norm_melt, sum)
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum <- aggregate(Abundance ~ Environment + Treatment + Day + Sample, phylo_main_no_cont_NoOutliers_Class_norm_melt, sum)
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum <- left_join(phylo_main_no_cont_NoOutliers_Class_norm_melt, phylo_main_no_cont_NoOutliers_Class_norm_melt_sum, by=c("Environment", "Treatment", "Day", "Sample"))
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct <- phylo_main_no_cont_NoOutliers_Class_norm_melt_sum %>% 
  mutate(Abundance_percentage = Abundance.x / Abundance.y *100)
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct$Abundance_percentage <- round(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct$Abundance_percentage, 2)
#Find the 11 most abundant classes for the whole dataset including Phaeobacter (Only 12 colors for plotting)
Top12_Class <- phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct %>% group_by(OTU) %>%
  dplyr::summarise('Abundance_percentage' = mean(Abundance_percentage)) %>% arrange(desc(Abundance_percentage)) %>% head(11)

Top12_Class[[1]]

#Change Class names to "Others" if not among the 12 most abundant 
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 <- mutate(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct, Class = ifelse(OTU != "P. inhibens OTU" & OTU != "Gammaproteobacteria"&
                                                                                                                                       OTU != "Alphaproteobacteria" & OTU != "Flavobacteriia"&
                                                                                                                                       OTU != "Sphingobacteriia" & OTU != "Actinobacteria"&
                                                                                                                                       OTU != "Betaproteobacteria" & OTU != "Epsilonproteobacteria"&
                                                                                                                                       OTU != "dTDAproteobacteria" & OTU != "Planctomycetia"& 
                                                                                                                                       OTU != "Bacteria_Unclassified_Unclassified", "Others", OTU))

#Clean up names
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- ifelse(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class=="Bacteria_Unclassified_Unclassified", "Unclassified",  phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class)
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment <- ifelse(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment=="Seawater at Day 0", "Seawater",  phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment)
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Environment_1 <- ifelse(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Environment=="Seawater at Day 0", "Planktonic suspension", phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Environment)
#Reorder and clean up table
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Day <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Day, levels=c("0","1","4","10"))
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class , levels=c("Actinobacteria","Flavobacteriia","Planctomycetia","Sphingobacteriia","Alphaproteobacteria", "Betaproteobacteria", "dTDAproteobacteria","Epsilonproteobacteria","Gammaproteobacteria", "Others","Unclassified", "P. inhibens OTU"))
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment , levels=c("Seawater", "Control", "WT", "dTDA"))
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 <- arrange(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, Day, Treatment)

#Arrange Sample order according to Days and Treatments
Sample_ID_order <- distinct(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, Sample, .keep_all = TRUE)
Sample_ID_order <- unclass(arrange(Sample_ID_order, Day, Treatment)[5])$Sample

phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Sample <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Sample, levels=Sample_ID_order)




#Figure S2B -----
library(ggh4x)

# bar_plot_nolg <- ggplot(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, aes(x=Sample, y = Abundance_percentage, fill = Class)) +
#   geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
#   scale_fill_brewer(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
#   facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
#   theme_bw() +
#   theme(#axis.line = element_line(color='black'),
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position="none",
#     text = element_text(size=10))#,


#Order classes according to network figure legend
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class, 
                                        levels=c("Flavobacteriia","Gammaproteobacteria",
                                                "Alphaproteobacteria",
                                                "Actinobacteria", "Sphingobacteriia",
                                                "Planctomycetia",
                                                "Betaproteobacteria", "Epsilonproteobacteria", "Unclassified",
                                                "Others", "P. inhibens OTU"))

Paired_colors_expand_norm_shortcol <- c("#FFCE8F","#E6759E","#F89181","#343077","#f9c629","#f29d46","#ea3333","#f4722b","#fefd7d","#4e7677","#9C599E", "#f9c629","#2ace82")

bar_plot_nolg <- ggplot(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, aes(x=Sample, y = Abundance_percentage, fill = Class)) +
  geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
  scale_fill_manual(
    values = Paired_colors_expand_norm_shortcol,
    aesthetics = c("fill") 
  ) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="none",
    text = element_text(size=10))#,




legend_bar_plot <- cowplot::get_legend(bar_plot_nolg + theme(legend.position="bottom") + guides(fill = guide_legend(title.position = "top", ncol = 5)))

#Figure S1 Combined ----


# pdf(file = "Figures/Figure_S1_Alpha_diversity.pdf",   # The directory you want to save the file in
#     width = 11,
#     height = 7)
# 
# plot_grid(Observed_Genera_rich, bar_plot_nolg, legend_bar_plot, 
#           ncol=1,  labels = c('A', 'B', NA), rel_heights = c(3, 3, 1.5))
# 
# dev.off()

tiff("Figures/Figure_S1_Alpha_diversity.tiff", units="in", res=300,
     width = 11,
     height = 7)

plot_grid(Observed_Genera_rich, legend_p_plot, bar_plot_nolg, legend_bar_plot, 
          ncol=1,  labels = c('A', NA, 'B', NA), rel_heights = c(3, 0.5, 3, 1.5))


dev.off()

png("Figures/Figure_S1_Alpha_diversity.png", units="in", res=300,
     width = 11,
     height = 7)

plot_grid(Observed_Genera_rich, legend_p_plot, bar_plot_nolg, legend_bar_plot, 
          ncol=1,  labels = c('A', NA, 'B', NA), rel_heights = c(3, 0.5, 3, 1.5))


dev.off()


##Specific proportions------
phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 %>% filter(Treatment!="Control", Day=="4") %>% 
  group_by(Class) %>%  dplyr::summarise(max(Abundance_percentage), min(Abundance_percentage))

phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 %>% filter(Day=="1") %>% 
  group_by(Class) %>%  dplyr::summarise(max(Abundance_percentage), min(Abundance_percentage))

phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 %>% filter(Day!="1") %>% 
  group_by(Class) %>%  dplyr::summarise(max(Abundance_percentage), min(Abundance_percentage))

#Stats of richness and P. inhibens relative abuhndances ----
Observed_Genera_tab$Day <- factor(Observed_Genera_tab$Day, levels=c("0","1", "4", "10"))
Observed_Genera_tab_stat <- dplyr::group_by(Observed_Genera_tab, Environment, Treatment, Day) %>%
  dplyr::summarise(n=n(),'Observed Richness' = mean(Obs), sd = sd(Obs))
# Simpson_Genera_tab_stat <- dplyr::group_by(Observed_Genera_tab, Environment, Treatment, Day) %>% filter(Treatment!="Control") %>%
#   dplyr::summarise(n=n(),'Inverse Simpson' = median(Simpson_Genera), sd = sd(Simpson_Genera))

PhaeoRelAbundance <- t(as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_Class_norm))[rownames(as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_Class_norm)))=="P. inhibens", ])
PhaeoRelAbundance <- cbind(PhaeoRelAbundance, Observed_Genera_tab[,2:4])

PhaeoRelAbundance_stat <- PhaeoRelAbundance %>% group_by(Environment, Treatment, Day) %>%
  dplyr::summarise(n=n(),'P. inhibens Relative abundance' = median(`P. inhibens`), sd = sd(`P. inhibens`))

Stat_all <- cbind(Observed_Genera_tab_stat, PhaeoRelAbundance_stat[,5:6])

Stat_all[,-1] %>% filter(Treatment=="Control") %>% dplyr::summarise(mean(`P. inhibens Relative abundance`), sd(`P. inhibens Relative abundance`))

#qPCR absolute abudances for P. inhibens
qPCR_260221_stat <- qPCR_260221_stat[1:6]

colnames(qPCR_260221_stat) <- c("Environment", "Treatment","Day", "n","P. inhibens abundance","sd")
#change capital S in suspension to s.
#Do not plot "ns"; replace with NA 
qPCR_260221_stat$Environment <- ifelse(qPCR_260221_stat$Environment=="Planktonic Suspension", "Planktonic suspension",qPCR_260221_stat$Environment)

#Table 1 ----

Stat_all <- Stat_all %>%
  full_join(qPCR_260221_stat, by = c("Environment","Treatment","Day"), keep = FALSE)

write.table(Stat_all, file = "Stat_all.csv", dec=",")










#############################################
          #qPCR analysis with stats#

##################################################


setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/qPCR/")
library(emmeans)
library(lme4)
library(dplyr)
library(ggbeeswarm)
library(ggh4x)
library(dplyr)
library(tidyverse)

#load data
qPCR_260221 <- read.table("qPCR_040521_50ulelute.txt",  header = TRUE, sep = "\t")
qPCR_260221 <- as.data.frame(qPCR_260221)
#add day
qPCR_260221$Day <- factor(qPCR_260221$Day, levels = c("0","1", "4", "10"))
qPCR_260221 %>% filter(Day == "0", Env == "Biofilm", Treatment != "Control")

#Make means for qPCR technical reps
qPCR_260221 <- qPCR_260221 %>% group_by(Rep_tec, Treatment, Env, Day) %>% dplyr::summarise(CFU.cm2 = mean(CFU.cm2), Log = mean(Log))

qPCR_260221 %>% filter(Treatment!="Control", Day == "0", Env =="Biofilm")  %>% group_by(Env, Treatment, Day) %>%
  dplyr::summarise(n=n(),'P. inhibens abundance' = mean(Log), sd = sd(Log), max=max(Log))

colnames(qPCR_260221)[5] <- "Phaeobacter_CFU.cm2"
colnames(qPCR_260221)[6] <- "Phaeobacter_Log_CFU.cm2"

# qPCR_260221_T0 <- qPCR_260221 %>% filter(Day == "0")
# qPCR_260221_T0 <- qPCR_260221_T0[,c(1,6)]
# write.table(qPCR_260221_T0, "qPCR_260221_T0",append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

qPCR_absolutte_16S <- read.table("16S_absolutte_qPCR_2.txt", header = TRUE, sep = "\t")
qPCR_absolutte_16S <- as.data.frame(qPCR_absolutte_16S)

#setdiff(qPCR_260221[1],qPCR_absolutte_16S[1])

#qPCR_absolutte_16S_NAs <- cbind(qPCR_absolutte_16S_NAs, c(rep(NA, nrow(qPCR_absolutte_16S_NAs))))

#colnames(qPCR_absolutte_16S_NAs)[2] <- colnames(qPCR_absolutte_16S)[2]
#qPCR_absolutte_16S <- rbind(qPCR_absolutte_16S, qPCR_absolutte_16S_NAs)


qPCR_260221 <- left_join(qPCR_260221, qPCR_absolutte_16S, by = 'Rep_tec')


qPCR_260221 <- as.data.frame(qPCR_260221)
qPCR_260221$CFU_ml_total <- as.numeric(qPCR_260221$CFU_ml_total)

str(qPCR_260221)
#qPCR_260221 %>% filter(Day == "0")

colnames(qPCR_260221)
####Absolute abundance#####
qPCR_260221 %>% filter(Treatment != "Control")




#################Emmeans on linear models#####################
## Subset for Biofilm samples
qPCR_260221_Biofilm <- qPCR_260221 %>% filter(Env=="Biofilm") 

#Add replicate ID to dataframe
qPCR_260221_Biofilm$ID=gsub("(.*)(T[0-9])(.*)","\\1",qPCR_260221_Biofilm$Rep_tec)

#LMER Phaeobacter abundance
qPCR_260221_Biofilm_LMER_phaeo=lmer(Phaeobacter_Log_CFU.cm2~Treatment*Day + (1|ID), qPCR_260221_Biofilm)
anova(qPCR_260221_Biofilm_LMER_phaeo)

qPCR_260221_Biofilm_LMER_Phaeo_EM=emmeans(qPCR_260221_Biofilm_LMER_phaeo, ~Treatment | Day)
qPCR_260221_Biofilm_LMER_Phaeo_EM_stat <- contrast(qPCR_260221_Biofilm_LMER_Phaeo_EM, interaction = "pairwise", adjust = "Bonferroni")
qPCR_260221_Biofilm_LMER_Phaeo_EM_stat

#LMER 16S absolutte abundance or microbial load
qPCR_260221_Biofilm_LMER_16S=lmer(Log_CFU_ml_total~Treatment*Day + (1|ID), qPCR_260221_Biofilm)
anova(qPCR_260221_Biofilm_LMER_16S)

qPCR_260221_Biofilm_LMER_16S_EM=emmeans(qPCR_260221_Biofilm_LMER_16S, ~Treatment | Day)
qPCR_260221_Biofilm_LMER_16S_EM_stat <- contrast(qPCR_260221_Biofilm_LMER_16S_EM, interaction = "pairwise", adjust = "Bonferroni")
qPCR_260221_Biofilm_LMER_16S_EM_stat


## Subset for Planktonic suspension samples
qPCR_260221_SW <- qPCR_260221 %>% filter(Env=="Planktonic Suspension")
#Add replicate ID to dataframe
qPCR_260221_SW$ID=gsub("(.*)(T[0-9])(.*)","\\1",qPCR_260221_SW$Rep_tec)

#LMER
qPCR_260221_SW_LMER_phaeo=lmer(Phaeobacter_Log_CFU.cm2~Treatment*Day + (1|ID), qPCR_260221_SW)
anova(qPCR_260221_SW_LMER_phaeo)

qPCR_260221_SW_LMER_phaeo_EM=emmeans(qPCR_260221_SW_LMER_phaeo, ~Treatment | Day)
qPCR_260221_SW_LMER_phaeo_EM_stat <- contrast(qPCR_260221_SW_LMER_phaeo_EM, interaction = "pairwise", adjust = "Bonferroni")
qPCR_260221_SW_LMER_phaeo_EM_stat

#LMER 16S absolutte abundance or microbial load
qPCR_260221_SW_LMER_16S=lmer(Log_CFU_ml_total~Treatment*Day + (1|ID), qPCR_260221_SW)
anova(qPCR_260221_SW_LMER_16S)

qPCR_260221_SW_LMER_16S_EM=emmeans(qPCR_260221_SW_LMER_16S, ~Treatment | Day)
qPCR_260221_SW_LMER_16S_EM_stat <- contrast(qPCR_260221_SW_LMER_16S_EM, interaction = "pairwise", adjust = "Bonferroni")
qPCR_260221_SW_LMER_16S_EM_stat



#########Means and standard deviations for Phaeobacter abundances #########

#reorder treatments
qPCR_260221$Treatment <- factor(qPCR_260221$Treatment, levels = c("Control","WT", "dTDA"))

#qPCR absolute abundances for P. inhibens
qPCR_260221_stat_Phaeo <- qPCR_260221 %>% group_by(Env, Treatment, Day) %>%
  dplyr::summarise(n=n(),'P. inhibens abundance' = mean(Phaeobacter_Log_CFU.cm2), sd = sd(Phaeobacter_Log_CFU.cm2), max=max(Phaeobacter_Log_CFU.cm2))

head(qPCR_260221_stat_Phaeo, 20)

#Plots
#Define colors
colors <- c("Control" = "#80cdc1", "WT" = "#a6611a", "dTDA" = "#dfc27d")

qPCR_260221$Treatment <- factor(qPCR_260221$Treatment, levels=c("Control", "WT", "dTDA"))


#########Means and standard deviations for Abs abundances #########

#reorder treatments
qPCR_260221$Treatment <- factor(qPCR_260221$Treatment, levels = c("Control","WT", "dTDA"))

#qPCR absolute abundances for P. inhibens
qPCR_260221_stat_total_abs <- na.omit(qPCR_260221) %>%  group_by(Env, Treatment, Day) %>%
  dplyr::summarise(n=n(),'MeanLog_CFU_ml_total' = mean(Log_CFU_ml_total), sd = sd(Log_CFU_ml_total), max=max(Log_CFU_ml_total))




na.omit(qPCR_260221) %>%  group_by(Env, Treatment, Day) %>% mutate(sdev = sd(Log_CFU_ml_total) ) 

str(qPCR_260221)
head(qPCR_260221_stat_Phaeo, 20)

#Plots
#Define colors
colors <- c("Control" = "#80cdc1", "WT" = "#a6611a", "dTDA" = "#dfc27d")

qPCR_260221$Treatment <- factor(qPCR_260221$Treatment, levels=c("Control", "WT", "dTDA"))



## Plot Phaeo -------
p_phaeo <- ggplot((qPCR_260221), aes(x = Day, y = Phaeobacter_Log_CFU.cm2, fill = Treatment)) + geom_quasirandom(aes(col = Treatment), dodge.width=.8, cex=2)  +
stat_summary(fun.data=mean_se, 
geom="errorbar", color="black", width=0.2,  position=position_dodge(0.8)) + 
  stat_summary(fun=mean, geom="point", color="black", position=position_dodge(0.8)) +
  facet_grid(cols = vars(Env)) + 
  xlab("Day") + ylab("Log(CFUs/cm2 or ml)") + guides(fill = FALSE) +
  scale_colour_manual(
    values = colors,
    aesthetics = c("colour") 
  ) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "bottom")

#Add emmeans significance to bars
qPCR_260221_Biofilm_LMER_Phaeo_EM_stat_dat <-as.data.frame(summary(qPCR_260221_Biofilm_LMER_Phaeo_EM_stat)[c("Treatment_pairwise","p.value")])
qPCR_260221_SW_LMER_Phaeo_EM_stat_dat <-as.data.frame(summary(qPCR_260221_SW_LMER_phaeo_EM_stat)[c("Treatment_pairwise","p.value")])
qPCR_260221_Biofilm_LMER_Phaeo_EM_stat_dat$Env <- rep("Biofilm")
qPCR_260221_SW_LMER_Phaeo_EM_stat_dat$Env <- rep("Planktonic Suspension")
emmeans_dat_all <- rbind(qPCR_260221_Biofilm_LMER_Phaeo_EM_stat_dat, qPCR_260221_SW_LMER_Phaeo_EM_stat_dat)
emmeans_dat_all$signif <- c("a", "b", "b",
                            "a", "b", "b",
                            "a", "b", "b",
                            "a", "b","c",
                            "a", "b", "b",
                            "a", "b", "b",
                            "a", "b", "b",
                            "a", "b","c")
emmeans_dat_all$Treatment <- c("Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA")

emmeans_dat_all$Treatment <- factor(emmeans_dat_all$Treatment, levels = c("Control", "WT", "dTDA"))
emmeans_dat_all$Day <- c(rep("0", 3), rep("1", 3), rep("4", 3), rep("10", 3), rep("0", 3), rep("1", 3), rep("4", 3), rep("10", 3))
#Clean up table
emmeans_dat_all <- emmeans_dat_all[,3:6]
#Add max values to define where labels should be plottet according to y axis
emmeans_dat_all_Phaeo <- merge(emmeans_dat_all, qPCR_260221_stat_Phaeo, by=c("Env", "Treatment", "Day"))[,c(1:4, 8)]



# pdf(file = "../phyloseq/Figures/Figure_S1_qPCR_plot121121.pdf",   
#     width = 6,
#     height = 3) 

#p_phaeo + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Env))

# dev.off()

# tiff("../phyloseq/Figures/Figure_S1_qPCR_plot121121.tiff", units="in", width=6, height=3, res=300)

#p + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Env))

# dev.off()









## Plot Bacterial absolute abundances -------
p_16S <- ggplot((qPCR_260221), aes(x = Day, y = Log_CFU_ml_total, fill = Treatment)) + geom_quasirandom(aes(col = Treatment), dodge.width=.8, cex=2)  +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", color="black", width=0.2,  position=position_dodge(0.8)) + 
  stat_summary(fun=mean, geom="point", color="black", position=position_dodge(0.8)) +
  facet_grid(cols = vars(Env)) + 
  xlab("Day") + ylab("Log(CFUs/cm2 or ml)") + guides(fill = FALSE) +
  scale_colour_manual(
    values = colors,
    aesthetics = c("colour") 
  ) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "bottom")

#Add emmeans significance to bars
qPCR_260221_Biofilm_LMER_16S_EM_stat_dat <-as.data.frame(summary(qPCR_260221_Biofilm_LMER_16S_EM_stat)[c("Treatment_pairwise","p.value")])
qPCR_260221_SW_LMER_16S_EM_stat_dat <-as.data.frame(summary(qPCR_260221_SW_LMER_16S_EM_stat)[c("Treatment_pairwise","p.value")])
qPCR_260221_Biofilm_LMER_16S_EM_stat_dat$Env <- rep("Biofilm")
qPCR_260221_SW_LMER_16S_EM_stat_dat$Env <- rep("Planktonic Suspension")
emmeans_dat_all <- rbind(qPCR_260221_Biofilm_LMER_16S_EM_stat_dat, qPCR_260221_SW_LMER_16S_EM_stat_dat)
emmeans_dat_all$signif <- c(NA, NA, NA,
                            "a", "b", "b",
                            "a", "b", "b",
                            "a", "b","b",
                            NA, NA, NA,
                            "a", "a", "a",
                            "a", "a", "a",
                            "a", "a","a")
emmeans_dat_all$Treatment <- c("Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA",
                               "Control", "WT", "dTDA")

emmeans_dat_all$Treatment <- factor(emmeans_dat_all$Treatment, levels = c("Control", "WT", "dTDA"))
emmeans_dat_all$Day <- c(rep("0", 3), rep("1", 3), rep("4", 3), rep("10", 3), rep("0", 3), rep("1", 3), rep("4", 3), rep("10", 3))
#Clean up table
emmeans_dat_all <- emmeans_dat_all[,3:6]
#Add max values to define where labels should be plottet according to y axis
emmeans_dat_all_16S <- merge(emmeans_dat_all, qPCR_260221_stat_total_abs, by=c("Env", "Treatment", "Day"))[,c(1:4, 8)]

# pdf(file = "../phyloseq/Figures/Figure_S1_qPCR_plot121121.pdf",   
#     width = 6,
#     height = 3) 

p_16S  + geom_text(data=emmeans_dat_all_16S, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Env))

# dev.off()

# tiff("../phyloseq/Figures/Figure_S1_qPCR_plot121121.tiff", units="in", width=6, height=3, res=300)
# 
# p + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Env))
# 
# dev.off()


#### All abundances plot combined #####

#Combine dataframes
qPCR_260221_phaeo <- qPCR_260221[,c(1,2,3,4,6)]
emmeans_dat_all_Phaeo
qPCR_260221_phaeo$qPCR <- rep("Phaeobacter abundance")
colnames(qPCR_260221_phaeo)[5] <- "Log_CFU.cm2"

emmeans_dat_all_Phaeo$qPCR <- rep("Phaeobacter abundance")

qPCR_260221_16S <- qPCR_260221[,c(1,2,3,4,8)]
qPCR_260221_16S$qPCR <- rep("Total bacterial abundance")
colnames(qPCR_260221_16S)[5] <- "Log_CFU.cm2"

emmeans_dat_all_16S$qPCR <- rep("Total bacterial abundance")

qPCR_260221_long <- rbind(qPCR_260221_phaeo, qPCR_260221_16S)
qPCR_260221_long$qPCR <- factor(qPCR_260221_long$qPCR, levels = c("Total bacterial abundance", "Phaeobacter abundance"))

emmeans_dat_all <- rbind(emmeans_dat_all_16S, emmeans_dat_all_Phaeo)
emmeans_dat_all$qPCR <- factor(emmeans_dat_all$qPCR, levels = c("Total bacterial abundance", "Phaeobacter abundance"))



labels <- c("Total bacerial abundance", "italic('Phaeobacter') \n abundance") 



p_all <- ggplot((qPCR_260221_long), aes(x = Day, y = Log_CFU.cm2, fill = Treatment)) + geom_quasirandom(aes(col = Treatment), dodge.width=.8, cex=2)  +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", color="black", width=0.2,  position=position_dodge(0.8)) + 
  stat_summary(fun=mean, geom="point", color="black", position=position_dodge(0.8)) +
  facet_grid(cols = vars(qPCR), rows = vars(Env)) + 
  xlab("Day") + ylab("Log(CFUs/cm2 or ml)") + guides(fill = FALSE) +
  scale_colour_manual(
    values = colors,
    aesthetics = c("colour") 
  ) + theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "bottom",  text = element_text(size=10))



tiff("../phyloseq/Figures/Figure_1_qPCR_plot230222.tiff", units="in", width=4, height=6, res=300)

p_all + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(qPCR), rows = vars(Env)) 

dev.off()


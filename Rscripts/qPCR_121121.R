#############################################
          #qPCR analysis with stats#

##################################################


setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/qPCR/")
library(emmeans)
library(lme4)
library(dplyr)
library(ggbeeswarm)

#load data
qPCR_260221 <- read.table("qPCR_040521_50ulelute.txt",  header = TRUE, sep = "\t")
qPCR_260221 <- as.data.frame(qPCR_260221)
#add day
qPCR_260221$Day <- factor(qPCR_260221$Day, levels = c("0","1", "4", "10"))
#Make means for qPCR technical reps
qPCR_260221 <- qPCR_260221 %>% group_by(Rep_tec, Treatment, Env, Day) %>% dplyr::summarise(CFU.cm2 = mean(CFU.cm2), Log = mean(Log))

qPCR_260221 %>% filter(Treatment!="Control", Day == "0", Env =="Biofilm")  %>% group_by(Env, Treatment, Day) %>%
  dplyr::summarise(n=n(),'P. inhibens abundance' = mean(Log), sd = sd(Log), max=max(Log))




#################TWO-WAY ANOVA by emmeans#####################
## Subset for Biofilm samples
qPCR_260221_Biofilm <- qPCR_260221 %>% filter(Env=="Biofilm") 

#Add replicate ID to dataframe
qPCR_260221_Biofilm$ID=gsub("(.*)(T[0-9])(.*)","\\1",qPCR_260221_Biofilm$Rep_tec)

#LMER
qPCR_260221_Biofilm_LMER=lmer(Log~Treatment*Day + (1|ID), qPCR_260221_Biofilm)
anova(qPCR_260221_Biofilm_LMER)

qPCR_260221_Biofilm_EM_LMER=emmeans(qPCR_260221_Biofilm_LMER, ~Treatment | Day)
qPCR_260221_Biofilm_EM_LMER_stat <- contrast(qPCR_260221_Biofilm_EM_LMER, interaction = "pairwise", adjust = "Bonferroni")
class(summary(qPCR_260221_Biofilm_EM_LMER_stat))


## Subset for Planktonic suspension samples
qPCR_260221_SW <- qPCR_260221 %>% filter(Env=="Planktonic Suspension")
#Add replicate ID to dataframe
qPCR_260221_SW$ID=gsub("(.*)(T[0-9])(.*)","\\1",qPCR_260221_SW$Rep_tec)

#LMER
qPCR_260221_SW_LMER=lmer(Log~Treatment*Day + (1|ID), qPCR_260221_SW)
anova(qPCR_260221_SW_LMER)

qPCR_260221_SW_EM_LMER=emmeans(qPCR_260221_SW_LMER, ~Treatment | Day)
qPCR_260221_SW_EM_LMER_stat <- contrast(qPCR_260221_SW_EM_LMER, interaction = "pairwise", adjust = "Bonferroni")


#########Means and standard deviations#########

#reorder treatments
qPCR_260221$Treatment <- factor(qPCR_260221$Treatment, levels = c("Control","WT", "dTDA"))

#qPCR absolute abudances for P. inhibens
qPCR_260221_stat <- qPCR_260221 %>% group_by(Env, Treatment, Day) %>%
  dplyr::summarise(n=n(),'P. inhibens abundance' = mean(Log), sd = sd(Log), max=max(Log))

head(qPCR_260221_stat, 20)

#Plots
#Define colors
colors <- c("Control" = "#80cdc1", "WT" = "#a6611a", "dTDA" = "#dfc27d")



## Plot -------
p <- ggplot(qPCR_260221, aes(x = Day, y = Log, fill = Treatment)) + geom_quasirandom(aes(col = Treatment), dodge.width=.8, cex=2)  +
stat_summary(fun.data=mean_sd, 
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
qPCR_260221_Biofilm_EM_LMER_stat_dat <-as.data.frame(summary(qPCR_260221_Biofilm_EM_LMER_stat)[c("Treatment_pairwise","p.value")])
qPCR_260221_SW_EM_LMER_stat_dat <-as.data.frame(summary(qPCR_260221_SW_EM_LMER_stat)[c("Treatment_pairwise","p.value")])
qPCR_260221_Biofilm_EM_LMER_stat_dat$Env <- rep("Biofilm")
qPCR_260221_SW_EM_LMER_stat_dat$Env <- rep("Planktonic Suspension")
emmeans_dat_all <- rbind(qPCR_260221_Biofilm_EM_LMER_stat_dat, qPCR_260221_SW_EM_LMER_stat_dat)
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
emmeans_dat_all <- merge(emmeans_dat_all, qPCR_260221_stat, by=c("Env", "Treatment", "Day"))[,c(1:4, 8)]

pdf(file = "../phyloseq/Figures/Figure_S1_qPCR_plot121121.pdf",   
    width = 6,
    height = 3) 

p + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Env))

dev.off()

tiff("../phyloseq/Figures/Figure_S1_qPCR_plot121121.tiff", units="in", width=6, height=3, res=300)

p + geom_text(data=emmeans_dat_all, aes(x = Day, y = max+1, label = signif, group=Treatment), position=position_dodge(0.8)) + facet_grid(cols = vars(Env))

dev.off()




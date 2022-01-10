setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/phyloseq/")

library("SpiecEasi")
library("phyloseq")
library("NetCoMi")
library("metagMisc") #Split phyloseq object into different groups
library(microbiome)
library("dplyr")
library(phyloseq)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(magick)

#Load in the phyloseq object 
load('phylo_main_no_cont_NoOutliers.RData')
#For saving time load in the netComp comparisons with 1000 permutations
load('netcomp_WTDelta_SparCC_ALL_COMP_100_1000L_221021.RData')
#load('NetCoMi270921_permuted.RData')


#Aggregate to genus level
phylo_main_no_cont_NoOutliers_genus <- aggregate_taxa(phylo_main_no_cont_NoOutliers, "Genus")


#Biofilm Day 1: Comparing WT and dTDA network (th = 0.6) -----

phylo_main_no_cont_NoOutliers_genus_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers_genus, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_genus_Biofilm_day1 <- subset_samples(phylo_main_no_cont_NoOutliers_genus_Biofilm, Day=="1")

#Divide phyloseq object into each sample pool
phylo_main_no_cont_NoOutliers_genus_Biofilm_day1_split <- metagMisc::phyloseq_sep_variable(phylo_main_no_cont_NoOutliers_genus_Biofilm_day1, "Treatment") 


net_WTdTDA_SparCC_Biofilm_day1_100 <- netConstruct(data = phylo_main_no_cont_NoOutliers_genus_Biofilm_day1_split$WT, #Samples from the WT treatment at day 1 for the Biofilm
                                                    data2 = phylo_main_no_cont_NoOutliers_genus_Biofilm_day1_split$dTDA, #Samples from the dTDA treatment at day 1 for the Biofilm
                                                    measure = "sparcc",  #Correlation method
                                                    filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                                                    filtTaxPar = list(totalReads = 100, numbSamp = 5),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                                                    zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                                                    sparsMethod = "threshold", thresh = 0.6, #Filter all edges with strength below 0.6
                                                    dissFunc = "signed")


props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6 <- netAnalyze(net_WTdTDA_SparCC_Biofilm_day1_100, clustMethod = "cluster_fast_greedy",
                                                            hubQuant = 0.95, hubPar = c("degree", "between", "closeness"))




summary(props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6, showCentr = "all",  numbNodes = 5L, digits = 5L,
        groupNames = c("WT","dTDA")) 



# #Comparing WT and dTDA network plot
# plot_WTdTDA_SparCC_Biofilm_day1_100 <- plot(props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6,
#                                              # nodeFilter = "highestDegree",
#                                              # nodeFilterPar = 60,
#                                              nodeColor = "cluster", 
#                                              nodeSize = "eigenvector",
#                                              edgeWidth = 0.3,
#                                              # edgeFilter = "highestWeight",
#                                              # edgeFilterPar = 300,
#                                              borderCol = "gray40", 
#                                              title1 = "Network on Genus level with SparCC measures", 
#                                              showTitle = TRUE,
#                                              groupNames = c("WT","dTDA"),
#                                              cexTitle = 1, 
#                                              cexLabels = 0.8,
#                                              #labelLength = 15,
#                                              #labelPattern = c(120,"'",4),
#                                              labelScale = FALSE, repulsion = 0.9, 
#                                              shortenLabels = "none", 
#                                              #charToRm = "Bacteria_",
#                                              nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
#                                              hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
#                                              edgeTranspHigh = 40,
#                                              mar = c(2,2,3,1), set.seed(1000))
# 
# labels1 <- plot_WTdTDA_SparCC_Biofilm_day1_100$labels$labels1
# labels2 <- plot_WTdTDA_SparCC_Biofilm_day1_100$labels$labels2
# 
# nodeNames=list(labels1=labels1, labels2=labels2)
# 
# names(nodeNames[[1]])=labels1
# names(nodeNames[[2]])=labels2
# 
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# # pdf(file = "props_net_WTdTDA_SparCC_Biofilm_day1_100.pdf",   # The directory you want to save the file in
# #     width = 13, # The width of the plot in inches
# #     height = 5.5) # The height of the plot in inches
# 
# plot(props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6, 
#      labels = nodeNames,
#      nodeColor = "cluster",
#      nodeSize = "eigenvector",
#      edgeWidth = 0.3,
#      # edgeFilter = "highestWeight",
#      # edgeFilterPar = 300,
#      borderCol = "gray40", 
#      title1 = "Network on Genus level with spearman correlations", 
#      showTitle = TRUE,
#      groupNames = c("WT","dTDA"),
#      cexTitle = 1, 
#      #labelLength = 15,
#      #labelPattern = c(120,"'",4),
#      labelScale = FALSE, repulsion = 0.9, 
#      shortenLabels = "none", 
#      #charToRm = "Bacteria_",
#      nodeSizeSpread = 3, nodeTransp = 65, hubTransp = 50,
#      cexNodes = 1, edgeTranspLow = 0,
#      hubBorderWidth = 2, 
#      mar = c(2,2,3,1), set.seed(1000))

# dev.off()

### Comparing the sample similarity networks
#?netCompare


# netcomp_WTdTDA_SparCC_Biofilm_day1_100_1000L_0.6 <- netCompare(props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6, nPerm = 1000L,
#                                                                 permTest = TRUE, cores = 2L, 
#                                                                 storeAssoPerm = FALSE)


#load('NetCoMi_netcomp_WTdTDA_SparCC_Biofilm_day1_100_0.6_1000L_161021.RData')

#netcomp_WTdTDA_SparCC_Biofilm_day1_100_1000L_0.6$pvalDiffGlobalLCC$pvalavPath <- 1

summary(netcomp_WTdTDA_SparCC_Biofilm_day1_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 25L  ) #

#P. inhibens OTU sum of (+/-) edge connections
assomat_WT_SparCC_Biofilm_day1 <- net_WTdTDA_SparCC_Biofilm_day1_100$assoMat1 #WT
assomat_dTDA_SparCC_Biofilm_day1 <- net_WTdTDA_SparCC_Biofilm_day1_100$assoMat2 #dTDA

assomat_WT_SparCC_Biofilm_day1_Phaeobacter <- data.frame(assomat_WT_SparCC_Biofilm_day1)[which(dimnames(assomat_WT_SparCC_Biofilm_day1)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_Biofilm_day1_Phaeobacter_clean <- filter_all(assomat_WT_SparCC_Biofilm_day1_Phaeobacter, any_vars(. != 0& . != 1))

assomat_WT_SparCC_Biofilm_day1_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


assomat_dTDA_SparCC_Biofilm_day1_Phaeobacter <- as.data.frame(assomat_dTDA_SparCC_Biofilm_day1)[which(dimnames(assomat_dTDA_SparCC_Biofilm_day1)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_Biofilm_day1_Phaeobacter_clean <- filter_all(assomat_dTDA_SparCC_Biofilm_day1_Phaeobacter, any_vars(. != 0 & . != 1))

assomat_dTDA_SparCC_Biofilm_day1_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))






#Planktonic suspension Day 1: Comparing WT and dTDA network (th = 0.6) -----

phylo_main_no_cont_NoOutliers_genus_SW <- subset_samples(phylo_main_no_cont_NoOutliers_genus, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_genus_SW_day1 <- subset_samples(phylo_main_no_cont_NoOutliers_genus_SW, Day=="1")


phylo_main_no_cont_NoOutliers_genus_SW_day1_split <- metagMisc::phyloseq_sep_variable(phylo_main_no_cont_NoOutliers_genus_SW_day1, "Treatment") 


net_WTdTDA_SparCC_SW_day1_100 <- netConstruct(data = phylo_main_no_cont_NoOutliers_genus_SW_day1_split$WT,
                                                    data2 = phylo_main_no_cont_NoOutliers_genus_SW_day1_split$dTDA,
                                                    measure = "sparcc",
                                                    filtTax = c("totalReads", "numbSamp"),
                                                    filtTaxPar = list(totalReads = 100, numbSamp = 5),
                                                    zeroMethod = "none", normMethod = "none",
                                                    sparsMethod = "threshold", thresh = 0.6,
                                                    dissFunc = "signed")


props_net_WTdTDA_SparCC_SW_day1_100_0.6 <- netAnalyze(net_WTdTDA_SparCC_SW_day1_100, clustMethod = "cluster_fast_greedy",
                                                            hubQuant = 0.95, hubPar = c("degree", "between", "closeness"))



summary(props_net_WTdTDA_SparCC_SW_day1_100_0.6, showCentr = "all",  numbNodes = 5L, digits = 5L,
        groupNames = c("WT","dTDA")) 

#Comparing WT and dTDA network plot
# plot_WTdTDA_SparCC_SW_day1_100 <- plot(props_net_WTdTDA_SparCC_SW_day1_100_0.6,
#                                              # nodeFilter = "highestDegree",
#                                              # nodeFilterPar = 60,
#                                              nodeColor = "cluster", 
#                                              nodeSize = "eigenvector",
#                                              edgeWidth = 0.3,
#                                              # edgeFilter = "highestWeight",
#                                              # edgeFilterPar = 300,
#                                              borderCol = "gray40", 
#                                              title1 = "Network on Genus level with SparCC measures", 
#                                              showTitle = TRUE,
#                                              groupNames = c("WT","dTDA"),
#                                              cexTitle = 1, 
#                                              cexLabels = 0.8,
#                                              #labelLength = 15,
#                                              #labelPattern = c(120,"'",4),
#                                              labelScale = FALSE, repulsion = 0.9, 
#                                              shortenLabels = "none", 
#                                              #charToRm = "Bacteria_",
#                                              nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
#                                              hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
#                                              edgeTranspHigh = 40,
#                                              mar = c(2,2,3,1), set.seed(1000))
# 
# labels1 <- plot_WTdTDA_SparCC_SW_day1_100$labels$labels1
# labels2 <- plot_WTdTDA_SparCC_SW_day1_100$labels$labels2
# 
# nodeNames=list(labels1=labels1, labels2=labels2)
# 
# names(nodeNames[[1]])=labels1
# names(nodeNames[[2]])=labels2
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# # pdf(file = "props_net_WTdTDA_SparCC_SW_day1_100.pdf",   # The directory you want to save the file in
# #     width = 13, # The width of the plot in inches
# #     height = 5.5) # The height of the plot in inches
# 
# plot(props_net_WTdTDA_SparCC_SW_day1_100_0.6, 
#      labels = nodeNames,
#      nodeColor = "cluster",
#      nodeSize = "eigenvector",
#      edgeWidth = 0.3,
#      # edgeFilter = "highestWeight",
#      # edgeFilterPar = 300,
#      borderCol = "gray40", 
#      title1 = "Network on Genus level with spearman correlations", 
#      showTitle = TRUE,
#      groupNames = c("WT","dTDA"),
#      cexTitle = 1, 
#      #labelLength = 15,
#      #labelPattern = c(120,"'",4),
#      labelScale = FALSE, repulsion = 0.9, 
#      shortenLabels = "none", 
#      #charToRm = "Bacteria_",
#      nodeSizeSpread = 3, nodeTransp = 65, hubTransp = 50,
#      cexNodes = 1, edgeTranspLow = 0,
#      hubBorderWidth = 2, 
#      mar = c(2,2,3,1), set.seed(1000))

# dev.off()

### Comparing the sample similarity networks
#?netCompare

# netcomp_WTdTDA_spear_SW_day1_1000_1000L <- netCompare(props_net_WTdTDA_SparCC_SW_day1_1000, nPerm = 1000L,
#                                                        permTest = TRUE, cores = 2L, 
#                                                        storeAssoPerm = FALSE)

#netcomp_WTdTDA_SparCC_SW_day1_100_1000L_0.6 <- netCompare(props_net_WTdTDA_SparCC_SW_day1_100_0.6, nPerm = 1000L,
#                                                                 permTest = TRUE, cores = 2L, 
#                                                                 storeAssoPerm = FALSE)


#load('NetCoMi_netcomp_WTdTDA_SparCC_SW_day1_100_0.6_1000L_161021.RData')
netcomp_WTdTDA_SparCC_SW_day1_100_1000L_0.6$pvalDiffGlobalLCC$pvalavPath <- 1

summary(netcomp_WTdTDA_SparCC_SW_day1_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 25L  ) #


#P. inhibens OTU sum of (+/-) edge connections
assomat_WT_SparCC_SW_day1 <- net_WTdTDA_SparCC_SW_day1_100$assoMat1 #WT
assomat_dTDA_SparCC_SW_day1 <- net_WTdTDA_SparCC_SW_day1_100$assoMat2 #dTDA

assomat_WT_SparCC_SW_day1_Phaeobacter <- data.frame(assomat_WT_SparCC_SW_day1)[which(dimnames(assomat_WT_SparCC_SW_day1)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_SW_day1_Phaeobacter_clean <- filter_all(assomat_WT_SparCC_SW_day1_Phaeobacter, any_vars(. != 0& . != 1))

assomat_WT_SparCC_SW_day1_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


assomat_dTDA_SparCC_SW_day1_Phaeobacter <- as.data.frame(assomat_dTDA_SparCC_SW_day1)[which(dimnames(assomat_dTDA_SparCC_SW_day1)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_SW_day1_Phaeobacter_clean <- filter_all(assomat_dTDA_SparCC_SW_day1_Phaeobacter, any_vars(. != 0 & . != 1))

assomat_dTDA_SparCC_SW_day1_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))








#Biofilm Day 4: Comparing WT and dTDA network (th = 0.6) -----

phylo_main_no_cont_NoOutliers_genus_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers_genus, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_genus_Biofilm_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_genus_Biofilm, Day=="4")


phylo_main_no_cont_NoOutliers_genus_Biofilm_day4_split <- metagMisc::phyloseq_sep_variable(phylo_main_no_cont_NoOutliers_genus_Biofilm_day4, "Treatment") 


net_WTdTDA_SparCC_Biofilm_day4_100 <- netConstruct(data = phylo_main_no_cont_NoOutliers_genus_Biofilm_day4_split$WT,
                                                    data2 = phylo_main_no_cont_NoOutliers_genus_Biofilm_day4_split$dTDA,
                                                    measure = "sparcc",
                                                    filtTax = c("totalReads", "numbSamp"),
                                                    filtTaxPar = list(totalReads = 100, numbSamp = 5),
                                                    zeroMethod = "none", normMethod = "none",
                                                    sparsMethod = "threshold", thresh = 0.6,
                                                    dissFunc = "signed")


props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6 <- netAnalyze(net_WTdTDA_SparCC_Biofilm_day4_100, clustMethod = "cluster_fast_greedy",
                                                        hubQuant = 0.95, hubPar = c("eigenvector"))

summary(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6, showCentr = "all",  numbNodes = 5L, digits = 5L,
        groupNames = c("WT","dTDA")) 



#Comparing WT and dTDA network plot
plot_WTdTDA_SparCC_Biofilm_day4_100 <- plot(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6,
                                             sameLayout = TRUE, 
                                             # nodeFilter = "highestDegree",
                                             # nodeFilterPar = 60,
                                             nodeColor = "cluster", 
                                             nodeSize = "eigenvector",
                                             edgeWidth = 0.3,
                                             # edgeFilter = "highestWeight",
                                             # edgeFilterPar = 300,
                                             borderCol = "gray40", 
                                             title1 = "Network on Genus level with SparCC measures", 
                                             showTitle = TRUE,
                                             groupNames = c("WT","dTDA"),
                                             cexTitle = 1, 
                                             cexLabels = 0.8,
                                             #labelLength = 15,
                                             #labelPattern = c(120,"'",4),
                                             labelScale = FALSE, repulsion = 0.9, 
                                             shortenLabels = "none", 
                                             #charToRm = "Bacteria_",
                                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                                             edgeTranspHigh = 40,
                                             mar = c(2,2,3,1), set.seed(1000))



labels1 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels1
labels2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2

nodeNames=list(labels1=labels1, labels2=labels2)

names(nodeNames[[1]])=labels1
names(nodeNames[[2]])=labels2

nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])

nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "Genus"

dim(nodecol1_dat)
length(nodecol1)


tax_order <- as.data.frame(tax_table(phylo_main_no_cont_NoOutliers_genus_Biofilm_day4))

nodecol1_dat_all <- left_join(nodecol1_dat, tax_order, by ='Genus')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$Class

names(nodecol_order) <- nodecol1

nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 

Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 

# show_col(Paired_colors_expand_norm)
# 
# Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))

Paired_colors_expand_norm <- (c(Sea_Sunset_colors,Tropical_colors,Mountain_colors, Peep_beach_colors))

show_col(Paired_colors_expand_norm)
#show_col(brewer.pal(n = 12, name = "Paired"))


# nodecol_order <- factor(nodecol_order, levels=c("Actinobacteria","Flavobacteriia",
#                                                 "Planctomycetia","Sphingobacteriia",
#                                                 "Alphaproteobacteria", "Caldilineae", 
#                                                 "Opitutae","Deinococci",
#                                                 "Gammaproteobacteria", "Cytophagia", 
#                                                 "Unclassified","Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria", 
#                                                 "Verrucomicrobiae"))
nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
                                                "Cytophagia","Alphaproteobacteria",
                                                "Actinobacteria", "Sphingobacteriia",
                                                "Caldilineae","Deltaproteobacteria",
                                                "Opitutae", "Planctomycetia",
                                                "Deinococci",
                                              "Holophagae", "Chlamydiia", "Phycisphaerae",
                                                "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
                                                "Verrucomicrobiae"))
#Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))

Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))


show_col(Paired_colors_expand_norm)

######## Figure 3A ----

tiff("Figures/Figure_3A_Network.tiff", units="in", width=8.5, height=5, res=300)

plot(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "feature", featVecCol = nodecol_order, colorVec = Paired_colors_expand_norm,
     nodeSize = "eigenvector",
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40", 
     showTitle = TRUE,
     groupNames = c("WT","dTDA"),
     cexTitle = 1, 
     # posCol = "#70e599", 
     # negCol = "red",
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 1.5, 
     shortenLabels = "none", 
     #charToRm = "Bacteria_",
     nodeSizeSpread = 3, nodeTransp = 10, edgeTransp = 50,
     cexNodes = 1, edgeTranspLow = 0,
     hubBorderWidth = 2, 
     mar = c(2,2,3,1), set.seed(1000))

dev.off()

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- NetCoMi:::colToTransp(Paired_colors_expand_norm, 10)
#Save legends for final plot

tiff("Figures/Figure_3A_Network_legends.tiff", units="in", width=1.7 ,height=2, res=300)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0,2), ylim=c(0,2))
par(mar = c(0, 0, 0, 0), ps=6)


legend(-3.2, 12, pt.cex = 1, title = "Class",title.adj=0,
       legend=levels(nodecol_order), col = phylcol_transp, bty = "n", pch=16, y.intersp = 3, x.intersp = 1.7)


legend(0, -3,  title = "Estimated correlation", title.adj=0,
       legend = c("+","-"), lty = 1, lwd = 1.5, col = c("darkgreen","red"), 
       bty = "n", horiz = TRUE)

dev.off()


  #P. inhibens OTU sum of (+/-) edge connections
assomat_WT_SparCC_Biofilm_day4 <- net_WTdTDA_SparCC_Biofilm_day4_100$assoMat1 #WT
assomat_dTDA_SparCC_Biofilm_day4 <- net_WTdTDA_SparCC_Biofilm_day4_100$assoMat2 #dTDA

assomat_WT_SparCC_Biofilm_day4_Phaeobacter <- data.frame(assomat_WT_SparCC_Biofilm_day4)[which(dimnames(assomat_WT_SparCC_Biofilm_day4)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_Biofilm_day4_Phaeobacter_clean <- filter_all(assomat_WT_SparCC_Biofilm_day4_Phaeobacter, any_vars(. != 0& . != 1))

assomat_WT_SparCC_Biofilm_day4_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


assomat_dTDA_SparCC_Biofilm_day4_Phaeobacter <- as.data.frame(assomat_dTDA_SparCC_Biofilm_day4)[which(dimnames(assomat_dTDA_SparCC_Biofilm_day4)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_Biofilm_day4_Phaeobacter_clean <- filter_all(assomat_dTDA_SparCC_Biofilm_day4_Phaeobacter, any_vars(. != 0 & . != 1))

assomat_dTDA_SparCC_Biofilm_day4_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))





# pdf(file = "props_net_WTdTDA_SparCC_Biofilm_day4_100.pdf",   # The directory you want to save the file in
#     width = 8, # The width of the plot in inches
#     height = 5.5) # The height of the plot in inches
# 
# 
# 
# plot(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6, 
#      sameLayout = FALSE, 
#      labels = nodeNames,
#      nodeColor = "feature", featVecCol = nodecol_order, colorVec = Paired_colors_expand,
#      nodeSize = "eigenvector",
#      edgeWidth = 0.3,
#      # edgeFilter = "highestWeight",
#      # edgeFilterPar = 300,
#      borderCol = "gray40", 
#      title1 = "Network on Genus level with spearman correlations", 
#      showTitle = TRUE,
#      groupNames = c("WT","dTDA"),
#      cexTitle = 1, 
#      #labelLength = 15,
#      #labelPattern = c(120,"'",4),
#      labelScale = FALSE, repulsion = 0.9, 
#      shortenLabels = "none", 
#      #charToRm = "Bacteria_",
#      nodeSizeSpread = 3, nodeTransp = 10, edgeTransp = 50,
#      cexNodes = 1, edgeTranspLow = 0,
#      hubBorderWidth = 2, 
#      mar = c(2,2,3,1), set.seed(1000))
# 
# 
# dev.off()

### Comparing the sample similarity networks

# #netcomp_WTdTDA_SparCC_Biofilm_day4_100_1000L_0.6 <- netCompare(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6, nPerm = 1000L,
#                                                             permTest = TRUE, cores = 2L, 
#                                                             storeAssoPerm = FALSE)

summary(netcomp_WTdTDA_SparCC_Biofilm_day4_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 25L  ) #

#Planktonic suspension Day 4: Comparing WT and dTDA network (th = 0.6) -----
phylo_main_no_cont_NoOutliers_genus_SW <- subset_samples(phylo_main_no_cont_NoOutliers_genus, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_genus_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_genus_SW, Day=="4")

phylo_main_no_cont_NoOutliers_genus_SW_day4_split <- metagMisc::phyloseq_sep_variable(phylo_main_no_cont_NoOutliers_genus_SW_day4, "Treatment") 

net_WTdTDA_SparCC_SW_day4_100 <- netConstruct(data = phylo_main_no_cont_NoOutliers_genus_SW_day4_split$WT,
                                               data2 = phylo_main_no_cont_NoOutliers_genus_SW_day4_split$dTDA,
                                               measure = "sparcc",
                                               filtTax = c("totalReads", "numbSamp"),
                                               filtTaxPar = list(totalReads = 100, numbSamp = 5),
                                               zeroMethod = "none", normMethod = "none",
                                               sparsMethod = "threshold", thresh = 0.6,
                                               dissFunc = "signed")


props_net_WTdTDA_SparCC_SW_day4_100 <- netAnalyze(net_WTdTDA_SparCC_SW_day4_100, clustMethod = "cluster_fast_greedy",
                                                   hubQuant = 0.95, hubPar = c("degree", "between", "closeness"))

summary(props_net_WTdTDA_SparCC_SW_day4_100, showCentr = "all",  numbNodes = 5L, digits = 5L,
        groupNames = c("WT","dTDA")) 

#Comparing WT and dTDA network plot
# plot_WTdTDA_SparCC_SW_day4_100 <- plot(props_net_WTdTDA_SparCC_SW_day4_100,
#                                         # nodeFilter = "highestDegree",
#                                         # nodeFilterPar = 60,
#                                         nodeColor = "cluster", 
#                                         nodeSize = "eigenvector",
#                                         edgeWidth = 0.3,
#                                         # edgeFilter = "highestWeight",
#                                         # edgeFilterPar = 300,
#                                         borderCol = "gray40", 
#                                         title1 = "Network on Genus level with SparCC measures", 
#                                         showTitle = TRUE,
#                                         groupNames = c("WT","dTDA"),
#                                         cexTitle = 1, 
#                                         cexLabels = 0.8,
#                                         #labelLength = 15,
#                                         #labelPattern = c(120,"'",4),
#                                         labelScale = FALSE, repulsion = 0.9, 
#                                         shortenLabels = "none", 
#                                         #charToRm = "Bacteria_",
#                                         nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
#                                         hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
#                                         edgeTranspHigh = 40,
#                                         mar = c(2,2,3,1), set.seed(1000))
# 
# labels1 <- plot_WTdTDA_SparCC_SW_day4_100$labels$labels1
# labels2 <- plot_WTdTDA_SparCC_SW_day4_100$labels$labels2
# 
# nodeNames=list(labels1=labels1, labels2=labels2)
# 
# names(nodeNames[[1]])=labels1
# names(nodeNames[[2]])=labels2
# 
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# # pdf(file = "props_net_WTdTDA_SparCC_SW_day4_100.pdf",   # The directory you want to save the file in
# #     width = 13, # The width of the plot in inches
# #     height = 5.5) # The height of the plot in inches
# 
# plot(props_net_WTdTDA_SparCC_SW_day4_100, 
#      labels = nodeNames,
#      nodeColor = "cluster",
#      nodeSize = "eigenvector",
#      edgeWidth = 0.3,
#      # edgeFilter = "highestWeight",
#      # edgeFilterPar = 300,
#      borderCol = "gray40", 
#      title1 = "Network on Genus level with spearman correlations", 
#      showTitle = TRUE,
#      groupNames = c("WT","dTDA"),
#      cexTitle = 1, 
#      #labelLength = 15,
#      #labelPattern = c(120,"'",4),
#      labelScale = FALSE, repulsion = 0.9, 
#      shortenLabels = "none", 
#      #charToRm = "Bacteria_",
#      nodeSizeSpread = 3, nodeTransp = 65, hubTransp = 50,
#      cexNodes = 1, edgeTranspLow = 0,
#      hubBorderWidth = 2, 
#      mar = c(2,2,3,1), set.seed(1000))

# dev.off()

### Comparing the sample similarity networks

# #netcomp_WTdTDA_SparCC_SW_day4_100_1000L_0.6 <- netCompare(props_net_WTdTDA_SparCC_SW_day4_100, nPerm = 1000L,
#                                                        permTest = TRUE, cores = 2L, 
#                                                        storeAssoPerm = FALSE)
# 


summary(netcomp_WTdTDA_SparCC_SW_day4_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 5L, digits = 3L,) #

#P. inhibens OTU sum of (+/-) edge connections
assomat_WT_SparCC_SW_day4 <- net_WTdTDA_SparCC_SW_day4_100$assoMat1 #WT
assomat_dTDA_SparCC_SW_day4 <- net_WTdTDA_SparCC_SW_day4_100$assoMat2 #dTDA

assomat_WT_SparCC_SW_day4_Phaeobacter <- data.frame(assomat_WT_SparCC_SW_day4)[which(dimnames(assomat_WT_SparCC_SW_day4)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_SW_day4_Phaeobacter_clean <- filter_all(assomat_WT_SparCC_SW_day4_Phaeobacter, any_vars(. != 0& . != 1))

assomat_WT_SparCC_SW_day4_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


assomat_dTDA_SparCC_SW_day4_Phaeobacter <- as.data.frame(assomat_dTDA_SparCC_SW_day4)[which(dimnames(assomat_dTDA_SparCC_SW_day4)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_SW_day4_Phaeobacter_clean <- filter_all(assomat_dTDA_SparCC_SW_day4_Phaeobacter, any_vars(. != 0 & . != 1))

assomat_dTDA_SparCC_SW_day4_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))







#Biofilm Day 10: Comparing WT and dTDA network (th = 0.6) -----

phylo_main_no_cont_NoOutliers_genus_Biofilm <- subset_samples(phylo_main_no_cont_NoOutliers_genus, Environment=="Biofilm")
phylo_main_no_cont_NoOutliers_genus_Biofilm_day10 <- subset_samples(phylo_main_no_cont_NoOutliers_genus_Biofilm, Day=="10")


phylo_main_no_cont_NoOutliers_genus_Biofilm_day10_split <- metagMisc::phyloseq_sep_variable(phylo_main_no_cont_NoOutliers_genus_Biofilm_day10, "Treatment") 


net_WTdTDA_SparCC_Biofilm_day10_100 <- netConstruct(data = phylo_main_no_cont_NoOutliers_genus_Biofilm_day10_split$WT,
                                                    data2 = phylo_main_no_cont_NoOutliers_genus_Biofilm_day10_split$dTDA,
                                                    measure = "sparcc",
                                                    filtTax = c("totalReads", "numbSamp"),
                                                    filtTaxPar = list(totalReads = 100, numbSamp = 5),
                                                    zeroMethod = "none", normMethod = "none",
                                                    sparsMethod = "threshold", thresh = 0.6,
                                                    dissFunc = "signed")


props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6 <- netAnalyze(net_WTdTDA_SparCC_Biofilm_day10_100, clustMethod = "cluster_fast_greedy",
                                                            hubQuant = 0.95, hubPar = c("degree", "between", "closeness"))



summary(props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6, showCentr = "all",  numbNodes = 5L, digits = 5L,
        groupNames = c("WT","dTDA")) 

#Comparing WT and dTDA network plot
# plot_WTdTDA_SparCC_Biofilm_day10_100 <- plot(props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6,
#                                              # nodeFilter = "highestDegree",
#                                              # nodeFilterPar = 60,
#                                              nodeColor = "cluster", 
#                                              nodeSize = "eigenvector",
#                                              edgeWidth = 0.3,
#                                              # edgeFilter = "highestWeight",
#                                              # edgeFilterPar = 300,
#                                              borderCol = "gray40", 
#                                              title1 = "Network on Genus level with SparCC measures", 
#                                              showTitle = TRUE,
#                                              groupNames = c("WT","dTDA"),
#                                              cexTitle = 1, 
#                                              cexLabels = 0.8,
#                                              #labelLength = 15,
#                                              #labelPattern = c(120,"'",4),
#                                              labelScale = FALSE, repulsion = 0.9, 
#                                              shortenLabels = "none", 
#                                              #charToRm = "Bacteria_",
#                                              nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
#                                              hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
#                                              edgeTranspHigh = 40,
#                                              mar = c(2,2,3,1), set.seed(1000))
# labels1 <- plot_WTdTDA_SparCC_Biofilm_day10_100$labels$labels1
# labels2 <- plot_WTdTDA_SparCC_Biofilm_day10_100$labels$labels2
# 
# nodeNames=list(labels1=labels1, labels2=labels2)
# 
# names(nodeNames[[1]])=labels1
# names(nodeNames[[2]])=labels2
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# pdf(file = "props_net_WTdTDA_SparCC_Biofilm_day10_100.pdf",   # The directory you want to save the file in
#     width = 13, # The width of the plot in inches
#     height = 5.5) # The height of the plot in inches
# 
# plot(props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6, 
#      labels = nodeNames,
#      nodeColor = "cluster",
#      nodeSize = "eigenvector",
#      edgeWidth = 0.3,
#      # edgeFilter = "highestWeight",
#      # edgeFilterPar = 300,
#      borderCol = "gray40", 
#      title1 = "Network on Genus level with spearman correlations", 
#      showTitle = TRUE,
#      groupNames = c("WT","dTDA"),
#      cexTitle = 1, 
#      #labelLength = 15,
#      #labelPattern = c(120,"'",4),
#      labelScale = FALSE, repulsion = 0.9, 
#      shortenLabels = "none", 
#      #charToRm = "Bacteria_",
#      nodeSizeSpread = 3, nodeTransp = 65, hubTransp = 50,
#      cexNodes = 1, edgeTranspLow = 0,
#      hubBorderWidth = 2, 
#      mar = c(2,2,3,1), set.seed(1000))
# 
# dev.off()


### Comparing the sample similarity networks

# netcomp_WTdTDA_SparCC_Biofilm_day10_100_1000L_0.6 <- netCompare(props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6, nPerm = 1000L,
#                                                                 permTest = TRUE, cores = 2L, 
#                                                                 storeAssoPerm = FALSE)
# 


summary(netcomp_WTdTDA_SparCC_Biofilm_day10_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 25L  ) #



#P. inhibens OTU sum of (+/-) edge connections
assomat_WT_SparCC_Biofilm_day10 <- net_WTdTDA_SparCC_Biofilm_day10_100$assoMat1 #WT
assomat_dTDA_SparCC_Biofilm_day10 <- net_WTdTDA_SparCC_Biofilm_day10_100$assoMat2 #dTDA

assomat_WT_SparCC_Biofilm_day10_Phaeobacter <- data.frame(assomat_WT_SparCC_Biofilm_day10)[which(dimnames(assomat_WT_SparCC_Biofilm_day10)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_Biofilm_day10_Phaeobacter_clean <- filter_all(assomat_WT_SparCC_Biofilm_day10_Phaeobacter, any_vars(. != 0& . != 1))

assomat_WT_SparCC_Biofilm_day10_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


assomat_dTDA_SparCC_Biofilm_day10_Phaeobacter <- as.data.frame(assomat_dTDA_SparCC_Biofilm_day10)[which(dimnames(assomat_dTDA_SparCC_Biofilm_day10)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_Biofilm_day10_Phaeobacter_clean <- filter_all(assomat_dTDA_SparCC_Biofilm_day10_Phaeobacter, any_vars(. != 0 & . != 1))

assomat_dTDA_SparCC_Biofilm_day10_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))



#SW Day 10: Comparing WT and dTDA network (th = 0.6) -----

phylo_main_no_cont_NoOutliers_genus_SW <- subset_samples(phylo_main_no_cont_NoOutliers_genus, Environment=="Planktonic suspension")
phylo_main_no_cont_NoOutliers_genus_SW_day10 <- subset_samples(phylo_main_no_cont_NoOutliers_genus_SW, Day=="10")


phylo_main_no_cont_NoOutliers_genus_SW_day10_split <- metagMisc::phyloseq_sep_variable(phylo_main_no_cont_NoOutliers_genus_SW_day10, "Treatment") 


net_WTdTDA_SparCC_SW_day10_100 <- netConstruct(data = phylo_main_no_cont_NoOutliers_genus_SW_day10_split$WT,
                                                     data2 = phylo_main_no_cont_NoOutliers_genus_SW_day10_split$dTDA,
                                                     measure = "sparcc",
                                                     filtTax = c("totalReads", "numbSamp"),
                                                     filtTaxPar = list(totalReads = 100, numbSamp = 5),
                                                     zeroMethod = "none", normMethod = "none",
                                                     sparsMethod = "threshold", thresh = 0.6,
                                                     dissFunc = "signed")


props_net_WTdTDA_SparCC_SW_day10_100_0.6 <- netAnalyze(net_WTdTDA_SparCC_SW_day10_100, clustMethod = "cluster_fast_greedy",
                                                             hubQuant = 0.95, hubPar = c("degree", "between", "closeness"))


summary(props_net_WTdTDA_SparCC_SW_day10_100_0.6, showCentr = "all",  numbNodes = 5L, digits = 5L,
        groupNames = c("WT","dTDA")) 

# plot_WTdTDA_SparCC_SW_day10_100 <- plot(props_net_WTdTDA_SparCC_SW_day10_100_0.6,
#                                               # nodeFilter = "highestDegree",
#                                               # nodeFilterPar = 60,
#                                               nodeColor = "cluster", 
#                                               nodeSize = "eigenvector",
#                                               edgeWidth = 0.3,
#                                               # edgeFilter = "highestWeight",
#                                               # edgeFilterPar = 300,
#                                               borderCol = "gray40", 
#                                               title1 = "Network on Genus level with SparCC measures", 
#                                               showTitle = TRUE,
#                                               groupNames = c("WT","dTDA"),
#                                               cexTitle = 1, 
#                                               cexLabels = 0.8,
#                                               #labelLength = 15,
#                                               #labelPattern = c(120,"'",4),
#                                               labelScale = FALSE, repulsion = 0.9, 
#                                               shortenLabels = "none", 
#                                               #charToRm = "Bacteria_",
#                                               nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
#                                               hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
#                                               edgeTranspHigh = 40,
#                                               mar = c(2,2,3,1), set.seed(1000))
# 
# labels1 <- plot_WTdTDA_SparCC_SW_day10_100$labels$labels1
# labels2 <- plot_WTdTDA_SparCC_SW_day10_100$labels$labels2
# 
# nodeNames=list(labels1=labels1, labels2=labels2)
# 
# names(nodeNames[[1]])=labels1
# names(nodeNames[[2]])=labels2
# 
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# # pdf(file = "props_net_WTdTDA_SparCC_SW_day10_100.pdf",   # The directory you want to save the file in
# #     width = 13, # The width of the plot in inches
# #     height = 5.5) # The height of the plot in inches
# 
# plot(props_net_WTdTDA_SparCC_SW_day10_100_0.6, 
#      labels = nodeNames,
#      nodeColor = "cluster",
#      nodeSize = "eigenvector",
#      edgeWidth = 0.3,
#      # edgeFilter = "highestWeight",
#      # edgeFilterPar = 300,
#      borderCol = "gray40", 
#      title1 = "Network on Genus level with spearman correlations", 
#      showTitle = TRUE,
#      groupNames = c("WT","dTDA"),
#      cexTitle = 1, 
#      #labelLength = 15,
#      #labelPattern = c(120,"'",4),
#      labelScale = FALSE, repulsion = 0.9, 
#      shortenLabels = "none", 
#      #charToRm = "Bacteria_",
#      nodeSizeSpread = 3, nodeTransp = 65, hubTransp = 50,
#      cexNodes = 1, edgeTranspLow = 0,
#      hubBorderWidth = 2, 
#      mar = c(2,2,3,1), set.seed(1000))

# dev.off()

# netcomp_WTdTDA_SparCC_SW_day10_100_1000L_0.6 <- netCompare(props_net_WTdTDA_SparCC_SW_day10_100_0.6, nPerm = 1000L,
#                                                                  permTest = TRUE, cores = 2L, 
#                                                                  storeAssoPerm = FALSE)
# 
summary(netcomp_WTdTDA_SparCC_SW_day10_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 25L  ) #


#P. inhibens OTU sum of (+/-) edge connections
assomat_WT_SparCC_SW_day10 <- net_WTdTDA_SparCC_SW_day10_100$assoMat1 #WT
assomat_dTDA_SparCC_SW_day10 <- net_WTdTDA_SparCC_SW_day10_100$assoMat2 #dTDA

assomat_WT_SparCC_SW_day10_Phaeobacter <- data.frame(assomat_WT_SparCC_SW_day10)[which(dimnames(assomat_WT_SparCC_SW_day10)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_SW_day10_Phaeobacter_clean <- filter_all(assomat_WT_SparCC_SW_day10_Phaeobacter, any_vars(. != 0& . != 1))

assomat_WT_SparCC_SW_day10_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


assomat_dTDA_SparCC_SW_day10_Phaeobacter <- as.data.frame(assomat_dTDA_SparCC_SW_day10)[which(dimnames(assomat_dTDA_SparCC_SW_day10)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_SW_day10_Phaeobacter_clean <- filter_all(assomat_dTDA_SparCC_SW_day10_Phaeobacter, any_vars(. != 0 & . != 1))

assomat_dTDA_SparCC_SW_day10_Phaeobacter_clean%>%mutate(positive=sum(Phaeobacter_inhibens>0),negative=sum(Phaeobacter_inhibens<0))


# Clean up env -----
#Clean up environment and save all comparisons

# rm(list=ls()) # remove everything from workspace
# tmp.env <- new.env() # create a temporary environment
# load("NetCoMi_netcomp_WTdTDA_SparCC_SW_day10_100_0.6_1000L_161021.RData", envir=tmp.env) # load workspace into temporary environment
# netcomp_WTdTDA_SparCC_Biofilm_day1_100_1000L_0.6 <- get("netcomp_WTdTDA_SparCC_Biofilm_day1_100_1000L_0.6", pos=tmp.env) # get the objects you need into your globalenv()
# netcomp_WTdTDA_SparCC_SW_day1_100_1000L_0.6 <- get("netcomp_WTdTDA_SparCC_SW_day1_100_1000L_0.6", pos=tmp.env) # get the objects you need into your globalenv()
# netcomp_WTdTDA_SparCC_Biofilm_day4_100_1000L_0.6 <- get("netcomp_WTdTDA_SparCC_Biofilm_day4_100_1000L_0.6", pos=tmp.env) # get the objects you need into your globalenv()
# netcomp_WTdTDA_SparCC_SW_day4_100_1000L_0.6 <- get("netcomp_WTdTDA_SparCC_SW_day4_100_1000L", pos=tmp.env) # get the objects you need into your globalenv()
# netcomp_WTdTDA_SparCC_Biofilm_day10_100_1000L_0.6 <- get("netcomp_WTdTDA_SparCC_Biofilm_day10_100_1000L_0.6", pos=tmp.env) # get the objects you need into your globalenv()
# netcomp_WTdTDA_SparCC_SW_day10_100_1000L_0.6 <- get("netcomp_WTdTDA_SparCC_SW_day10_100_1000L_0.6", pos=tmp.env) # get the objects you need into your globalenv()
# 
# #x <- tmp.env$x # equivalent to previous line
# rm(tmp.env) # remove the temporary environment to free up memory
# 
# save.image(file='netcomp_WTdTDA_SparCC_ALL_COMP_100_1000L_221021.RData')

### --- Summary comparisons at all times ----
#Day 1
sink("netcomp_WTdTDA_SparCC_all.txt")
summary(netcomp_WTdTDA_SparCC_Biofilm_day1_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 1000L  )
summary(netcomp_WTdTDA_SparCC_SW_day1_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 1000L  )
#Day 4
summary(netcomp_WTdTDA_SparCC_Biofilm_day4_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 1000L  )
summary(netcomp_WTdTDA_SparCC_SW_day4_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 1000L  )

#Day 10
summary(netcomp_WTdTDA_SparCC_Biofilm_day10_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 1000L  )
summary(netcomp_WTdTDA_SparCC_SW_day10_100_1000L_0.6,  groupNames = c("WT","dTDA"), numbNodes = 1000L  )
sink()


#Strongest positive and negative correlations Phaeobacter -----
#Day 1 Biofilm
assomat_net_WT_SparCC_Biofilm_day1_100 <- net_WTdTDA_SparCC_Biofilm_day1_100$assoMat1 #WT
assomat_net_dTDA_SparCC_Biofilm_day1_100 <- net_WTdTDA_SparCC_Biofilm_day1_100$assoMat2 #dTDA

assomat_net_WT_SparCC_Biofilm_day1_100_Phaeobacter <- data.frame(assomat_net_WT_SparCC_Biofilm_day1_100)[which(dimnames(assomat_net_WT_SparCC_Biofilm_day1_100)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean <- filter_all(assomat_net_WT_SparCC_Biofilm_day1_100_Phaeobacter, any_vars(. != 0& . != 1))

assomat_net_dTDA_SparCC_Biofilm_day1_100_Phaeobacter <- data.frame(assomat_net_dTDA_SparCC_Biofilm_day1_100)[which(dimnames(assomat_net_dTDA_SparCC_Biofilm_day1_100)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean <- filter_all(assomat_net_dTDA_SparCC_Biofilm_day1_100_Phaeobacter, any_vars(. != 0& . != 1))


#Day 1 SW
assomat_net_WT_SparCC_SW_day1_100 <- net_WTdTDA_SparCC_SW_day1_100$assoMat1 #WT
assomat_net_dTDA_SparCC_SW_day1_100 <- net_WTdTDA_SparCC_SW_day1_100$assoMat2 #dTDA

assomat_net_WT_SparCC_SW_day1_100_Phaeobacter <- data.frame(assomat_net_WT_SparCC_SW_day1_100)[which(dimnames(assomat_net_WT_SparCC_SW_day1_100)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean <- filter_all(assomat_net_WT_SparCC_SW_day1_100_Phaeobacter, any_vars(. != 0& . != 1))

assomat_net_dTDA_SparCC_SW_day1_100_Phaeobacter <- data.frame(assomat_net_dTDA_SparCC_SW_day1_100)[which(dimnames(assomat_net_dTDA_SparCC_SW_day1_100)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean <- filter_all(assomat_net_dTDA_SparCC_SW_day1_100_Phaeobacter, any_vars(. != 0& . != 1))


#Day 4 Biofilm
assomat_net_WT_SparCC_Biofilm_day4_100 <- net_WTdTDA_SparCC_Biofilm_day4_100$assoMat1 #WT
assomat_net_dTDA_SparCC_Biofilm_day4_100 <- net_WTdTDA_SparCC_Biofilm_day4_100$assoMat2 #dTDA

assomat_net_WT_SparCC_Biofilm_day4_100_Phaeobacter <- data.frame(assomat_net_WT_SparCC_Biofilm_day4_100)[which(dimnames(assomat_net_WT_SparCC_Biofilm_day4_100)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean <- filter_all(assomat_net_WT_SparCC_Biofilm_day4_100_Phaeobacter, any_vars(. != 0& . != 1))

assomat_net_dTDA_SparCC_Biofilm_day4_100_Phaeobacter <- data.frame(assomat_net_dTDA_SparCC_Biofilm_day4_100)[which(dimnames(assomat_net_dTDA_SparCC_Biofilm_day4_100)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean <- filter_all(assomat_net_dTDA_SparCC_Biofilm_day4_100_Phaeobacter, any_vars(. != 0& . != 1))


#Day 4 SW
assomat_net_WT_SparCC_SW_day4_100 <- net_WTdTDA_SparCC_SW_day4_100$assoMat1 #WT
assomat_net_dTDA_SparCC_SW_day4_100 <- net_WTdTDA_SparCC_SW_day4_100$assoMat2 #dTDA

assomat_net_WT_SparCC_SW_day4_100_Phaeobacter <- data.frame(assomat_net_WT_SparCC_SW_day4_100)[which(dimnames(assomat_net_WT_SparCC_SW_day4_100)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean <- filter_all(assomat_net_WT_SparCC_SW_day4_100_Phaeobacter, any_vars(. != 0& . != 1))

assomat_net_dTDA_SparCC_SW_day4_100_Phaeobacter <- data.frame(assomat_net_dTDA_SparCC_SW_day4_100)[which(dimnames(assomat_net_dTDA_SparCC_SW_day4_100)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean <- filter_all(assomat_net_dTDA_SparCC_SW_day4_100_Phaeobacter, any_vars(. != 0& . != 1))

#Day 10 Biofilm
assomat_net_WT_SparCC_Biofilm_day10_100 <- net_WTdTDA_SparCC_Biofilm_day10_100$assoMat1 #WT
assomat_net_dTDA_SparCC_Biofilm_day10_100 <- net_WTdTDA_SparCC_Biofilm_day10_100$assoMat2 #dTDA

assomat_net_WT_SparCC_Biofilm_day10_100_Phaeobacter <- data.frame(assomat_net_WT_SparCC_Biofilm_day10_100)[which(dimnames(assomat_net_WT_SparCC_Biofilm_day10_100)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean <- filter_all(assomat_net_WT_SparCC_Biofilm_day10_100_Phaeobacter, any_vars(. != 0& . != 1))

assomat_net_dTDA_SparCC_Biofilm_day10_100_Phaeobacter <- data.frame(assomat_net_dTDA_SparCC_Biofilm_day10_100)[which(dimnames(assomat_net_dTDA_SparCC_Biofilm_day10_100)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean <- filter_all(assomat_net_dTDA_SparCC_Biofilm_day10_100_Phaeobacter, any_vars(. != 0& . != 1))




#Day 10 SW
assomat_net_WT_SparCC_SW_day10_100 <- net_WTdTDA_SparCC_SW_day10_100$assoMat1 #WT
assomat_net_dTDA_SparCC_SW_day10_100 <- net_WTdTDA_SparCC_SW_day10_100$assoMat2 #dTDA

assomat_net_WT_SparCC_SW_day10_100_Phaeobacter <- data.frame(assomat_net_WT_SparCC_SW_day10_100)[which(dimnames(assomat_net_WT_SparCC_SW_day10_100)[[1]] == "Phaeobacter_inhibens")]
assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean <- filter_all(assomat_net_WT_SparCC_SW_day10_100_Phaeobacter, any_vars(. != 0& . != 1))

assomat_net_dTDA_SparCC_SW_day10_100_Phaeobacter <- data.frame(assomat_net_dTDA_SparCC_SW_day10_100)[which(dimnames(assomat_net_dTDA_SparCC_SW_day10_100)[[1]] == "Phaeobacter_inhibens")]
assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean <- filter_all(assomat_net_dTDA_SparCC_SW_day10_100_Phaeobacter, any_vars(. != 0& . != 1))


#Generating big dataframe with eigenvector values of each node Phaeobacter is connected to


assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean$Day <- rep("1")
assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean$Environment <- rep("Biofilm")
assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean$Treatment <- rep("WT")
assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean <- merge(assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6$centralities$eigenv1), by = 0)
colnames(assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")

assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean$Day <- rep("1")
assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean$Environment <- rep("Biofilm")
assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean$Treatment <- rep("dTDA")
assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean <- merge(assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_Biofilm_day1_100_0.6$centralities$eigenv2), by = 0)
colnames(assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")

assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean$Day <- rep("1")
assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean$Environment <- rep("Planktonic suspension")
assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean$Treatment <- rep("WT")
assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean <- merge(assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_SW_day1_100_0.6$centralities$eigenv1), by = 0)
colnames(assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")

assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean$Day <- rep("1")
assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean$Environment <- rep("Planktonic suspension")
assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean$Treatment <- rep("dTDA")
assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean <- merge(assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_SW_day1_100_0.6$centralities$eigenv2), by = 0)
colnames(assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean$Day <- rep("4")
assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean$Environment <- rep("Biofilm")
assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean$Treatment <- rep("WT")
assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean <- merge(assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6$centralities$eigenv1), by = 0)
colnames(assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")

#Sum of positive and negatives
assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean %>% filter(Degree>0) %>% dplyr::summarise(n())
assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean %>% filter(Degree<0) %>% dplyr::summarise(n())


assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean$Day <- rep("4")
assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean$Environment <- rep("Biofilm")
assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean$Treatment <- rep("dTDA")
assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean <- merge(assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_Biofilm_day4_100_0.6$centralities$eigenv2), by = 0)
colnames(assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean$Day <- rep("4")
assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean$Environment <- rep("Planktonic suspension")
assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean$Treatment <- rep("WT")
assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean <- merge(assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_SW_day4_100$centralities$eigenv1), by = 0)
colnames(assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")

assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean$Day <- rep("4")
assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean$Environment <- rep("Planktonic suspension")
assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean$Treatment <- rep("dTDA")
assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean <- merge(assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_SW_day4_100$centralities$eigenv2), by = 0)
colnames(assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean$Day <- rep("10")
assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean$Environment <- rep("Biofilm")
assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean$Treatment <- rep("WT")
assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean <- merge(assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6$centralities$eigenv1), by = 0)
colnames(assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean$Day <- rep("10")
assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean$Environment <- rep("Biofilm")
assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean$Treatment <- rep("dTDA")
assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean <- merge(assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_Biofilm_day10_100_0.6$centralities$eigenv2), by = 0)
colnames(assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean$Day <- rep("10")
assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean$Environment <- rep("Planktonic suspension")
assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean$Treatment <- rep("WT")
assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean <- merge(assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_SW_day10_100_0.6$centralities$eigenv1), by = 0)
colnames(assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean$Day <- rep("10")
assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean$Environment <- rep("Planktonic suspension")
assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean$Treatment <- rep("dTDA")
assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean <- merge(assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean,as.data.frame(props_net_WTdTDA_SparCC_SW_day10_100_0.6$centralities$eigenv2), by = 0)
colnames(assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean) <- c("Node","Degree","Day","Environment", "Treatment","Eigenvector")


degree_Phaeobacter <- rbind(assomat_WT_SparCC_Biofilm_day1_100_Phaeobacter_clean, assomat_dTDA_SparCC_Biofilm_day1_100_Phaeobacter_clean,
                            assomat_WT_SparCC_SW_day1_100_Phaeobacter_clean,assomat_dTDA_SparCC_SW_day1_100_Phaeobacter_clean, 
                            assomat_WT_SparCC_Biofilm_day4_100_Phaeobacter_clean, assomat_dTDA_SparCC_Biofilm_day4_100_Phaeobacter_clean,
                            assomat_WT_SparCC_SW_day4_100_Phaeobacter_clean, assomat_dTDA_SparCC_SW_day4_100_Phaeobacter_clean,
                            assomat_WT_SparCC_Biofilm_day10_100_Phaeobacter_clean, assomat_dTDA_SparCC_Biofilm_day10_100_Phaeobacter_clean,
                            assomat_WT_SparCC_SW_day10_100_Phaeobacter_clean, assomat_dTDA_SparCC_SW_day10_100_Phaeobacter_clean)

colnames(degree_Phaeobacter) <- c( "Genus","Correlation strength","Day","Environment","Treatment","Eigenvector score")
degree_Phaeobacter_tax <- merge(degree_Phaeobacter, as.data.frame(tax_table(phylo_main_no_cont_NoOutliers_genus)), by = "Genus") 



#Ordering the Classes
#degree_Phaeobacter_tax$Class <- factor(degree_Phaeobacter_tax$Class , levels=c("Actinobacteria", "Alphaproteobacteria", "Betaproteobacteria", "Chloroplast", "Cyanobacteria", "Deltaproteobacteria","Epsilonproteobacteria","Gammaproteobacteria","Flavobacteriia","Planctomycetia","Sphingobacteriia","Verrucomicrobiae"))
degree_Phaeobacter_tax$Treatment <- factor(degree_Phaeobacter_tax$Treatment, levels=c("WT", "dTDA"))

degree_Phaeobacter_tax_Biofilm_Day4 <- degree_Phaeobacter_tax %>% filter(Day=="4", Environment=="Biofilm")


degree_Phaeobacter_tax_Biofilm_Day4 <- degree_Phaeobacter_tax %>% filter(Day=="4", Environment=="Biofilm")


#Join LogFC size from ANCOM analysis and plot as size
load(file="df_fig1_Biofilm_day4")
df_fig1_Biofilm_day4
#Load df_fig1_Biofilm_day4 from ANCOM-BC analysis

df_fig1_Biofilm_day4_WT.vs.dTDA <- df_fig1_Biofilm_day4[,c(1,3)] %>% filter(taxon_id %in% degree_Phaeobacter_tax_Biofilm_Day4$Genus) %>% filter(dTDA!=0) 
colnames(df_fig1_Biofilm_day4_WT.vs.dTDA) <- c("Genus","LogFC_WTvsdTDA")
df_fig1_Biofilm_day4_WT.vs.dTDA$Treatment <- rep("WT")
degree_Phaeobacter_tax_Biofilm_Day4 <- left_join(degree_Phaeobacter_tax_Biofilm_Day4, df_fig1_Biofilm_day4_WT.vs.dTDA, by=c("Genus","Treatment"))

degree_Phaeobacter_tax_Biofilm_Day4$LogFC_WTvsdTDA <- abs(degree_Phaeobacter_tax_Biofilm_Day4$LogFC_WTvsdTDA)
#Change NAs to value 0
degree_Phaeobacter_tax_Biofilm_Day4[is.na(degree_Phaeobacter_tax_Biofilm_Day4)] <- 0


####
df_fig1_Biofilm_day4

#Filtering all taxa wiuth read sum < 100
phylo_main_no_cont_NoOutliers_Genus_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_genus, function (x) {sum(x > 100) > 0}, prune=TRUE)

#summary(colSums(otu_table(phylo_main_no_cont_NoOutliers_Genus_filter)))

#Remove Phaeobacter before making differential analyses
phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo <- subset_taxa(phylo_main_no_cont_NoOutliers_Genus_filter, Genus != "Phaeobacter_inhibens")


NoOutliers_Genus_filter_noPhaeo_tax_dat <- as.data.frame(tax_table(phylo_main_no_cont_NoOutliers_Genus_filter_noPhaeo))
NoOutliers_Genus_filter_noPhaeo_tax_dat$taxon_id <- NoOutliers_Genus_filter_noPhaeo_tax_dat$unique

#Sort by drecreasing order
df_fig1_Biofilm_day4_Order_Control.vs.dTDA <- df_fig1_Biofilm_day4 %>% left_join(NoOutliers_Genus_filter_noPhaeo_tax_dat, by = "taxon_id") 

#How many of 33 positive connections was confermed by ANCOM-BC
degree_Phaeobacter_tax_Biofilm_Day4 %>% filter(LogFC_WTvsdTDA>0)
#How many of 22 negative connections was confermed by ANCOM-BC
degree_Phaeobacter_tax_Biofilm_Day4 %>% filter(LogFC_WTvsdTDA<0)


#Sort Class in same order as network plot with the same colors 
# degree_Phaeobacter_tax_Biofilm_Day4$Class <- ifelse(is.na(degree_Phaeobacter_tax_Biofilm_Day4$Class), "Deltaproteobacteria", as.character(degree_Phaeobacter_tax_Biofilm_Day4$Class))
unique((degree_Phaeobacter_tax_Biofilm_Day4 %>% filter(Day=="4", Environment=="Biofilm"))$Class)
levels(degree_Phaeobacter_tax$Class)
degree_Phaeobacter_tax_Biofilm_Day4$Class <- factor(degree_Phaeobacter_tax_Biofilm_Day4$Class, levels=c("Flavobacteriia","Gammaproteobacteria",
                                                                                "Cytophagia","Alphaproteobacteria",
                                                                                "Actinobacteria", "Sphingobacteriia",
                                                                                "Caldilineae","Deltaproteobacteria",
                                                                                "Opitutae", "Planctomycetia",
                                                                                "Deinococci"))



# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia", 
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia", 
#                                                 "Deinococci",
#                                               "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria", 
#                                                 "Verrucomicrobiae"))


degree_Phaeobacter_tax_plot <- degree_Phaeobacter_tax_Biofilm_Day4 %>%
    ggplot(aes(x=`Correlation strength`, y=`Eigenvector score`)) +
    geom_point(aes(size=LogFC_WTvsdTDA), color="black",pch=21) +
     geom_point(aes(colour=Class, size=LogFC_WTvsdTDA), alpha=0.7) + scale_size(range = c(3, 14)) +
    scale_colour_manual(values = Paired_colors_expand_norm) +
    facet_grid(cols=vars(factor(Treatment, levels=c("WT", "dTDA"))))  + theme_bw() +
    theme(#axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  legend.position = "none") + geom_vline(xintercept = 0, colour="gray", size = 0.3)


degree_Phaeobacter_tax_plot_lg <- degree_Phaeobacter_tax_Biofilm_Day4 %>%
    ggplot(aes(x=`Correlation strength`, y=`Eigenvector score`)) +
    geom_point(aes(size=LogFC_WTvsdTDA), color="black",pch=21,)  +
    facet_grid(cols=vars(factor(Treatment, levels=c("WT", "dTDA"))))  + theme_bw() +
    theme(#axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  legend.position = "none") 


legend_degree_Phaeobacter_tax_plot <- cowplot::get_legend(degree_Phaeobacter_tax_plot + theme(legend.position="bottom") + 
                                                              guides(colour = guide_legend(title.position = "top", ncol = 2), size = guide_legend(title.position = "top", nrow = 1)))


library(gridExtra)
library(magick)


#ggdraw(degree_Phaeobacter_tax_plot) + draw_image("Figures/Figure_3A_Network.tiff", scale = 0.9)

#FIGURE 3 -----

degree_Phaeobacter_tax_plot_B <- plot_grid(NULL,NULL,
                         degree_Phaeobacter_tax_plot, NULL, 
                         rel_widths = c(2,1.5), rel_heights = c(2,1))




degree_Phaeobacter_tax_plot_A_B <- ggdraw(degree_Phaeobacter_tax_plot_B)  + 
    draw_image("Figures/Figure_3A_Network.tiff", 
               y = 0.16) + draw_label("B", x = 0.025, y = 0.35, fontface = "bold") + draw_label("A", x = 0.025, y = 0.97, fontface = "bold") 
    

   
tiff("Figures/Figure_3.tiff", units="in", width=8.5, height=8, res=300)

ggdraw(degree_Phaeobacter_tax_plot_B)  + 
    draw_image("Figures/Figure_3A_Network.tiff", 
               y = 0.16) + draw_label("B", x = 0.025, y = 0.35, fontface = "bold") + draw_label("A", x = 0.025, y = 0.97, fontface = "bold") 
    
dev.off() 

tiff("Figures/Figure_3.tiff", units="in", width=8.5, height=8, res=300)

ggdraw(degree_Phaeobacter_tax_plot_A_B)  + 
    draw_image("Figures/Figure_3A_Network_legends.tiff", 
               x = 0.25, y=-0.28, scale = 0.43) 

dev.off() 



  # pdf(file = "degree_Phaeobacter_tax.pdf",   # The directory you want to save the file in
#     width = 7.5, # The width of the plot in inches
#     height = 5) # The height of the plot in inches
# 
# grid.arrange(degree_Phaeobacter_tax_plot, legend_degree_Phaeobacter_tax_plot, ncol=1, nrow = 2, 
#                               layout_matrix = rbind(c(1), c(2)),
#                                 heights = c(1.6, 0.7))
# dev.off()


# png(filename = "degree_Phaeobacter_tax.png",
#     width = 6, height = 4, units = "px", pointsize = 12,
#     bg = "white")
# 
# grid.arrange(degree_Phaeobacter_tax_plot, legend_degree_Phaeobacter_tax_plot, ncol=1, nrow = 2, 
#              layout_matrix = rbind(c(1), c(2)),
#              heights = c(1.6, 0.7))
# dev.off()





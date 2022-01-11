setwd("/Users/pernillekjersgaardbech/Documents/QIIME2/PChem0002/MACem00002/phyloseq/")
library(rRDPData)
library(debar)
library(tidyverse)
library("devtools")
library(Rcpp)
library("dada2")
library(plyr)
library("Biostrings")
library("phyloseq")
library(ggplot2)
library( "DESeq2" )
library(vegan)
library(venn)


## ASV Tables Created in QIIME2 ----
asv <- read.table("ASV_table.txt", row.names = 1, header = TRUE, sep = "\t")

#read average 
mean(colSums(asv))
sd(colSums(asv))
#read average
load('TrainedRdp021020_qiime.RData')
# str(colSums(asv))
# asv_200000 <- cols_timeums(asv/t(cols_timeums(asv))*200000)
# str(asv/cols_timeums(asv)*200000)
#raremax <- (raremax <- min(rowSums(asv_rarefy)))
#rarecurve(t(asv), step = 10000, sample = 20000, col = "blue", cex = 0.6, label=F)

ddenoised_seq_qiime <- Biostrings::readDNAStringSet("dna-sequences.fasta", format = "fasta")

asv.table <- phyloseq::otu_table(asv, taxa_are_rows = TRUE)
phylo_main <- phyloseq::phyloseq(asv.table, denoised_seq_qiime)
#asv/asv table

#refseq(phylo_main)
#asv_table(phylo_main)
metadata <- data.frame(t(read.table("metadata.txt", row.names = 1, header = TRUE, sep = "\t")))
str(metadata)
#str(asv_table(phylo_main))
#Inspect distribution of sequence lengths
#****plot(table(nchar(getSequences(denoised_seq_qiime))))
#table(nchar(getSequences(denoised_seq_qiime)))

#Classifier rdp through dada2 in R on Sif to save time
#library(rRDPData)
#library("devtools")
#library(Rcpp)
#library("dada2")
#library(plyr)
#library("Biostrings")
#denoised_seq_qiime = readDNAStringSet("dna-sequences.fasta")
#tt_qiime_both_strands <- assignTaxonomy(denoised_seq_qiime, "silva_nr_v132_train_set.fa", multithread = 20, tryRC =TRUE)
#tt.plus_qiime_both_strands <- addSpecies(tt_qiime_both_strands, "silva_species_assignment_v132.fa", verbose=TRUE)

#set.seed(100) # Initialize random number generator for reproducibility
#taxa_genus <- assignTaxonomy(denoised_seq_qiime, "rdp_train_set_16.fa", multithread=20, tryRC =TRUE)
#genus.species <- assignSpecies(denoised_seq_qiime, "rdp_species_assignment_16.fa", tryRC =TRUE)
#genus.species_allowMultiple <- assignSpecies(denoised_seq_qiime, "rdp_species_assignment_16.fa", tryRC =TRUE, allowMultiple=TRUE)
#save.image(file='TrainedRdp021020_qiime.RData')
#genus.species.plus <- addSpecies(genus.species, "rdp_species_assignment_16.fa", verbose=TRUE)
#603 out of 39917 were assigned to the species level.
#Of which 603 had genera consistent with the input table
#By default the assignSpecies method only returns species assignments if there is no ambiguity, i.e. all exact matches were to the same species. However, given that we are generally working with fragments of the 16S gene, it is common that exact matches are made to multiple sequences that are identical over the sequenced region.
#genus.species_allowMultiple.plus <- addSpecies(genus.species_allowMultiple, "rdp_species_assignment_16.fa", verbose=TRUE)
#603 out of 39917 were assigned to the species level.
#Of which 966 had genera consistent with the input table.
#save.image(file='TrainedRdp021020_qiime.RData')
#taxa_genus.plus <- addSpecies(taxa_genus, "rdp_species_assignment_16.fa", verbose=TRUE)
#save.image(file='TrainedRdp021020_qiime.RData')

#load('TrainedRdp021020_qiime.RData')
# 
# head(unname(taxa_genus.plus))
# 
taxa.print <- taxa_genus.plus # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# 
rownames(taxa_genus.plus) <- rownames(otu_table(phylo_main)) #naming rownames same of asv_table


#decomtam ----
library(decontam)
packageVersion("decontam") 

colnames(asv) # our blank are column 46
vector_for_decontam <- c(rep(FALSE, c(45)), rep(TRUE, c(1)), rep(FALSE, c(123)))

contam_df <- isContaminant(t(asv.table), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 31 as contaminants

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

# taking a peak at the two identified contaminants
taxa_genus.plus[row.names(taxa_genus.plus) %in% contam_asvs, ]

# making new fasta file
#contam_indices <- which(asv %in% paste0(">", contam_asvs))
#dont_want <- sort(c(contam_indices, contam_indices + 1))
#asv_fasta_no_contam <- asv[- dont_want]

# making new count table
asv_tab_no_contam <- asv[!row.names(asv.table) %in% contam_asvs, ]
asv_tab_no_contam <- phyloseq::otu_table(asv_tab_no_contam, taxa_are_rows = TRUE)

# making new taxonomy table
asv_tax_no_contam <- taxa_genus.plus[!row.names(taxa_genus.plus) %in% contam_asvs, ]
#How many ASvs were classified
colSums(!is.na(as(asv_tax_no_contam, "matrix")))/(39330+556)*100
#replace NA values with Unclassified
asv_tax_no_contam[is.na(asv_tax_no_contam)] <- "Unclassified"


#No contam Phylo-object -----
str(metadata)
metadata_no_neg <- metadata[-c(46), ]
phylo_main_no_cont <- phyloseq(otu_table(asv_tab_no_contam), 
                               tax_table(asv_tax_no_contam), sample_data(metadata_no_neg))


#Remove outliers
#Remove outliers dTDA_3SurT1_1 and dTDA_2SurT4_1 and WT(only P. inhibens ASVs) and X_1WatT0_2 replicates since they were amplified by 38 cycles

phylo_main_no_cont_no_outliers <- subset_samples(phylo_main_no_cont, Name!="dTDA_3SurT1_1" & Name!= "dTDA_2SurT4_1" & Name!= "WT_1SurT0" & Name!= "WT_2WatT4_2" & Name!= "_1WatT0_2_1" & Name!= "_1WatT0_2_2"& Name!= "_1WatT0_2_3")
#Phaeobacter ASVs
phylo_main_no_cont_PhaeoASVs <- subset_samples(phylo_main_no_cont, Name== "WT_1SurT0")
#Remove rowsums = 0
phylo_main_no_cont_PhaeoASVs <- otu_table(phylo_main_no_cont_PhaeoASVs)[rowSums(otu_table(phylo_main_no_cont_PhaeoASVs))!=0,]
#Phaeobacter ASV names
PhaeoASV_names <- as.array(c("7933f9d855a75dc3b415409eb5125d5c","f989b73965cf9a81e2e4f8b87d223dc6"))

#Remove rowsums = 0
ASV_table_clean_NoOutliers <- otu_table(phylo_main_no_cont_no_outliers)[rowSums(otu_table(phylo_main_no_cont_no_outliers))!=0,]

#make clean tax table
ASV_table_clean_NoOutliers_tax <- (asv_tax_no_contam[rownames(asv_tax_no_contam) %in% rownames(ASV_table_clean_NoOutliers), ])
ASV_table_clean_NoOutliers_tax_dat <- as.data.frame(ASV_table_clean_NoOutliers_tax)
ASV_table_clean_NoOutliers_tax_dat$ASV_id <- rownames(ASV_table_clean_NoOutliers_tax)
ASV_table_clean_NoOutliers_tax <- (as.matrix(ASV_table_clean_NoOutliers_tax_dat))
dim(asv_tax_no_contam)
str(ASV_table_clean_NoOutliers_tax)
str(ASV_table_clean_NoOutliers_tax_dat)
#Change wrong classification of Phaeobacter (Tropicibacter)
Phaeo_tax_indx <- which(ASV_table_clean_NoOutliers_tax_dat$ASV_id  %in%  PhaeoASV_names)
ASV_table_clean_NoOutliers_tax_dat[Phaeo_tax_indx,]

ASV_table_clean_NoOutliers_tax_dat$Genus[Phaeo_tax_indx] = rep("Phaeobacter_inhibens")
ASV_table_clean_NoOutliers_tax <- (as.matrix(ASV_table_clean_NoOutliers_tax_dat))

#Rarecurves ----
#(raremax <- min(rowSums(otu)))
phylo_main_no_cont_no_outliers_dat <- as.data.frame((as(t(otu_table(phylo_main_no_cont_no_outliers)), "matrix")))

#rarecurve(phylo_main_no_cont_no_outliers_dat, step = 1000, col = "blue", cex = 0.6, label=FALSE, xlim=c(0, 25000), ylim=c(0, 3000))


Group_ind_metadata <- rownames(phylo_main_no_cont_no_outliers_dat) #individual sample names
Treatment <- c(rep("Control", 6), rep("dTDA", 1), rep("Control", 2), rep("dTDA", 5), 
                             rep("Control", 1), rep("dTDA", 2), rep("WT", 9), rep("Control", 9), rep("dTDA", 8), 
                             rep("WT", 9), rep("Seawater at Day 0", 3), rep("Control", 9), rep("dTDA", 9), rep("WT", 9), rep("Control", 9),
                             rep("dTDA", 9), rep("WT", 9), rep("Control", 9), rep("dTDA", 9), rep("WT", 8), 
                             rep("Control", 9), rep("dTDA", 8), rep("WT", 9))


Environment <- c(rep("Biofilm", 52), rep("Seawater at Day 0", 3), rep("Biofilm", 27), rep("Planktonic suspension", 79))

Day <- c(rep("1", 26),rep("4", 26), rep("0", 3), rep("10", 27), rep("1", 27), rep("4", 26), rep("10", 26))
# c(rep("dTDA", 6), "dTDA", rep("dTDA", 2), rep("dTDA", 6), rep("WT", 9))
ASV_metadata_all <- as.data.frame(cbind(Group_ind_metadata, Environment, Treatment, Day))
rownames(ASV_metadata_all) <- ASV_metadata_all$Group_ind_metadata

#Phylo object with no contaminants and no outliers -----
phylo_main_no_cont_NoOutliers <- phyloseq(otu_table(phylo_main_no_cont_no_outliers, taxa_are_rows = TRUE), 
                                          tax_table(ASV_table_clean_NoOutliers_tax), sample_data(ASV_metadata_all))

sample_data(phylo_main_no_cont_NoOutliers)

save.image(file='phylo_main_no_cont_NoOutliers.RData')



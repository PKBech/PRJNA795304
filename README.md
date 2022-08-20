# Data and scripts to reproduce results from the manuscript with data in BioProject PRJNA795304

### Data
All data used in the study can be found in the folder: 'Data/'.

Two files for the creation of the phyloseq object in the R script "phylo_main_no_cont_NoOutliers.R" 
- ASV_table.txt created by the QIIME2 and DADA2 pipeline
- metadata.txt all metadata explaining each sample

Statistical analysis of the qPCR of the absolutte population and Phaeobacter spp. abundances
- 16S_absolutte_qPCR_2.txt
- qPCR_040521_50ulelute.txt

### Scripts
All scripts used for analysis can be found in 'Rscripts/'. 

All downstream analysis loads phyloseq created in: 
- phylo_main_no_cont_NoOutliers.R

The R environment "phylo_main_no_cont_NoOutliers.RData" is loaded into following scripts:
- alpha_161121.R (Table 1; Figure S1)
- Beta-diversity_030122.R (Table 1 and 2; Figure 2)
- ANCOMBC_040422.R (Figure 3)
- NetCoMi_030122.R (Table S1; Figure 4)

To determine the low abundant read cut-off we used:
- Filtration210921.R

Absolutte population and Phaeobacter spp. abundances
- qPCR_040422.R (Table 1; Figure 1)

![plot](./Figures/Figure_1.tiff)




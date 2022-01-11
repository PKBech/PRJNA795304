#Filtration 

#Aggregate to genus level
phylo_main_no_cont_NoOutliers_Genus <- aggregate_taxa(phylo_main_no_cont_NoOutliers, "Genus")
#read sum per genus 
genus_sum <- rowSums(data.frame(otu_table(phylo_main_no_cont_NoOutliers_Genus)))
#data frame
genus_sum_dat=data.frame(N=names(genus_sum), val=genus_sum)
#Ordering in decreasing order
genus_sum_dat=genus_sum_dat[order(genus_sum_dat$val,decreasing = T),]

#Plot log10 of all genus sums
plot(log10(genus_sum_dat$val))

#Plot differential of log10 for all genus sums
plot(diff(diff(log10(genus_sum_dat$val))))

length(genus_sum_dat$val)

100*sum(genus_sum_dat$val[genus_sum_dat$val<100])/sum(genus_sum_dat$val) #0.02948133 is filtered out

#read average 
mean(colSums(asv))
sd(colSums(asv))



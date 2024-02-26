rm(list = ls())
library(mambo)

df_16s <- read.csv("../../Data/merged_data/no_duplicates/Merged2018_16S_otu_filtered.csv", row.names=1)
df_18s <- read.csv("../../Data/merged_data/no_duplicates/Merged2018_18S_otu_filtered.csv", row.names=1)

pca.test <- mambo('18S', df_18s, '16S', df_16s, nrep = 20, chains = 5)

plotPCs(pca.test, '18S')
plotPCs(pca.test, '18S', 3, 4)

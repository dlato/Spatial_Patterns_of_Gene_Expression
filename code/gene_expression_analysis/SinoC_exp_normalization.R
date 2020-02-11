###########################
#SINO CHROM 
###########################
library(ggplot2)
require(ggplot2)
require("ggplot2")
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)
library("edgeR")




chrom_datafile <- "GSE69880_CHROM_normalized_reads_counts_abcdef_GEO_subset.txt"
data_chrom <- read.delim(chrom_datafile, header = TRUE, row.names = "Gene_ID")
data_frame_exp <- as.data.frame(data_chrom)

#all samples should be apart of the same group
group <- c(1,1,1)

#normalize all replicates and samples
exp_list <- DGEList(counts=data_frame_exp, group=group)
exp_normalized <- calcNormFactors(exp_list)
normalized_exp_df <- cpm(exp_normalized, normalized.lib.sizes=FALSE)

write.csv(normalized_exp_df, file = "gene_expression_normalized_gene_name.txt", row.names = FALSE)

#re-import the gene exp data to add the gene names to it.
exp_dat <- "gene_expression_normalized_gene_name.txt"
exp_df <- read.csv(exp_dat)

#combining replicates by taking the average
exp_df$Median <- apply(exp_df[,c(1,2,3)], 1, median.default)
head(exp_df)
#med_exp_graph <- ggplot(mean_df_melt, aes(x = gene_name, y = exp_var, color = sample_num)) + 
#  ggtitle(label = "Streptomyces") +
#  #  geom_point(position = position_dodge(width = 0.4)) +
#  geom_point(position = position_dodge(width = 0.4)) 
##  theme(legend.position = "none")
#pdf("streo_median_mean_gene_exp_check.pdf")
#med_exp_graph
#dev.off()


median_exp_df <- exp_df[,4, drop=FALSE]

geneID <- row.names(normalized_exp_df)
exp_dat_ID <- cbind(geneID, median_exp_df)

  

write.csv(exp_dat_ID, file = "Sino_chrom_gene_expression_combined_normalized.txt", row.names = FALSE)

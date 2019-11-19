###########################
#SINO pSymA
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


chrom_datafile_1 <- "GSE69880_PSYMA_normalized_reads_counts_abcdef_GEO_subset.txt"
data_chrom_1 <- read.delim(chrom_datafile_1, header = TRUE, row.names = "Gene_ID")
exp_dat_all_reps = as.data.frame(data_chrom_1)

#chrom_datafile_1 <- "NC_020527_10.count"
#data_chrom_1 <- read.delim(chrom_datafile_1, header = TRUE)
#rownames(data_chrom_1) <- data_chrom_1$gene
#sub_1 <- data_chrom_1[,2, drop=FALSE]
#
#chrom_datafile_2 <- "NC_020527_5.count"
#data_chrom_2 <- read.delim(chrom_datafile_2, header = TRUE)
#rownames(data_chrom_2) <- data_chrom_2$gene
#sub_2 <- data_chrom_2[,2, drop=FALSE]
#
##merging each dataset by row name into one large dataset
#exp_dat_all_reps <- Reduce(merge, lapply(list(sub_1, sub_2), function(x) data.frame(x, rn = row.names(x))))
#row.names(exp_dat_all_reps) <- exp_dat_all_reps$rn
#exp_dat_all_reps <- exp_dat_all_reps[,-1]

#all samples should be apart of the same group
num_samples <- length(colnames(exp_dat_all_reps))
group <- rep(1, num_samples)

#normalize all replicates and samples
exp_list <- DGEList(counts=exp_dat_all_reps, group=group)
exp_normalized <- calcNormFactors(exp_list)
normalized_exp_df <- cpm(exp_normalized, normalized.lib.sizes=FALSE)

write.csv(normalized_exp_df, file = "gene_expression_normalized_gene_name.txt", row.names = FALSE)

#re-import the gene exp data to add the gene names to it.
exp_dat <- "gene_expression_normalized_gene_name.txt"
exp_df <- read.csv(exp_dat)

#combining replicates by taking the average
exp_df$Median <- apply(exp_df[,c(1,2,3)], 1, median.default)
#below is unecessary bc there is only one dataset that we have.
#mean_df <- exp_df[,c(1,4)]
#head(mean_df)
#
#
##making dataframe with max - mean
#mean_df <- abs(mean_df - rowMeans(mean_df)) 
#head(mean_df)
#geneID <- row.names(normalized_exp_df)
#mean_df <- cbind(mean_df,geneID)
#mean_df_melt <- melt(mean_df, id.vars = "geneID", value.name = "sample_num")
#head(mean_df_melt)
#colnames(mean_df_melt) <- c("gene_name", "sample_num", "exp_var") 
#head(mean_df_melt)
#class(mean_df_melt$gene_name)
#class(mean_df_melt$sample_num)
#class(mean_df_melt$exp_var)
#
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

  

write.csv(exp_dat_ID, file = "Sino_pSymA_gene_expression_combined_normalized.txt", row.names = FALSE)

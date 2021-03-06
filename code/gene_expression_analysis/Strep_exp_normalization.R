###########################
#Strep CHROM 
###########################
library(ggplot2)
require(ggplot2)
require("ggplot2")
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)
library("edgeR")



chrom_datafile_2 <- "GSE57268_GSM1378113_WT_Day_3a_counts.txt"
data_chrom_2 <- read.delim(chrom_datafile_2, header = FALSE)
colnames(data_chrom_2) <- c("Gene_ID", "2_1")
rownames(data_chrom_2) <- data_chrom_2$Gene_ID
sub_2 <- data_chrom_2[,2, drop=FALSE]

chrom_datafile_3 <- "GSE57268_GSM1378114_WT_Day_3b_counts.txt"
data_chrom_3 <- read.delim(chrom_datafile_3, header = FALSE)
colnames(data_chrom_3) <- c("Gene_ID", "2_2")
rownames(data_chrom_3) <- data_chrom_3$Gene_ID
sub_3 <- data_chrom_3[,2, drop=FALSE]

#merging each dataset by row name into one large dataset
exp_dat_all_reps <- Reduce(merge, lapply(list(sub_2, sub_3), function(x) data.frame(x, rn = row.names(x))))
row.names(exp_dat_all_reps) <- exp_dat_all_reps$rn
exp_dat_all_reps <- exp_dat_all_reps[,-1]

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

#combining replicates by finding the median for each group of replicates
exp_df$Sample2_median <- apply(exp_df[,c(1,2)], 1, median.default)

##combining each dataset by finding median of the medians of each dataset
#exp_df$Median <- apply(exp_df[,c(1,4)], 1, median.default)
head(exp_df)

#mytitle <- expression(paste(italic("Streptomyces"), " Chromosome"))
#med_exp_graph <- ggplot(mean_df_melt, aes(x = gene_name, y = log10(exp_var), color = sample_num)) + 
#  ggtitle(mytitle) +
#  xlab("Genes") +
#  ylab("|CPM - Mean CPM| (log10)") +
#  labs(color = "Datasets")+
#  scale_color_hue(labels = c("T999", "T888"))+
#  #  geom_point(position = position_dodge(width = 0.4)) +
#  geom_point(position = position_dodge(width = 0.4),
##             alpha = 1/10,
#             shape = 1) +
#  theme(axis.text.x = element_blank())
#pdf("strep_median_mean_gene_exp_check.pdf")
#med_exp_graph
#dev.off()


median_exp_df <- exp_df[,-c(1:2)]

geneID <- row.names(normalized_exp_df)
exp_dat_ID <- cbind(geneID, median_exp_df)

  

write.csv(exp_dat_ID, file = "Strep_chrom_gene_expression_combined_normalized.txt", row.names = FALSE)

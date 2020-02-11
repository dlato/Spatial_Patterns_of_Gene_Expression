###########################
#BASS CHROM 
###########################
library(ggplot2)
require(ggplot2)
require("ggplot2")
library(reshape2)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library("edgeR")




chrom_datafile_1 <- "GSE104816_GSM2808311_CHAST1_S1_all_R1_001_count.txt"
data_chrom_1 <- read.delim(chrom_datafile_1, header = FALSE)
colnames(data_chrom_1) <- c("Gene_ID", "1_1")
rownames(data_chrom_1) <- data_chrom_1$Gene_ID
sub_1 <- data_chrom_1[,2, drop=FALSE]

chrom_datafile_2 <- "GSE104816_GSM2808312_CHAST2_S2_all_R1_001_count.txt"
data_chrom_2 <- read.delim(chrom_datafile_2, header = FALSE)
colnames(data_chrom_2) <- c("Gene_ID", "1_2")
rownames(data_chrom_2) <- data_chrom_2$Gene_ID
sub_2 <- data_chrom_2[,2, drop=FALSE]

chrom_datafile_3 <- "GSE104816_GSM2808313_CHAST3_S3_all_R1_001_count.txt"
data_chrom_3 <- read.delim(chrom_datafile_3, header = FALSE)
colnames(data_chrom_3) <- c("Gene_ID", "1_3")
rownames(data_chrom_3) <- data_chrom_3$Gene_ID
sub_3 <- data_chrom_3[,2, drop=FALSE]

chrom_datafile_4 <- "GSE67058_edgePro_counts_subset.txt"
data_chrom_4 <- read.delim(chrom_datafile_4, header = TRUE)
rownames(data_chrom_4) <- data_chrom_4$gene
sub_4 <- data_chrom_4[,-1]

chrom_datafile_5 <- "GSE93894_SLX085-02-v4.1002vs375_subset_gene_expression.txt"
data_chrom_5 <- read.delim(chrom_datafile_5, header =TRUE)
colnames(data_chrom_5) <- c("Gene_ID", "Sample3_median")
rownames(data_chrom_5) <- data_chrom_5$Gene_ID
sub_5 <- data_chrom_5[,2, drop=FALSE]

#this data set needs to have gene names changed to gene locus
chrom_datafile_6 <- "GSE80786_counts_matrix.txt"
data_chrom_6 <- read.delim(chrom_datafile_6, header = TRUE)
rownames(data_chrom_6) <- data_chrom_6$GeneName
sub_6 <- data_chrom_6[,c(-1,-5)]
sub_6[is.na(sub_6)] <- 0

#below is the file that links the gene name to locus tag
IDs <- "Bacillus_subtilis_168_NC_000964_gene_ID_and_locus_tag.txt"
IDs_df <- read.delim(IDs, header = FALSE)
colnames(IDs_df) <- c("locus_tag","Gene_ID")
#remove dups
IDs_df <- unique(IDs_df)
row.names(IDs_df) <- IDs_df$locus_tag
gene_id <- as.character(IDs_df$Gene_ID)
gene_id <- make.unique(gene_id)
IDs_df$Gene_ID <- gene_id
row.names(IDs_df) <- IDs_df$Gene_ID

#merging the ID data with the samples
ID_exp_merged <- Reduce(merge, lapply(list(IDs_df, sub_6), function(x) data.frame(x, rn = row.names(x))))
row.names(ID_exp_merged) <- ID_exp_merged$locus_tag
ID_exp_merged <- ID_exp_merged[,-c(1,2,3)]

#merging each dataset by row name into one large dataset
exp_dat_all_reps <- Reduce(merge, lapply(list(sub_1, sub_2, sub_3, sub_4,sub_5,ID_exp_merged), function(x) data.frame(x, rn = row.names(x))))
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

#combining replicates by finding the median or each group of replicates
exp_df$Sample1_median <- apply(exp_df[,c(1,2,3)], 1, median.default)
exp_df$Sample2_median <- apply(exp_df[,c(4,5,6,7,8,9)], 1, median.default)
exp_df$Sample4_median <- apply(exp_df[,c(11,12,13)], 1, median.default)

#combining each dataset by finding median of the medians of each dataset
exp_df$Median <- apply(exp_df[,c(10,14,15,16)], 1, median.default)
head(exp_df)
mean_df <- exp_df[,c(10,14,15,16)]


#making dataframe with max - mean
mean_df <- abs(mean_df - rowMeans(mean_df)) 
head(mean_df)
geneID <- row.names(normalized_exp_df)
mean_df <- cbind(mean_df,geneID)
mean_df_melt <- melt(mean_df, id.vars = "geneID", value.name = "sample_num")
head(mean_df_melt)
colnames(mean_df_melt) <- c("gene_name", "sample_num", "exp_var") 
head(mean_df_melt)
class(mean_df_melt$gene_name)
class(mean_df_melt$sample_num)
class(mean_df_melt$exp_var)

mytitle <- expression(paste(italic("B.subtilis"), " Chromosome"))
med_exp_graph <- ggplot(mean_df_melt, aes(x = gene_name, y = log10(exp_var), color = sample_num)) + 
  ggtitle(mytitle) +
  xlab("Genes") +
  ylab("|CPM - Mean CPM| (log10)") +
  labs(color = "GEO Dataset")+
  scale_color_hue(labels = c("GSE93894", "GSE104816", "GSE67058", "GSE80786")) +
#  scale_color_manual(labels = c("G", "S", "E", "6")) +
  #  geom_point(position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4),
             #             alpha = 1/10,
             shape = 1) +
  theme(axis.text.x = element_blank())
pdf("bass_median_mean_gene_exp_check.pdf")
med_exp_graph
dev.off()

median_exp_df <- exp_df[,-c(1:16)]

geneID <- row.names(normalized_exp_df)
exp_dat_ID <- cbind(geneID, median_exp_df)

  

write.csv(exp_dat_ID, file = "Bass_chrom_gene_expression_combined_normalized.txt", row.names = FALSE)

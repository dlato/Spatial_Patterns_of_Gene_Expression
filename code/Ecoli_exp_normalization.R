###########################
#ECOLI CHROM 
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




chrom_datafile_1 <- "GSE60522_RNAseq_ProcessedData_subset.txt"
data_chrom_1 <- read.delim(chrom_datafile_1, header = TRUE)
rownames(data_chrom_1) <- data_chrom_1$gene
sub_1 <- data_chrom_1[,-1]
head(sub_1)

chrom_datafile_10 <- "GSE114917_count_table.txt"
data_chrom_10 <- read.delim(chrom_datafile_10, header = TRUE)
rownames(data_chrom_10) <- data_chrom_10$GeneName
sub_11 <- data_chrom_10[,-1]

#below 5 samples need to have their IDs converted to gene name
chrom_datafile_5 <- "GSE73673_GSM1900506_WT_1.htcount.txt"
data_chrom_5 <- read.delim(chrom_datafile_5, header =FALSE)
colnames(data_chrom_5) <- c("Gene_ID", "3_1")
rownames(data_chrom_5) <- data_chrom_5$Gene_ID
sub_5 <- data_chrom_5[,2, drop=FALSE]

chrom_datafile_6 <- "GSE73673_GSM1900507_WT_2.htcount.txt"
data_chrom_6 <- read.delim(chrom_datafile_6, header =FALSE)
colnames(data_chrom_6) <- c("Gene_ID", "3_2")
rownames(data_chrom_6) <- data_chrom_6$Gene_ID
sub_6 <- data_chrom_6[,2, drop=FALSE]

chrom_datafile_7 <- "GSE73673_GSM1900508_WT_3.htcount.txt"
data_chrom_7 <- read.delim(chrom_datafile_7, header =FALSE)
colnames(data_chrom_7) <- c("Gene_ID", "3_3")
rownames(data_chrom_7) <- data_chrom_7$Gene_ID
sub_7 <- data_chrom_7[,2, drop=FALSE]

chrom_datafile_2 <- "GSE40313_WT1_GSM991216_AW1.fastq.sam.hits.expression.txt"
data_chrom_2 <- read.delim(chrom_datafile_2, header = FALSE)
colnames(data_chrom_2) <- c("Gene_ID", "2_1")
rownames(data_chrom_2) <- data_chrom_2$Gene_ID
sub_2 <- data_chrom_2[,2, drop=FALSE]

chrom_datafile_3 <- "GSE40313_WT2_GSM991217_AW2.fastq.sam.hits.expression.txt"
data_chrom_3 <- read.delim(chrom_datafile_3, header = FALSE)
colnames(data_chrom_3) <- c("Gene_ID", "2_2")
rownames(data_chrom_3) <- data_chrom_3$Gene_ID
sub_3 <- data_chrom_3[,2, drop=FALSE]

chrom_datafile_8 <- "GSE54199_EC_Cont_DNA.txt"
data_chrom_8 <- read.delim(chrom_datafile_8, sep = "\t", fileEncoding="utf-16",header =TRUE)
data_chrom_8 <- data_chrom_8[,c(-2,-3,-4,-5,-8)]
colnames(data_chrom_8) <- c("Gene_ID","4_1", "4_2")
rownames(data_chrom_8) <- data_chrom_8$Gene_ID
sub_8 <- data_chrom_8[,-1]

chrom_datafile_9 <- "GSE54199_EC_Cont_RNA.txt"
data_chrom_9 <- read.delim(chrom_datafile_9, header = TRUE)
#data_chrom_9 <- read.delim(chrom_datafile_9, sep = "\t", fileEncoding="utf-16",header =FALSE)
data_chrom_9 <- data_chrom_9[,c(-2,-3,-4,-5,-8)]
colnames(data_chrom_9) <- c("Gene_ID","5_1", "5_2")
rownames(data_chrom_9) <- data_chrom_9$Gene_ID
sub_12 <- data_chrom_9[,-1]

#merging this sample bc the geneIDs have to be edited
exp_dat_all_reps <- Reduce(merge, lapply(list(sub_5, sub_6, sub_7, sub_2, sub_3, sub_8, sub_12), function(x) data.frame(x, rn = row.names(x))))
row.names(exp_dat_all_reps) <- exp_dat_all_reps$rn
sample3_all_reps <- exp_dat_all_reps[,-1]

#below is the file that links the gene name to locus tag
IDs <- "Ecoli_K12_MG1655_chrom_U00096_geneID_locus_tag.txt"
IDs_df <- read.delim(IDs, header = FALSE)
colnames(IDs_df) <- c("locus_tag","Gene_ID")
#remove dups
IDs_df <- unique(IDs_df)
row.names(IDs_df) <- IDs_df$locus_tag
gene_id <- as.character(IDs_df$Gene_ID)
gene_id <- make.unique(gene_id)
IDs_df$Gene_ID <- gene_id

#merging the ID data with the samples
ID_exp_merged <- Reduce(merge, lapply(list(IDs_df, sample3_all_reps), function(x) data.frame(x, rn = row.names(x))))
row.names(ID_exp_merged) <- ID_exp_merged$Gene_ID
ID_exp_merged <- ID_exp_merged[,-c(1,2,3)]

#below dataset has link btwn gene name and gene ID
chrom_datafile_9 <- "GSE85914_expression_subset.txt"
data_chrom_9 <- read.delim(chrom_datafile_9, header =TRUE)
colnames(data_chrom_9) <- c("Gene_ID", "Syn", "Exp")
Gene_ID <- as.character(data_chrom_9$Gene_ID)
Gene_ID <- make.unique(Gene_ID)
data_chrom_9$Gene_ID <- Gene_ID
rownames(data_chrom_9) <- data_chrom_9$Gene_ID
sub_9 <- data_chrom_9[,3, drop=FALSE]

chrom_datafile_10 <- "GSE98890_gene_read_counts_subset.txt"
data_chrom_10 <- read.delim(chrom_datafile_10, header =TRUE)
colnames(data_chrom_10) <- c("Gene_ID", "10_1", "10_2")
data_chrom_10 <- unique(data_chrom_10)
Gene_ID <- as.character(data_chrom_10$Gene_ID)
Gene_ID <- make.unique(Gene_ID)
data_chrom_10$Gene_ID <- Gene_ID
rownames(data_chrom_10) <- data_chrom_10$Gene_ID
sub_10 <- data_chrom_10[,-1]


#merging each dataset by row name into one large dataset
exp_dat_all_reps <- Reduce(merge, lapply(list(sub_1, sub_11, ID_exp_merged, sub_9, sub_10), function(x) data.frame(x, rn = row.names(x))))
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


######################################
# checking if all samples have similar expression levels for each gene
######################################
#first restructuring the dataframe
#getting the gene names into the df
geneID <- row.names(normalized_exp_df)
per_gene_exp_test_df <- exp_df
rownames(per_gene_exp_test_df) <- geneID
head(per_gene_exp_test_df)
#picking 100 random genes to see if they cluster
gene_nums <- sample(1:length(rownames(per_gene_exp_test_df)),100,replace=F)
sub_sample_exp <- per_gene_exp_test_df[gene_nums,]
#sub_sample_exp <- t(sub_sample_exp)
head(sub_sample_exp)

#making dataframe with max - mean
test_df <- abs(sub_sample_exp - rowMeans(sub_sample_exp)) 
sample_num <- colnames(test_df)
head(test_df)
#making dataframe ggplot will like
gene_name <- rownames(test_df)
df <- cbind(test_df, gene_name)
class(df)
df_melt <- melt(df, id.vars = "gene_name", value.name = "sample_num")
colnames(df_melt) <- c("gene_name", "sample_num", "exp_var") 
head(df_melt)
class(df_melt$gene_name)
class(df_melt$sample_num)
class(df_melt$exp_var)

ggplot(df_melt, aes(x = gene_name, y = exp_var, color = sample_num)) + 
  geom_point(position = position_dodge(width = 0.4)) +
  theme(legend.position = "none")



#combining replicates by finding the median for each group of replicates
exp_df$Sample1_median <- apply(exp_df[,c(1,2,3)], 1, median.default)
exp_df$Sample2_median <- apply(exp_df[,c(4,5,6,7)], 1, median.default)
exp_df$Sample3_median <- apply(exp_df[,c(8,9,10)], 1, median.default)
exp_df$Sample4_median <- apply(exp_df[,c(11,12)], 1, median.default)
exp_df$Sample5_median <- apply(exp_df[,c(13,14)], 1, median.default)
exp_df$Sample6_median <- apply(exp_df[,c(15,16)], 1, median.default)
exp_df$Sample8_median <- apply(exp_df[,c(18,19)], 1, median.default)

#combining each dataset by finding median of the medians of each dataset
exp_df$Median <- apply(exp_df[,c(17,20,21,22,23,24,25,26)], 1, median.default)
head(exp_df)
mean_df <- exp_df[,c(20,21,22,23,24,25,26,17)]
head(mean_df)

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


mytitle <- expression(paste(italic("E.coli"), " Chromosome"))
med_exp_graph <- ggplot(mean_df_melt, aes(x = gene_name, y = log10(exp_var), color = sample_num)) + 
  ggtitle(mytitle) +
  xlab("Genes") +
  ylab("|CPM - Mean CPM| (log10)") +
  labs(color = "GEO Dataset")+
  scale_color_hue(labels = c("GSE60522", "GSE114917", "GSE73673", "GSE40313", "GSE54199", "GSE54199", "GSE98890", "GSE85914")) +
  #  geom_point(position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4),
             #             alpha = 1/10,
             shape = 1) +
  theme(axis.text.x = element_blank())
pdf("ecoli_median_mean_gene_exp_check.pdf")
med_exp_graph
dev.off()


median_exp_df <- exp_df[,-c(1:26)]

geneID <- row.names(normalized_exp_df)



exp_dat_ID <- cbind(geneID, median_exp_df)

  

write.csv(exp_dat_ID, file = "Ecoli_chrom_gene_expression_combined_normalized.txt", row.names = FALSE)


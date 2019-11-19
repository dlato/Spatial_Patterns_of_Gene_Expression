#will calculate the weighted mean for dN and dS
#to run: Rscript average_dN_dS NAME_OF_DATAFRAME_FILE BACTERIA_NAME REPLICON_NAME GENOME_POSITION_DATAFRAME_FILE LENGTH_OF_GENOME ORIC_POS TERMINUS_POS
#to run: Rscript average_dN_dS NAME_OF_DATAFRAME_FILE BACTERIA_NAME REPLICON_NAME LENGTH_OF_GENOME ORIC_POS TERMINUS_POS
# order of args:
# 1) expression file (csv)
# 2) bacteria name
# 3) replicon name 
# 4) max genome position
# 5) origin of replication location
# 6) terminus of replication location
# 7) length of chunks
# 8) out file bac name
#install.packages('ggplot2', dep = TRUE)
#install.packages('proto', dep = TRUE)
library("RColorBrewer")
loadNamespace("ggplot2")
library("gplots")
library("ggplot2")
library(ggplot2)
library(plyr)
library(dplyr)
library(proto)
library("proto")
is.proto.proto <- is.proto
library(reshape)
library(testthat)
library(gridExtra)
library(grid)
library(lattice)

#getting command line arguments
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
bac_name <- as.character(args[2])
replicon <- as.character(args[3])
out_file_name <- as.character(args[8])


options("scipen"=100, "digits"=5)
datafile <- file_name
merged_data <- read.table(datafile, sep = ",", header = TRUE)
head(merged_data)

mytitle <- substitute(italic(bac_name)~replicon, list(bac_name=bac_name, replicon=replicon))
options("scipen"=100, "digits"=10)
################################################################################
#ORIGIN SCALING AND BIDIRECTIONALITY
################################################################################
#first scaling things to the origin (if necessary)
max_pos <- as.numeric(args[4])
print("max_pos")
max_pos
oriC_pos <- as.numeric(args[5])
print("oriC")
oriC_pos
terminus <- as.numeric(args[6])
print("ter")
terminus
new_pos <- merged_data$Midpoint
tmp_pos <- merged_data$Midpoint
print("MIN POS")
min(merged_data$Midpoint)
head(merged_data)
 
if (bac_name == "E.coli" | replicon == "pSymA") {
  to_shift_ter <- max_pos - oriC_pos
  shifted_ter <-terminus + to_shift_ter
  terminus <- shifted_ter
}
print("shifted ter")
terminus

if (replicon == "pSymB") {
  shifted_ter <- terminus - oriC_pos
  terminus <- shifted_ter
}

if (bac_name == "E.coli" | replicon == "pSymA" | replicon == "pSymB") {
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    } else {
      tmp_end <- max_pos - oriC_pos
      new_pos[i] <- tmp_pos[i] + tmp_end
    }
  }
  tmp_pos <- new_pos
}
 
 
#now accounting for the bidirectionality. if things are between the start pos and
#the terminus then they will stay as the same position. If not, then they will be
#changed to a new position starting at 1 and going to the terminus
new_pos2 <- tmp_pos
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name == "S.meliloti") {
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] > terminus) {
      new_pos2[i] <- max_pos - tmp_pos[i]
    } else {
    }
  }
  tmp_pos <- new_pos2
   
  print("max tmp_pos")
  max(tmp_pos)
}
 

if (bac_name == "Streptomyces") {
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    }#if btwn origin and end of genome 
    if (tmp_pos[i] <= oriC_pos) {
      new_pos[i] <- -1 * (oriC_pos - tmp_pos[i])
    }#if btwn origin and beginning of genome
    if (tmp_pos[i] == oriC_pos) {
      new_pos[i] <- 0
    }#if equal to origin
  }
  tmp_pos <- new_pos
}


merged_data$Midpoint <- tmp_pos
#merged_data <- as.data.frame(cbind(merged_data$block, merged_data$gene, merged_data$sec, tmp_pos, merged_data$dS, merged_data$dN, merged_data$omega, merged_data$sec_len))
colnames(merged_data) <- c("X", "Id", "Start","End", "tmp_pos","Exp")
head(merged_data)
max(merged_data$tmp_pos)
min(merged_data$tmp_pos)
merged_data
 
#################################################################################
#################################################################################
## GENE EXPRESSION
#################################################################################
#################################################################################
mean(merged_data$Exp)

#new min and max for the scaled positions
min(merged_data$tmp_pos)
nmin_pos <- round_any(min(merged_data$tmp_pos), 10000, f=floor)
nmin_pos
nmax_pos <- max(merged_data$tmp_pos)
nmax_pos

#expty vec to hold the rows that split the dat into X kb chunks
exp_rows_to_split_dat <- vector()
#create vector of chunks
chunklen <- as.numeric(args[7])
chunklen
chunk_len_of_genome <- round_any(nmax_pos, 10000, f=ceiling) 
chunk_len_of_genome
chunks <- seq(nmin_pos, chunk_len_of_genome, 10000)
chunks
print("merged_data")
head(merged_data)
#order the data by positions, only if its NOT strep
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name == "S.meliloti") {
	ord_pos <- order(merged_data$tmp_pos)
	merged_data$tmp_pos <- merged_data$tmp_pos[ord_pos]
	merged_data$X <- merged_data$X[ord_pos]
	merged_data$Id <- merged_data$Id[ord_pos]
	merged_data$Start <- merged_data$Start[ord_pos]
	merged_data$End <- merged_data$End[ord_pos]
	merged_data$Exp <- merged_data$Exp[ord_pos]
}
print("BIDIRECTIONAL DATA")
merged_data
write.table(merged_data, 'bidirectional_data.csv', sep = "\t")

for (i in chunks) {
  exp_rows <-
which(abs(merged_data$tmp_pos-i)==min(abs(merged_data$tmp_pos-i)))
# finding the closest number to each 10kb without going over it
  exp_max_row <- max(exp_rows)
#  actual_pos <- data_chrom_ordered$midpoint[max_row]
#  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1
exp_rows_to_split_dat

#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
data_just_exp <- merged_data$Exp
#list_dat_sets <- split(data_chrom_ordered, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))
list_dat_sets_just_exp <- split(data_just_exp, findInterval(1:nrow(merged_data), exp_rows_to_split_dat))
list_dat_sets_just_exp

#############################################################
# TOTAL EXPRESSION
#############################################################
print("#############################################################")
print("# TOTAL EXPRESSION")
print("#############################################################")
##########
#add up all expression values in each 10kb section
##########
list_tot_add_exp_10kb <- lapply(list_dat_sets_just_exp, sum)
print("list_total_exp_add_10kb")
list_tot_add_exp_10kb
tot_gene_exp_10kb <- as.data.frame(matrix(unlist(list_tot_add_exp_10kb), byrow = F))
head(tot_gene_exp_10kb)
tail(tot_gene_exp_10kb)

#checking if the lengths are equal for the later bind
if (length(unlist(tot_gene_exp_10kb)) == length(chunks)) {
    chunks_pos_NOzero <- chunks   
} else { 
    chunks_pos_NOzero <- chunks[which(chunks != 0)]
}
chunks_pos_NOzero
print("LEN")
length(chunks_pos_NOzero)
length(unlist(tot_gene_exp_10kb))
tot_gene_exp_10kb <- cbind(tot_gene_exp_10kb, chunks_pos_NOzero)
head(tot_gene_exp_10kb)
tail(tot_gene_exp_10kb)
tot_gene_exp_10kb <- tot_gene_exp_10kb[-nrow(tot_gene_exp_10kb),]
write.csv(tot_gene_exp_10kb)
######################################################
# OUTLIERS (colouring them differntly in the graph)
######################################################
#calculating outliers (function)
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  #response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
 # if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
#  } else{
 #   cat("Nothing changed", "n")
  #  return(invisible(var_name))
  #}
}

print("##################################")
print("## CALCULATING OUTLIERS FOR EXP ADD 10KB ##")
print("##################################")
#setting new df to remove outliers
exp_add_10kb_outliers <- tot_gene_exp_10kb
#to run function:
#outlierKD(data_frame, column_with_data)
#finding outlier bars in the weighted subs data
outlierKD(exp_add_10kb_outliers, V1)
print("checking outliers df")
head(exp_add_10kb_outliers)
tail(exp_add_10kb_outliers)

exp_add_10kb_outlier_rows <- which(is.na(exp_add_10kb_outliers$V1))
print("rows where the bars are outliers")
exp_add_10kb_outlier_rows

exp_add_10kb_out_pos <- vector()
for (i in exp_add_10kb_outlier_rows) {
  exp_add_10kb_out_pos <- c(exp_add_10kb_out_pos, exp_add_10kb_outliers$chunks_pos_NOzero[i])
}
print("outlier genomic positions (bar position)")
exp_add_10kb_out_pos




######################################################
#################
#do lm on absoloute value to account for both replichores
print("### LM ON ABSOLOUTE VALUE OF TOTAL (added) GENE EXP ###")
lm_tot_exp <- lm(abs(V1) ~ chunks_pos_NOzero, data=exp_add_10kb_outliers)
summary(lm_tot_exp)

#################
##doing lm on each of the replichores (-tive and +tive numbers)
#### Left replichore ###
#print("### Left replichore ###")
#tot_exp_left_rep <- tot_gene_exp_10kb[which(tot_gene_exp_10kb$chunks_pos_NOzero <= 0 ),]
#head(tot_exp_left_rep)
#tail(tot_exp_left_rep)
#print("REMEMBER FOR THIS ONE THE ORI IS ON THE RIGHT!!! (SO EXPECT LM RESULTS TO BE OPPOSITE SIGN")
#lm_tot_exp_left_rep <- lm(V1 ~ chunks_pos_NOzero, data=tot_exp_left_rep)
#summary(lm_tot_exp_left_rep)
#### Right replichore ###
#print("### Right replichore ###")
#tot_exp_right_rep <- tot_gene_exp_10kb[which(tot_gene_exp_10kb$chunks_pos_NOzero > 0 ),]
#head(tot_exp_right_rep)
#tail(tot_exp_right_rep)
#lm_tot_exp_right_rep <- lm(V1 ~ chunks_pos_NOzero, data=tot_exp_right_rep)
#summary(lm_tot_exp_right_rep)
#################
print("#############################################################")
print("# MEDIAN")
print("#############################################################")
#########
# 10 kb linear reg for gene expression of median over 10kb
#########
#get median expression in each of the sections/ 10kb chunks
list_total_exp_level <- lapply(list_dat_sets_just_exp, median)
print("list_total_exp_level")
list_total_exp_level
#as df
median_gene_exp_10kb <- as.data.frame(matrix(unlist(list_total_exp_level), byrow = F))
median_gene_exp_10kb <- cbind(median_gene_exp_10kb, chunks_pos_NOzero)
write.csv(median_gene_exp_10kb)
write.table(median_gene_exp_10kb, 'median_10kb_exp_data.csv', sep = "\t")

print("##################################")
print("## CALCULATING OUTLIERS FOR MEDIAN EXP over 10KB ##")
print("##################################")
#setting new df to remove outliers
median_outliers <- median_gene_exp_10kb
#to run function:
#outlierKD(data_frame, column_with_data)
#finding outlier bars in the weighted subs data
outlierKD(median_outliers, V1)
print("checking outliers df")
head(median_outliers)
tail(median_outliers)

median_outlier_rows <- which(is.na(median_outliers$V1))
print("rows where the bars are outliers")
median_outlier_rows

median_out_pos <- vector()
for (i in median_outlier_rows) {
  median_out_pos <- c(median_out_pos, median_outliers$chunks_pos_NOzero[i])
}
print("outlier genomic positions (bar position)")
median_out_pos


##colours for histogram (outliers and reg bars)
#################
#do lm on absoloute value to account for both replichores
print("### LM ON ABSOLOUTE VALUE OF MEDIANs###")
lm_med_exp <- lm(abs(V1) ~ chunks_pos_NOzero, data=median_outliers)
summary(lm_med_exp)
#################
print("####################")
print("high median exp values")
print("####################")
median_gene_exp_10kb[which(median_gene_exp_10kb$V1 >= 400),]
#################
##doing lm on each of the replichores (-tive and +tive numbers)
#### Left replichore ###
#print("### Left replichore ###")
#median_exp_left_rep <- median_gene_exp_10kb[which(median_gene_exp_10kb$chunks_pos_NOzero <= 0 ),]
#print("REMEMBER FOR THIS ONE THE ORI IS ON THE RIGHT!!! (SO EXPECT LM RESULTS TO BE OPPOSITE SIGN")
#lm_median_exp_left_rep <- lm(V1 ~ chunks_pos_NOzero, data=median_exp_left_rep)
#summary(lm_median_exp_left_rep)
#### Right replichore ###
#print("### Right replichore ###")
#median_exp_right_rep <- median_gene_exp_10kb[which(median_gene_exp_10kb$chunks_pos_NOzero > 0 ),]
#lm_median_exp_right_rep <- lm(V1 ~ chunks_pos_NOzero, data=median_exp_right_rep)
#summary(lm_median_exp_right_rep)
################
#############################################################
#############################################################
##############LOGISTIC REGRESSION FOR EXP DAT###########
#############################################################
print("#############################################################")
print("##############LOGISTIC REGRESSION FOR EXP DAT ALL POINTS (no average over 10kb)###########")
print("#############################################################")
#################
#do lm on absoloute value to account for both replichores
print("### LM ON ABSOLOUTE VALUE OF ALL EXP DAT PTS###")
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=merged_data)
summary(lm_all_exp)
#################
## run a reg linear regression on the exp data
#### Left replichore ###
#print("### Left replichore ###")
#all_exp_left_rep <- merged_data[which(merged_data$tmp_pos <= 0 ),]
#print("REMEMBER FOR THIS ONE THE ORI IS ON THE RIGHT!!! (SO EXPECT LM RESULTS TO BE OPPOSITE SIGN")
#lm_all_exp_left_rep <- lm(Exp ~ tmp_pos, data=all_exp_left_rep)
#summary(lm_all_exp_left_rep)
#### Right replichore ###
#print("### Right replichore ###")
#all_exp_right_rep <- merged_data[which(merged_data$tmp_pos > 0 ),]
#lm_all_exp_right_rep <- lm(Exp ~ tmp_pos, data=all_exp_right_rep)
#summary(lm_all_exp_right_rep)
print("#############################################################")
print("##############AVERAGE GENE EXPRESSION ACROSS GENOME ###########")
print("#############################################################")
avg_gene_exp <- mean(abs(merged_data$Exp))
avg_gene_exp
merged_data[which(merged_data$Exp >= 100),]
length(which(merged_data$Exp >= 100))
#############################################################
#############################################################
# TOTAL NUMBER OF GENES
##########################################################################
print("#############################################################")
print("# TOTAL NUMBER OF GENES")
print("##########################################################################")
#count the total number of genes in each of the sections/ 10kb chunks
list_tot_gene_num <- lapply(list_dat_sets_just_exp, length)
#put list of num of genes into a new dataframe the ggplot can recognize
gene_num <- as.data.frame(matrix(unlist(list_tot_gene_num), byrow = F))
gene_num_data <- cbind(gene_num, chunks_pos_NOzero)
write.csv(gene_num_data)
print("MEAN GENE NUM PER 10KB")
mean(gene_num_data$V1)
gene_num_data$V1[which(gene_num_data$V1 ==0)]
########### LINEAR REG ON NUM OF GENES ##########
#################
#do lm on absoloute value to account for both replichores
print("### Lin reg on num of genes###")
lm_num_genes <- lm(abs(V1) ~ chunks_pos_NOzero, data=gene_num_data)
summary(lm_num_genes)
#################
#print("### Left replichore ###")
#gene_num_left_rep <- gene_num_data[which(gene_num_data$chunks_pos_NOzero <= 0 ),]
#print("REMEMBER FOR THIS ONE THE ORI IS ON THE RIGHT!!! (SO EXPECT LM RESULTS TO BE OPPOSITE SIGN")
#lm_gene_num_left_rep <- lm(V1 ~ chunks_pos_NOzero, data=gene_num_left_rep)
#summary(lm_gene_num_left_rep)
#print("### Right replichore ###")
#gene_num_right_rep <- gene_num_data[which(gene_num_data$chunks_pos_NOzero > 0 ),]
#lm_gene_num_right_rep <- lm(V1 ~ chunks_pos_NOzero, data=gene_num_right_rep)
#summary(lm_gene_num_right_rep)


#############################################################
# NORMALIZING TOTAL EXP BY NUMBER OF GENES
##########################################################################
print("#############################################################")
print("# NORMALIZING TOTAL EXP BY NUMBER OF GENES")
print("##########################################################################")
#putting the exp dat and gene num dat into one df
tot_exp_gene_num_dat <- merge(tot_gene_exp_10kb, gene_num_data, by="chunks_pos_NOzero")
head(tot_exp_gene_num_dat)
tot_exp_gene_num_dat
######## normalizing by num of genes
tot_exp_gene_num_dat$norm <- tot_exp_gene_num_dat$V1.x / tot_exp_gene_num_dat$V1.y
head(tot_exp_gene_num_dat)

print("##################################")
print("## CALCULATING OUTLIERS FOR NORMALIZED EXP 10KB ##")
print("##################################")
#setting new df to remove outliers
normalized_outliers <- tot_exp_gene_num_dat
#to run function:
#outlierKD(data_frame, column_with_data)
#finding outlier bars in the weighted subs data
outlierKD(normalized_outliers, norm)
print("checking outliers df")
head(normalized_outliers)
tail(normalized_outliers)

normalized_outlier_rows <- which(is.na(normalized_outliers$norm))
print("rows where the bars are outliers")
normalized_outlier_rows

normalized_out_pos <- vector()
for (i in normalized_outlier_rows) {
  normalized_out_pos <- c(normalized_out_pos, normalized_outliers$chunks_pos_NOzero[i])
}
print("outlier genomic positions (bar position)")
normalized_out_pos


#################
#do lm on absoloute value to account for both replichores
print("### Lin reg on normalized expression###")
lm_norm_exp <- lm(abs(norm) ~ chunks_pos_NOzero, data=normalized_outliers)
summary(lm_norm_exp)
#################
#print("### Left replichore ###")
#normalized_left_rep <- tot_exp_gene_num_dat[which(tot_exp_gene_num_dat$chunks_pos_NOzero <= 0 ),]
#print("REMEMBER FOR THIS ONE THE ORI IS ON THE RIGHT!!! (SO EXPECT LM RESULTS TO BE OPPOSITE SIGN")
#lm_normalized_left_rep <- lm(norm ~ chunks_pos_NOzero, data=normalized_left_rep)
#summary(lm_normalized_left_rep)
#print("### Right replichore ###")
#normalized_right_rep <- tot_exp_gene_num_dat[which(tot_exp_gene_num_dat$chunks_pos_NOzero > 0 ),]
#lm_normalized_right_rep <- lm(norm ~ chunks_pos_NOzero, data=normalized_right_rep)
#summary(lm_normalized_right_rep)


######################################################
# PLOT
######################################################
#colours for histogram (outliers and reg bars)
#out_pos_colours <- rep("#BD93BD", length(rm_outliers$weighted_subs))
#out_pos_colours <- rep("black", length(median_outliers$V1))
# COLOURS 
#out_pos_colours <- rep("#709F83", length(median_outliers$V1))
#B&W PRINTING
out_pos_colours <- rep("#868686", length(median_outliers$V1))
for (i in median_outlier_rows) {
    #colour for an outlier
#  out_pos_colours[i] <- "#F2EDEB"
#  out_pos_colours[i] <- "#CAD2C5"
  out_pos_colours[i] <- "#D1D1D1"
#  out_pos_colours[i] <- "#E8D7F1"
#  out_pos_colours[i] <- "#C1CAD6"
}

out_pos_colours

#rm_outliers$outlier <- out_pos
#head(rm_outliers)
#tail(rm_outliers)
###########################
options(scipen=3)
print("gene num hist")
gene_num_data$V1 <- as.numeric(as.character(gene_num_data$V1))
class(gene_num_data$chunks_pos_NOzero)
gene_num_hist <- (ggplot(gene_num_data, aes(x = chunks_pos_NOzero, y = V1)) 
#  + geom_histogram(stat = "identity", fill= "black",boundary=0) 
#  + geom_histogram(fill= "black", binwidth = 10) 
#  + geom_bar(stat = "identity", fill= "black") 
#  + geom_bar(stat = "identity", width = 1) 
# COLOURS 
#  + geom_bar(stat = "identity", color = "#709F83", fill = "#709F83") 
#B&W PRINTING
  + geom_bar(stat = "identity", color = "#868686", fill = "#868686") 
#  geom_histogram(stat = "identity", fill= "#FE5F55") 
  + labs(x = "Distance from the Origin of Replication (bp)", y = "Number of Genes") 
  + ggtitle(mytitle) 
  + theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =
#                                     0.5,color = "black"),
        axis.text.y = element_text(),
        axis.title.y = element_text(),
        axis.ticks.y = element_line(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), legend.position = c(0.9,0.9),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_line(colour = "grey98", size = 0.5),
        panel.spacing =      unit(0.25, "lines")) 
#   + geom_vline(xintercept = 0, colour = "red")
   + geom_vline(xintercept = 0, colour = "black")
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =
  #  0.5)) +
  #scale_x_continuous(expand = c(0, 0), breaks = my_labs_mill, labels = my_labs_mill)
)
###########################

###########################
median_gene_exp_10kb$chunks_pos_NOzero <- median_gene_exp_10kb$chunks_pos_NOzero / 1000000
head(median_gene_exp_10kb)
class(median_gene_exp_10kb$V1)
median_gene_exp_10kb$V1 <- as.numeric(as.character(median_gene_exp_10kb$V1))
#median_gene_exp_10kb$chunks_pos_NOzero <- as.factor(median_gene_exp_10kb$chunks_pos_NOzero)
#exp_bar_top <- (ggplot(median_gene_exp_10kb)  
exp_bar_top <- (ggplot(median_gene_exp_10kb, aes(x = chunks_pos_NOzero, y = V1))  
#exp_bar_top <- (ggplot(text_dat, aes(x = chunks_pos_NOzero, y = V1))  
   +  geom_bar(stat = "identity", fill= out_pos_colours, colour = out_pos_colours) 
#   +  geom_bar(stat = "identity", fill= out_pos_colours) 
#   +  geom_bar(stat = "identity", fill= out_pos_colours, width = 1) 
#   +  geom_bar(aes(x = chunks_pos_NOzero, y = V1, fill= out_pos_colours), width=1) 
#  + geom_bar(stat = "identity", fill= "#FE5F55") 
   + scale_colour_manual(values = out_pos_colours)
   + scale_fill_manual(values = out_pos_colours)
   + labs(x = "Distance from the Origin of Replication (million bp)", y = "Median
CPM") 
#  + ggtitle(mytitle) 
  #+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,color = "black"), 
  #+ theme(axis.text.x = element_text(hjust = 1, vjust = 0.5,color = "black"), 
  + theme(axis.text.x = element_text(color = "black"), 
      axis.text.y = element_text(),
      axis.title.x = element_text(),
      axis.title.y = element_text(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      plot.title = element_text(),
      legend.title = element_blank(), legend.position = c(0.9,0.9),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "grey90", size = 0.2),
      panel.grid.minor = element_line(colour = "grey98", size = 0.5),
      panel.spacing =      unit(0.25, "lines")
      )
   + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
   + geom_vline(xintercept = 0, colour = "black")
#   + geom_vline(xintercept = 0, colour = "red")
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =
  #  0.5)) +
  #scale_x_continuous(expand = c(0, 0), breaks = my_labs_mill, labels = my_labs_mill)
)
###########################

###########################
##below puts all the above graphs together
## 1 is gene_num_hist
## 2 is exp_bar_bottom
lay <- rbind(c(1),
             c(2))

pdf(paste(out_file_name,"_exp_tot_gene_num_bidirectionality_colour.pdf", sep=""))
grid.newpage()
grid.draw(rbind(ggplotGrob(gene_num_hist), ggplotGrob(exp_bar_top), size = "last"))
# next line adds border
grid.rect(width = 0.99, height = 0.99, gp = gpar(lwd = 2, col = "black", fill = NA))
dev.off()
###########################


pdf("TEST_norm_exp.pdf")
ggplot(tot_exp_gene_num_dat, aes(x=chunks_pos_NOzero, y=norm)) +
  geom_point()
dev.off()


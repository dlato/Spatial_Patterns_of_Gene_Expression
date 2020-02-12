setwd("C:/Users/synch/Documents/PhD/Gene_Expression/COG_30Jan20/COG/R_Code_Strep_9Jun17/")
#path="C:/Users/Daniella/Documents/Sinorhizobium2015/COG/R_Code_Strep_9Jun17/"

###########################
#STREP CHROM 
###########################
library(ggplot2)
require(ggplot2)
require("ggplot2")
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)



#chrom_datafile <- "C:/Users/synch/Documents/COG/R_Code_Strep_9Jun17/Strep_bingchenggensis_BCW1_COG_category_and_positions_multiple_COGs_dealt_with.txt"
chrom_datafile <- "Strep_bingchenggensis_BCW1_COG_category_and_positions_multiple_COGs_dealt_with.txt"
data_chrom <- read.table(chrom_datafile)
colnames(data_chrom) <- c("COG","start","end")
data_chrom$midpoint <- apply(data_chrom[,c(2,3)], 1, mean)
data_chrom <- data_chrom %>% mutate(midpoint = round(midpoint, 0))

#accounting for the bidirectionality in strep and that both
#ends of the chromosomes (bc it is linear) are considered "far"
#from the origin of rep. slightly different than the other bacteria
##max_pos <- max(chrom_gaps_data$tmp_pos)
tmp_pos <- data_chrom$midpoint
new_pos <- tmp_pos
oriC_pos <- 3419363
##terminus <- 1678398
#for(i in 1:length(tmp_pos)) {
#  if (tmp_pos[i] >= oriC_pos) {
#    new_pos[i] <- tmp_pos[i] - oriC_pos
#  }#if btwn origin and end of genome 
#  if (tmp_pos[i] <= oriC_pos) {
#    new_pos[i] <- oriC_pos - tmp_pos[i]
#  }#if btwn origin and beginning of genome
#  if (tmp_pos[i] == oriC_pos) {
#    new_pos[i] <- 0
#  }#if equal to origin
#}
#
#tmp_pos <- new_pos
#
#data_chrom$midpoint <- tmp_pos
#
#data_chrom_ordered <- data_chrom[order(data_chrom$midpoint),]


for(i in 1:length(tmp_pos)) {
  # right replichore
  if (tmp_pos[i] >= oriC_pos) {
    new_pos[i] <- tmp_pos[i] - oriC_pos
  }#if btwn origin and end of genome
  # left replichore
  if (tmp_pos[i] <= oriC_pos) {
    new_pos[i] <- -1 * (oriC_pos - tmp_pos[i])
    ## making sure the strand column accounts for bidirectional rep
    #if (merged2$strand[i] == 0) {
    #  new_strand[i] <- 1
    #} else {
    #  new_strand[i] <- 0
    #}
  }#if btwn origin and beginning of genome
  if (tmp_pos[i] == oriC_pos) {
    new_pos[i] <- 0
  }#if equal to origin
}
tmp_pos <- new_pos


data_chrom$midpoint <- tmp_pos
##merged2 <- as.data.frame(cbind(merged2$block, merged2$gene, merged2$sec, tmp_pos, merged2$dS, merged2$dN, merged2$omega, merged2$sec_len))
#colnames(merged2) <- c("Id","X",  "Start","End", "tmp_pos","Exp", "strand")
head(data_chrom)

#EDIT FROM HERE THIS LITTLE CHUNK
#new min and max for the scaled positions
min(data_chrom$midpoint)
nmin_pos <- round_any(min(data_chrom$midpoint), 500000, f=floor)
nmin_pos
nmax_pos <- max(data_chrom$midpoint)
nmax_pos

#expty vec to hold the rows that split the dat into X kb chunks
rows_to_split_dat <- vector()
#create vector of chunks
chunklen <- 500000
chunklen
chunk_len_of_genome <- round_any(nmax_pos, 500000, f=ceiling)
chunk_len_of_genome
chunks <- seq(nmin_pos, chunk_len_of_genome, 500000)
chunks
head(data_chrom)

for (i in chunks) {
  exp_rows <-
    which(abs(data_chrom$midpoint-i)==min(abs(data_chrom$midpoint-i)))
  # finding the closest number to each 10kb without going over it
  exp_max_row <- max(exp_rows)
  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  rows_to_split_dat <- c(rows_to_split_dat, exp_max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1
rows_to_split_dat




#
#rows_to_split_dat <- vector() #empty vec  to hold the rows that split the dat into 10kb chunks
#len_of_genome <- 8000000 #rounded up to closest 10,000 so that all of the points are included
#chunks <- seq(500000, len_of_genome, 500000)
#for (i in chunks) {
#  rows <- which(abs(data_chrom_ordered$midpoint-i)==min(abs(data_chrom_ordered$midpoint-i))) # finding the closest number to each 10kb without going over it
#  max_row <- max(rows)
##  actual_pos <- data_chrom_ordered$midpoint[max_row]
##  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
#  rows_to_split_dat <- c(rows_to_split_dat, max_row)
#}#for
##rows_to_split_dat <- rows_to_split_dat + 1

#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
##data_just_COG <- data_chrom_ordered$COG
data_just_COG <- data_chrom$COG
#list_dat_sets <- split(data_chrom_ordered, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))
##list_dat_sets_just_COG <- split(data_just_COG, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))
list_dat_sets_just_COG <- split(data_just_COG, findInterval(1:nrow(data_chrom), rows_to_split_dat))

#count the number of times each COG appears in each dataset
#frequency wise
list_freq_of_COGS <- lapply(list_dat_sets_just_COG, table)

#put list of freqs into a new dataframe the ggplot can recognize
COG_count <- as.data.frame(matrix(unlist(list_freq_of_COGS), byrow = F))
COG_cat <- levels(data_chrom$COG)
COG_cat_col <- rep(COG_cat, length(list_freq_of_COGS))
##below is for having the first chunk labled 0
##chunks_pos_zero <- c(0, chunks)
#below is for having the first chunk labled 100000
##chunks_pos_zero <- c(chunks, 8500000)
chunks_pos_zero <- chunks
chunk_pos <- rep(chunks_pos_zero, each = length(COG_cat))
freq_dat <- cbind(COG_count, COG_cat_col, chunk_pos)

#tick lables
my_labs_mill <- seq(nmin_pos,nmax_pos, by = 500000)


#make stacked bar graph
options(scipen=5)
p_stacked_bar <- ggplot(freq_dat, aes(x = chunk_pos, y = V1, fill = COG_cat_col)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() +
  labs(x = "Position in Genome (bp)", y = "% of COG Categories") +
  scale_fill_discrete(name="COG Category") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0.01)) +
  scale_x_continuous(expand = c(0, 0), breaks = my_labs_mill, labels = my_labs_mill)

##make bar graph of COGs (no freq)
p_bar_top <- ggplot(freq_dat, aes(x = chunk_pos, y = V1)) + 
  geom_bar(stat = "identity") + theme_bw() +
  ggtitle(expression(paste(italic("Streptomyces"), " Chromosome")))+
  labs(x = "Position in Genome (bp)", y = "Total Number of Genes") +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank() ) +
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0), breaks = my_labs_mill, labels = my_labs_mill)

##get and save the legend from the COG stacked
##bar graph
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(p_stacked_bar)

##remove the legend from the stacked bar graph
p_stacked_bar <- p_stacked_bar + theme(legend.position="none")

##below puts all the above graphs together
## 1 is p_bar_top
## 2 is the legend
## 3 is p_stacked_bar
lay <- rbind(c(1,2),
             c(3,2))

pdf("strep_chrom_COG_bargraph_and_hist_TEST.pdf")
grid.arrange(p_bar_top, legend, p_stacked_bar, layout_matrix = lay, widths=c(4,1), heights=c(1.5, 3.5))
#             top="Streptomyces Chromosome")
#grid.arrange(p_bar_top, legend, p_stacked_bar, ncol=2, nrow=2, widths=c(4.5,0.5), heights=c(1.5, 3.5), top="S.meliloti Chromosome")
dev.off()


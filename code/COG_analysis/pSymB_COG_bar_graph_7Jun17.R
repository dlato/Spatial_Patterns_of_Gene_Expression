setwd("C:/Users/synch/Documents/COG/R_Code_pSymB_7Jun17/")
path="C:/Users/Daniella/Documents/Sinorhizobium2015/COG/R_Code_pSymB_7Jun17/"

###########################
#PSYMB 
###########################
library(ggplot2)
require(ggplot2)
require("ggplot2")
library(reshape2)
library(dplyr)



chrom_datafile <- "C:/Users/synch/Documents/COG/R_Code_pSymB_7Jun17/Smel_1021_pSymB_COG_category_and_positions_multiple_COGs_dealt_with.txt"
data_chrom <- read.table(chrom_datafile)
colnames(data_chrom) <- c("COG","start","end")
data_chrom$midpoint <- apply(data_chrom[,c(2,3)], 1, mean)
data_chrom <- data_chrom %>% mutate(midpoint = round(midpoint, 0))

##account for bidirectionality of replication.
##using only the midpoint as the position. not scaling the 
#first scaling things to the origin (if necessary)
tmp_pos <- data_chrom$midpoint
max_pos <- max(data_chrom$midpoint)
new_pos <- tmp_pos
oriC_pos <- 55090
terminus <- 896756
for(i in 1:length(tmp_pos)) {
  if (tmp_pos[i] >= oriC_pos) {
    new_pos[i] <- tmp_pos[i] - oriC_pos
  } else {
    tmp_end <- max_pos - oriC_pos
    new_pos[i] <- tmp_pos[i] + tmp_end
  }
}

tmp_pos <- new_pos

#now accounting for the bidirectionality. if things are between the start pos and
#the terminus then they will stay as the same position. If not, then they will be
#changed to a new position starting at 1 and going to the terminus
new_pos2 <- tmp_pos
for(i in 1:length(tmp_pos)) {
  if (tmp_pos[i] > terminus) {
    new_pos2[i] <- max_pos - tmp_pos[i]
  } else {
  }
}

tmp_pos <- new_pos2


data_chrom$midpoint <- tmp_pos 


data_chrom_ordered <- data_chrom[order(data_chrom$midpoint),]

rows_to_split_dat <- vector() #empty vec  to hold the rows that split the dat into 10kb chunks
len_of_genome <- 800000 #rounded up to closest 10,000 so that all of the points are included
chunks <- seq(50000, len_of_genome, 50000)
for (i in chunks) {
  rows <- which(abs(data_chrom_ordered$midpoint-i)==min(abs(data_chrom_ordered$midpoint-i))) # finding the closest number to each 10kb without going over it
  max_row <- max(rows)
#  actual_pos <- data_chrom_ordered$midpoint[max_row]
#  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  rows_to_split_dat <- c(rows_to_split_dat, max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1

#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
data_just_COG <- data_chrom_ordered$COG
#list_dat_sets <- split(data_chrom_ordered, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))
list_dat_sets_just_COG <- split(data_just_COG, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))

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
chunks_pos_zero <- c(chunks, 850000)
chunk_pos <- rep(chunks_pos_zero, each = length(COG_cat))
freq_dat <- cbind(COG_count, COG_cat_col, chunk_pos)

#tick lables
my_labs_mill <- seq(0,850000, by = 50000)


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
  ggtitle(expression(paste(italic("S.meliloti"), " pSymB")))+
  labs(x = "Position in Genome (bp)", y = "Total Number of Genes") +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank() ) +
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
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

pdf("pSymB_COG_bargraph_and_hist_30Jan20.pdf")
grid.arrange(p_bar_top, legend, p_stacked_bar, layout_matrix = lay, widths=c(4,1), heights=c(1.5, 3.5))
#             top="S.meliloti pSymB")
#grid.arrange(p_bar_top, legend, p_stacked_bar, ncol=2, nrow=2, widths=c(4.5,0.5), heights=c(1.5, 3.5), top="S.meliloti Chromosome")
dev.off()


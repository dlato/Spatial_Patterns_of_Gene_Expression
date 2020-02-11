set.seed(1738)
setwd("/home/dlato/Documents/GoldingSummer14/COG_info_17Jan17/Strep/")
path="/home/dlato/Documents/GoldingSummer14/COG_info_17Jan17/Strep/"

###########################
#STREP
###########################

chrom_datafile <- "/home/dlato/Documents/GoldingSummer14/COG_info_17Jan17/Strep/Strep_chrom_COG_binary_data.csv"
data_chrom <- read.csv(chrom_datafile)

#accounting for bidirectionality of replication

#accounting for the bidirectionality in strep and that both
#ends of the chromosomes (bc it is linear) are considered "far"
#from the origin of rep. slightly different than the other bacteria
##max_pos <- max(chrom_gaps_data$tmp_pos)
tmp_pos <- data_chrom$midpoint
new_pos <- tmp_pos
oriC_pos <- 3419363
for(i in 1:length(tmp_pos)) {
      if (tmp_pos[i] >= oriC_pos) {
	      new_pos[i] <- tmp_pos[i] - oriC_pos
  }#if btwn origin and end of genome 
  if (tmp_pos[i] <= oriC_pos) {
          new_pos[i] <- oriC_pos - tmp_pos[i]
    }#if btwn origin and beginning of genome
    if (tmp_pos[i] == oriC_pos) {
	    new_pos[i] <- 0
      }#if equal to origin
}

tmp_pos <- new_pos


data_chrom$midpoint <- tmp_pos 

#logistic regression for each COG
COGs <- levels(data_chrom$COG)

for (cog in COGs) {
  print(cog)
  dat_COG <- data_chrom[which(data_chrom$COG == cog),]
  logistic_reg_model=glm(value~midpoint,family=binomial,dat_COG, control = list(maxit = 1000)) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
  print(summary(logistic_reg_model))
}

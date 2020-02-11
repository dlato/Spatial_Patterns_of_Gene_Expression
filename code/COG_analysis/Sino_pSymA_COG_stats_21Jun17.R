set.seed(1738)
setwd("/home/dlato/Documents/GoldingSummer14/COG_info_17Jan17/Sino/pSymA/")
path="/home/dlato/Documents/GoldingSummer14/COG_info_17Jan17/Sino/pSymA/"

###########################
#SINO pSymA
###########################

chrom_datafile <- "/home/dlato/Documents/GoldingSummer14/COG_info_17Jan17/Sino/pSymA/pSymA_COG_binary_data.csv"
data_chrom <- read.csv(chrom_datafile)

##account for bidirectionality of replication.
##using only the midpoint as the position. not scaling the 
#first scaling things to the origin (if necessary)
tmp_pos <- data_chrom$midpoint
max_pos <- max(data_chrom$midpoint)
new_pos <- tmp_pos
oriC_pos <- 0
terminus <- 816660
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

#logistic regression for each COG
COGs <- levels(data_chrom$COG)

for (cog in COGs) {
  print(cog)
  dat_COG <- data_chrom[which(data_chrom$COG == cog),]
  logistic_reg_model=glm(value~midpoint,family=binomial,dat_COG, control = list(maxit = 1000)) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
  print(summary(logistic_reg_model))
}

############################################
#interactive graphs for gene expression data
############################################
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(strip.text = element_text(size =10)) +
            #theme(plot.title = element_text(hjust = 0.5), 
            theme(plot.title = element_text(), 
                  panel.background = element_rect(fill = "white", colour = NA), 
                  panel.grid.major = element_line(colour = "grey90", size = 0.2), 
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
                  panel.spacing = unit(0.25, "lines"), 
                  axis.text=element_text(size=10),
                  legend.position="top") 
)
options(scipen=10000)
################################################################################################################################
################################################################################################################################
# ECOLI
################################################################################################################################
################################################################################################################################
############################################
#read in data
exp_data <- read.table("ecoli_bidirectional_data.csv", sep = "\t")
#re-scale position
exp_data$tmp_pos <- exp_data$tmp_pos / 1000000
############################################
#lin reg on points
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=exp_data)
summary(lm_all_exp)
############################################
#mytitle <- substitute(italic(bac_name)~replicon, list(bac_name=bac_name, replicon=replicon))
#bac_name <- "E.coli"
mytitle <- substitute(italic(bac_name)~replicon, list(bac_name="E.coli", replicon="Chromosome"))
#mytitle <- substitute(<i>bac_name</i>~replicon, list(bac_name="E.coli", replicon="Chromosome"))
test_title <- "<i>E.coli</i> Chromosome"
#scatter plot
scat <- (exp_data
         %>% highlight_key(~Id)
         %>% ggplot(aes(tmp_pos, Exp, group=Id)) 
         + ggtitle(test_title)
#         + ggtitle(mytitle)
         + geom_point(alpha = 0.5)
#         + labs(title = mytitle)
         + ylab("Expression (CPM)")
         + xlab("Bidirectional Genomic Position (Mbp)")
         + scale_y_log10()
)
#warning because of zero expression values changed on a log scale
scat

#interactive
inter_scat <- (scat
  %>% ggplotly(tooltip="group")
       #%>%  layout(title = "<i>E.coli</i> Chromosome",
       %>%  layout(
                   margin = list(b = 50, l = 50)) # to fully display the x and y axis labels
  %>% highlight(on="plotly_hover",
                off=NULL,
                color="#0083FF")
)
inter_scat
##save interactive
htmlwidgets::saveWidget(as_widget(inter_scat), "ecoli_gene_exp.html")

################################################################################################################################
################################################################################################################################
# BASS
################################################################################################################################
################################################################################################################################
#read in data
exp_data <- read.table("bass_bidirectional_data.csv", sep = "\t")
#re-scale position
exp_data$tmp_pos <- exp_data$tmp_pos / 1000000
############################################
#lin reg on points
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=exp_data)
summary(lm_all_exp)
############################################
test_title <- "<i>B.subtilis</i> Chromosome"
#scatter plot
scat <- (exp_data
         %>% highlight_key(~Id)
         %>% ggplot(aes(tmp_pos, Exp, group=Id)) 
         + ggtitle(test_title)
         #         + ggtitle(mytitle)
         + geom_point(alpha = 0.5)
         #         + labs(title = mytitle)
         + ylab("Expression (CPM)")
         + xlab("Bidirectional Genomic Position (Mbp)")
         + scale_y_log10()
)
#warning because of zero expression values changed on a log scale
scat

#interactive
inter_scat <- (scat
               %>% ggplotly(tooltip="group")
               #%>%  layout(title = "<i>E.coli</i> Chromosome",
               %>%  layout(
                 margin = list(b = 50, l = 50)) # to fully display the x and y axis labels
               %>% highlight(on="plotly_hover",
                             off=NULL,
                             color="#0083FF")
)
inter_scat
##save interactive
htmlwidgets::saveWidget(as_widget(inter_scat), "bsubtilis_gene_exp.html")

################################################################################################################################
################################################################################################################################
# Strep
################################################################################################################################
################################################################################################################################
#read in data
exp_data <- read.table("strep_bidirectional_data.csv", sep = "\t")
#re-scale position
exp_data$tmp_pos <- exp_data$tmp_pos / 1000000
############################################
#lin reg on points
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=exp_data)
summary(lm_all_exp)
############################################
test_title <- "<i>Streptomyces</i> Chromosome"
#scatter plot
scat <- (exp_data
         %>% highlight_key(~Id)
         %>% ggplot(aes(tmp_pos, Exp, group=Id)) 
         + ggtitle(test_title)
         #         + ggtitle(mytitle)
         + geom_point(alpha = 0.5)
         #         + labs(title = mytitle)
         + ylab("Expression (CPM)")
         + xlab("Genomic Position From The Origin of Replication (Mbp)")
         + scale_y_log10()
)
#warning because of zero expression values changed on a log scale
scat

#interactive
inter_scat <- (scat
               %>% ggplotly(tooltip="group")
               %>%  layout(
                 margin = list(b = 50, l = 50)) # to fully display the x and y axis labels
               %>% highlight(on="plotly_hover",
                             off=NULL,
                             color="#0083FF")
)
inter_scat
##save interactive
htmlwidgets::saveWidget(as_widget(inter_scat), "streptomyces_gene_exp.html")

################################################################################################################################
################################################################################################################################
# SinoC
################################################################################################################################
################################################################################################################################
#read in data
exp_data <- read.table("sinoC_bidirectional_data.csv", sep = "\t")
#re-scale position
exp_data$tmp_pos <- exp_data$tmp_pos / 1000000
############################################
#lin reg on points
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=exp_data)
summary(lm_all_exp)
############################################
test_title <- "<i>S.meliloti</i> Chromosome"
#scatter plot
scat <- (exp_data
         %>% highlight_key(~Id)
         %>% ggplot(aes(tmp_pos, Exp, group=Id)) 
         + ggtitle(test_title)
         #         + ggtitle(mytitle)
         + geom_point(alpha = 0.5)
         #         + labs(title = mytitle)
         + ylab("Expression (CPM)")
         + xlab("Bidirectional Genomic Position (Mbp)")
         + scale_y_log10()
)
#warning because of zero expression values changed on a log scale
scat

#interactive
inter_scat <- (scat
               %>% ggplotly(tooltip="group")
               #%>%  layout(title = "<i>E.coli</i> Chromosome",
               %>%  layout(
                 margin = list(b = 50, l = 50)) # to fully display the x and y axis labels
               %>% highlight(on="plotly_hover",
                             off=NULL,
                             color="#0083FF")
)
inter_scat
##save interactive
htmlwidgets::saveWidget(as_widget(inter_scat), "smeliloti_chromosome_gene_exp.html")
################################################################################################################################
################################################################################################################################
# pSymA
################################################################################################################################
################################################################################################################################
#read in data
exp_data <- read.table("sinoC_bidirectional_data.csv", sep = "\t")
#re-scale position
exp_data$tmp_pos <- exp_data$tmp_pos / 1000000
############################################
#lin reg on points
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=exp_data)
summary(lm_all_exp)
############################################
test_title <- "<i>S.meliloti</i> pSymA"
#scatter plot
scat <- (exp_data
         %>% highlight_key(~Id)
         %>% ggplot(aes(tmp_pos, Exp, group=Id)) 
         + ggtitle(test_title)
         #         + ggtitle(mytitle)
         + geom_point(alpha = 0.5)
         #         + labs(title = mytitle)
         + ylab("Expression (CPM)")
         + xlab("Bidirectional Genomic Position (Mbp)")
         + scale_y_log10()
)
#warning because of zero expression values changed on a log scale
scat

#interactive
inter_scat <- (scat
               %>% ggplotly(tooltip="group")
               #%>%  layout(title = "<i>E.coli</i> Chromosome",
               %>%  layout(
                 margin = list(b = 50, l = 50)) # to fully display the x and y axis labels
               %>% highlight(on="plotly_hover",
                             off=NULL,
                             color="#0083FF")
)
inter_scat
##save interactive
htmlwidgets::saveWidget(as_widget(inter_scat), "smeliloti_pSymA_gene_exp.html")

################################################################################################################################
################################################################################################################################
# pSymB
################################################################################################################################
################################################################################################################################
#read in data
exp_data <- read.table("sinoC_bidirectional_data.csv", sep = "\t")
#re-scale position
exp_data$tmp_pos <- exp_data$tmp_pos / 1000000
############################################
#lin reg on points
lm_all_exp <- lm(abs(Exp) ~ tmp_pos, data=exp_data)
summary(lm_all_exp)
############################################
test_title <- "<i>S.meliloti</i> pSymB"
#scatter plot
scat <- (exp_data
         %>% highlight_key(~Id)
         %>% ggplot(aes(tmp_pos, Exp, group=Id)) 
         + ggtitle(test_title)
         #         + ggtitle(mytitle)
         + geom_point(alpha = 0.5)
         #         + labs(title = mytitle)
         + ylab("Expression (CPM)")
         + xlab("Bidirectional Genomic Position (Mbp)")
         + scale_y_log10()
)
#warning because of zero expression values changed on a log scale
scat

#interactive
inter_scat <- (scat
               %>% ggplotly(tooltip="group")
               #%>%  layout(title = "<i>E.coli</i> Chromosome",
               %>%  layout(
                 margin = list(b = 50, l = 50)) # to fully display the x and y axis labels
               %>% highlight(on="plotly_hover",
                             off=NULL,
                             color="#0083FF")
)
inter_scat
##save interactive
htmlwidgets::saveWidget(as_widget(inter_scat), "smeliloti_pSymB_gene_exp.html")
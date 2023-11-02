#############################################################################
############### Extended Data Figure 6 ######################################

### Jana Ordon
### 230314

rm(list=ls())

# set working directory, change it accordingly
setwd("//fs-bio.mpipz.mpg.de/biodata/dep_psl/grp_psl/Jana/Results/Barcoding/WCS358_proof/liquid_cultures_PQQ_Rep2")

library(tidyverse)
library(ggrepel)
library(ggplot2)
library(Rcpp)

# open file
Kin <- read.delim("OD600_WCS358_PQQ_extended_data_fig_7.txt", stringsAsFactors = F)

# add combined column
Kin$treatment <- paste0(Kin$strain,"_", Kin$glucose)
Kin$treatment_PQQ <- paste0(Kin$treatment,"_", Kin$PQQ)


#set color code
colors <- c(  "wt_none_none" = "#626693",
              "wt_110mM_glucose_none" = "#414462",
              "wt_110mM_glucose_PQQ" = "#A9ACC6",
              "pqqF_none_none" = "#e07a5f",
              "pqqF_110mM_glucose_none" = "#8F5D5D",
              "pqqF_110mM_glucose_PQQ" = "#EAB69F",
              "cyoB_none_none" = "#81b29a",
              "cyoB_110mM_glucose_none" = "#5F797B",
              "cyoB_110mM_glucose_PQQ" = "#BABF95")

# plot XVM2_succinate
ggplot(Kin_XVM2_succinate) + 
  geom_ribbon(aes(x = Time..h., y = OD_mean,  fill=treatment_PQQ, ymin = OD_mean - OD_se/sqrt(3), ymax = OD_mean + OD_se/sqrt(3)), alpha=0.2) +
  geom_line(aes(x = Time..h., y = OD_mean,  color=treatment_PQQ)) +
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black"),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(color="Black", hjust = 1))+
  ylab("OD600")+
  xlab("Time [h]")+
  ggsave(paste("Kin_WCS358_XVM2_succinate.pdf", sep=""), width=7, height=5)


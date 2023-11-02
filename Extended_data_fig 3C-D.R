#############################################################################
############### Extended Data Fig 3c,d ####################################################


### Jana Ordon
### 230926
### Correlate PCR bias genomically integrated MoBacTag with PCR bias for genomic sequence 


rm(list=ls())

library("ggplot2")
library("tidyverse")
library("rcompanion")
library("ggpubr")
library("scales")
library("grid")
library("vegan")
library("ggrepel")

# set working directory, change it accordingly
setwd("//fs-bio.mpipz.mpg.de/biodata/dep_psl/grp_psl/Jana/Results/Barcoding/Primer_bias_qPCR")


# load data
qPCR <- read.delim("Extended_Data_Fig_3cd.txt", stringsAsFactors = F)

qPCR$primer_strain <- paste0(qPCR$primer, "_", qPCR$strain)

# set colors
color <- c( "MoBacTag_BC_wt" = "#3d405b",
            "plant_ITS_wt" = "#3d405b",
            "MoBacTag_BC_pqqF" = "#e07a5f",
            "plant_ITS_pqqF" = "#e07a5f",
            "MoBacTag_BC_cyoB" = "#81b29a",
            "plant_ITS_cyoB" = "#81b29a",
            "WCS358_qPCR_1_BC_wt"="#6C719D",
            "WCS358_qPCR_1_plant_ITS_wt"="#6C719D",
            "WCS358_qPCR_1_BC_pqqF"="#EFB8A9",
            "WCS358_qPCR_1_plant_ITS_pqqF"="#EFB8A9",
            "WCS358_qPCR_1_BC_cyoB"="#C0D8CC",
            "WCS358_qPCR_1_plant_ITS_cyoB"="#C0D8CC")

# filter
idx <- qPCR$primer %in% c("MoBacTag_BC", "WCS358_qPCR_1_BC") &
  TRUE
qPCR_BC <- qPCR[idx, ]

idx <- qPCR$primer %in% c("plant_ITS", "WCS358_qPCR_1_plant_ITS") &
  TRUE
qPCR_pITS <- qPCR[idx, ]



######### plot #############

# BC
ggplot((data = qPCR_BC), aes(x=dilution, y=ct, group=primer_strain, color=primer_strain))+
  geom_line(stat="smooth", method='lm', formula= y~x, se=FALSE, alpha=0.4) +
  geom_point(aes(color=primer_strain), alpha = 0.8)+
  scale_x_continuous(trans='log10') +
  scale_color_manual(values = color) +
  stat_regline_equation(size = 2, label.x = log10(0.01)) +
  stat_cor(size = 2, label.x = log10(0.1)) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  ggsave(paste("qPCR_correlation_BC.pdf", sep=""), width=3.2, height=3.2)

# pITS
ggplot((data = qPCR_pITS), aes(x=dilution, y=ct, group=primer_strain, color=primer_strain))+
  geom_line(stat="smooth", method='lm', formula= y~x, se=FALSE, alpha=0.4) +
  geom_point(aes(color=primer_strain), alpha = 0.8)+
  scale_x_continuous(trans='log10') +
  scale_color_manual(values = color) +
  stat_regline_equation(size = 2, label.x = log10(0.01)) +
  stat_cor(size = 2, label.x = log10(0.1)) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.size = "none") +
  ggsave(paste("qPCR_correlation_pITS.pdf", sep=""), width=3.2, height=3.2)


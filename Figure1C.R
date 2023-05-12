#############################################################################
############### Figure 1 c ##################################################

### Jana Ordon
### 220411
### Heatmap frequency primer per position


options(warn=-1)

# cleanup

rm(list=ls())


library("ggplot2")
library("tidyverse")


# set working directory, change it accordingly
setwd("C:/Users/Jana Ordon/Dropbox/R/Barcodes_temp")


# load data
BC_per_position <- read.delim("frequency_primer_per_position.txt", stringsAsFactors = F)

# plot
ggplot(BC_per_position, aes(fill = frequency_percent, x = position, y = ID)) +
  geom_tile() +
  scale_fill_gradient(low = "#386480", high = "#B6E0E2") +
  ggsave(paste("BC_frequency_per_position_heatmap.pdf", sep = ""), width = 6, height =4)
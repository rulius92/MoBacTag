#############################################################################
############### Extended Data Figure 7 ######################################



### Jana Ordon
### 230904

rm(list=ls())

options(warn=-1)


# set working directory, change it accordingly
setwd("C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_soil/")

library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("dplyr")
library("RColorBrewer")

# load plotting functions
source("C:/Users/Jana Ordon/Dropbox/R/MiSeq/plotting_functions.R")
source("C:/Users/Jana Ordon/Dropbox/R/MiSeq/cpcoa.func.R")

# directories

results.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_soil/results/"
data.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_soil/data/"
figures.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_soil/figures/"

# files

design.file <- paste(data.dir, "MoBacTag_soil_design.txt", sep="")
otu_table.file <- paste(results.dir, "MoBacTag_soil_ASV_table_.txt", sep="")
taxonomy.file <- paste(data.dir, "taxonomy.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

############# Extended Data Figure 7a ######################################
# remove spike
idx <- rownames(otu_table) %in% c("Spike_5f2874b1231b8117f0abec2e27a4b20c")
spike <- otu_table[idx, ]
otu_table_rel <- otu_table[!idx,  ]

# normalize otu tables

design$depth <- colSums(otu_table)
otu_table_rel <- apply(otu_table_rel, 2, function(x) x/sum(x))


### compute beta diversity
colors <- data.frame(group=c("Cas", "Golm","Jify"),
                     color=c(c_black, c_dark_brown, c_orange))

########## PCoA Bray-Curtis
########## subset samples of interest

### Compartment effect

idx <- design$soil %in% c("Cas", "Golm", "Jify") &
  design$strain %in% c("R13DBC_R13CBC")
TRUE

design <- design[idx, ]
otu_table_rel <- otu_table_rel[, idx]

otu_table_rel[is.na(otu_table_rel)] <- 0

bray_curtis <- vegdist(t(otu_table_rel), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

colors_pcoa <- colors[colors$group %in% points$soil, ]
points$soil <- factor(points$soil, levels=colors_pcoa$group)

# plot PCo 1 and 2

ggplot(points, aes(x=x, y=y, color=soil, shape=compartment)) +
  geom_point(alpha=.6, size=2.0) +
  scale_colour_manual(values=as.character(colors_pcoa$color)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_rel_16S_soils.pdf", sep=""), width=7, height=3)

############# Extended Data Figure 7b-d ######################################

# spike 16S normalized otu tables
idx <- rownames(otu_table) == "Spike_5f2874b1231b8117f0abec2e27a4b20c"
Spike <- otu_table[idx, ]
Spike <- unlist(Spike)
design <- cbind(design, Spike)
otu_table <- otu_table[!idx,  ]
otu_table_norm <- sweep(otu_table, 2, Spike, `/`)


# dataframe with information plus RA
df_16S <- melt(as.matrix(otu_table_norm))
colnames(df_16S) <- c("bac", "SampleID", "abundance")

#match sampleID with the corresponding WCS4174_deriv, compartment
df_16S$compartment <- design$compartment[match(df_16S$SampleID, design$SampleID)]
df_16S$soil <- design$soil[match(df_16S$SampleID, design$SampleID)]
df_16S$bio_replicate <- design$bio_replicate[match(df_16S$SampleID, design$SampleID)]
df_16S$tech_replicate <- design$tech_replicate[match(df_16S$SampleID, design$SampleID)]
df_16S$strain <- design$strain[match(df_16S$SampleID, design$SampleID)]

df_16S_t<-left_join(df_16S, taxonomy)

##subset samples
idx <- df_16S_t$strain %in% c("R13DBC_R13CBC") &
  df_16S_t$compartment %in% c("Root") &
  T
df_wt_root <- df_16S_t[idx, ]

idx <- df_16S_t$strain %in% c("R13DBC_R13CBC") &
  df_16S_t$compartment %in% c("Rhizosphere") &
  T
df_wt_rhizosphere <- df_16S_t[idx, ]

idx <- df_16S_t$strain %in% c("R13DBC_R13CBC") &
  df_16S_t$compartment %in% c("Soil") &
  T
df_wt_soil <- df_16S_t[idx, ]

nb.cols=39
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    axis.text.x=element_text(colour="black", size=14, hjust = 1,vjust = 1,angle=45),
                    #use none to remove legend, can choose top, right, left
                    legend.position="right",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans", size=15))

##plot root 16S
ggplot(df_wt_root, aes(x=soil, y=abundance, fill = phylum)) +
  geom_bar(position="stack", stat="summary", width = 0.8) +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  main_theme +
  labs(y="abundance")+
  ggsave(paste(figures.dir, "stacked_bar_root_wt_abs_phylum.pdf", sep=""), width=7.5, height=10)

ggplot(df_wt_rhizosphere, aes(x=soil, y=abundance, fill = phylum)) +
  geom_bar(position="stack", stat="summary", width = 0.8) +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  main_theme +
  labs(y="abundance")+
  ggsave(paste(figures.dir, "stacked_bar_rhizosphere_wt_abs_phylum.pdf", sep=""), width=7.5, height=2.8)

ggplot(df_wt_soil, aes(x=soil, y=abundance, fill = phylum)) +
  geom_bar(position="stack", stat="summary", width = 0.8) +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  main_theme +
  labs(y="abundance")+
  ggsave(paste(figures.dir, "stacked_bar_soil_wt_abs_phylum.pdf", sep=""), width=7.5, height=2.8)

############# Extended Data Figure 7e-g ######################################

# error correct BCs
idx <- df_16S_t$bac %in% c("R13DBC_5c90df4ca6a75d0f740bcd063fb13e88") &
  df_16S_t$strain %in% c("R13DBC_R13CBC") &
  T
df_R13Dbc <- df_16S_t[idx,  ]
df_R13Dbc$abundance <- df_R13Dbc$abundance / 1.2

idx <- df_16S_t$bac %in% c("R13CBC_31a93e5d5755852647bfbcdb7f55ae9d") &
  df_16S_t$strain %in% c("R13DBC_R13CBC") &
  T
df_R13Cbc <- df_16S_t[idx,  ]
df_R13Cbc$abundance <- df_R13Cbc$abundance / 1.2

df_BCs_cor <- rbind(df_R13Cbc, df_R13Dbc)

### subset samples
#BCs
idx <- df_BCs_cor$compartment %in% c("Root") &
  T
df_root_BC <- df_BCs_cor[idx, ]

idx <- df_BCs_cor$compartment %in% c("Rhizosphere") &
  T
df_rhizosphere_BC <- df_BCs_cor[idx, ]

idx <- df_BCs_cor$compartment %in% c("Soil") &
  T
df_soil_BC <- df_BCs_cor[idx, ]

# 16S
idx <- df_16S_t$bac %in% c("R13D_335241f9cdeb4a16e3f395cb7809d30e", 
                           "R13C_9282a921c4b198bda0d18f901d3ac07a") &
  df_16S_t$strain %in% c("R13DBC_R13CBC") &
  df_16S_t$compartment %in% c("Root") &
  T
df_root_16S <- df_16S_t[idx, ]

idx <- df_16S_t$bac %in% c("R13D_335241f9cdeb4a16e3f395cb7809d30e", 
                           "R13C_9282a921c4b198bda0d18f901d3ac07a") &
  df_16S_t$strain %in% c("R13DBC_R13CBC") &
  df_16S_t$compartment %in% c("Rhizosphere") &
  T
df_rhizosphere_16S <- df_16S_t[idx, ]

idx <- df_16S_t$bac %in% c("R13D_335241f9cdeb4a16e3f395cb7809d30e", 
                           "R13C_9282a921c4b198bda0d18f901d3ac07a") &
  df_16S_t$strain %in% c("R13DBC_R13CBC") &
  df_16S_t$compartment %in% c("Soil") &
  T
df_soil_16S <- df_16S_t[idx, ]

# define colors
colors <- c(  "R13D_335241f9cdeb4a16e3f395cb7809d30e" = "#B43E8F",
              "R13DBC_5c90df4ca6a75d0f740bcd063fb13e88" = "#B43E8F",
              "R13C_9282a921c4b198bda0d18f901d3ac07a" = "#5A00A3",
              "R13CBC_31a93e5d5755852647bfbcdb7f55ae9d" = "#5A00A3")

##plot BC
ggplot(df_root_BC, aes(x=soil, y=abundance, fill = bac)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA, size=.4, width=.85) +
  geom_jitter(aes(color=bac, shape = bio_replicate), position=position_jitterdodge(0.2), size=.9, alpha=0.9) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=14))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.position="top") +
  #ylim (0, 0.23) +
  labs(y="spike normalized")+
  ggsave(paste(figures.dir, "box_R13C_R13D_root_BC_abs.pdf", sep=""), width=2, height=4)

ggplot(df_rhizosphere_BC, aes(x=soil, y=abundance, fill = bac)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA, size=.4, width=.85) +
  geom_jitter(aes(color=bac, shape = bio_replicate), position=position_jitterdodge(0.2), size=.9, alpha=0.9) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=14))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.position="top") +
  #ylim (0, 0.23) +
  labs(y="spike normalized")+
  ggsave(paste(figures.dir, "box_R13C_R13D_rhizosphere_BC_abs.pdf", sep=""), width=2, height=4)

ggplot(df_soil_BC, aes(x=soil, y=abundance, fill = bac)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA, size=.4, width=.85) +
  geom_jitter(aes(color=bac, shape = bio_replicate), position=position_jitterdodge(0.2), size=.9, alpha=0.9) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=14))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.position="top") +
  #ylim (0, 0.23) +
  labs(y="spike normalized")+
  ggsave(paste(figures.dir, "box_R13C_R13D_soil_BC_abs.pdf", sep=""), width=2, height=4)


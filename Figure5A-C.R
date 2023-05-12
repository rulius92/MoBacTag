#############################################################################
############### Figure 5 ####################################################


############## Figure 5 a - c ###############################################

### Jana Ordon
### 221208
### stacked bar plot 16 read counts, different normalization


rm(list=ls())

# set working directory, change it accordingly
setwd("C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/")


# set directories
results.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/results/"
data.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/data/"
figures.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/figures/"

library("EnvStats")
library("RColorBrewer")
library("viridis")

# files
design.file <- paste(data.dir, "MoBacTag_V2_design.txt", sep="")
otu_table.file <- paste(results.dir, "Roots-16S_merged_asv_table.txt", sep="")
otu_table_plant_file <- paste(results.dir, "Roots-pITS_merged_asv_table.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
otu_table_pITS <- read.table(otu_table_plant_file, sep="\t", header=T, check.names=F)

############
# re-order data matrices 
idx <- design$SampleID %in% colnames(otu_table_pITS)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table_pITS))
otu_table <- otu_table[, idx]
otu_table_pITS <- otu_table_pITS[, idx]

# remove all strain except for 
idx <- rownames(otu_table) %in% c("Root101", "Root123D2", "Root1310", "Root131", "Root142",
                                  "Root265", "Root480", "Root61", "Root685", "Root68", "Root695",
                                  "Root77", "Root83", "Root935", "pBCC084_Spike2_bc72") 
otu_table <- otu_table[idx, ]

idx <- rownames(otu_table_pITS) %in% c("pITS_Arabidopsis", "pBCC084_Spike2_bc72") 
otu_table_pITS <- otu_table_pITS[idx, ]

# spike 16S normalized otu tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC084_Spike2_bc72"
Spike <- otu_table[idx, ]
Spike <- unlist(Spike)
design <- cbind(design, Spike)
otu_table <- otu_table[!idx,  ]
otu_table_16Snorm <- sweep(otu_table, 2, Spike, `/`)

# spike pITS normalized Arabidopsis ITS reads
idx <- rownames(otu_table_pITS) == "pBCC084_Spike2_bc72"
Spike_pITS <- otu_table_pITS[idx, ]
Spike_pITS <- unlist(Spike_pITS)
design <- cbind(design, Spike_pITS)
otu_table_pITS <- otu_table_pITS[!idx,  ]
otu_table_pITSnorm <- sweep(otu_table_pITS, 2, Spike_pITS, `/`)

# 16Sreads normalized by Arabidopsis ITS reads
At_ITS <- unlist(otu_table_pITSnorm)
otu_table_16S_ITSnorm <- sweep(otu_table_16Snorm, 2, At_ITS, `/`)

# rel abundance otu tables
otu_table_rel <- apply(otu_table, 2, function(x) x/sum(x)*100)

# dataframe with information plus RA
df_16S_norm <- melt(as.matrix(otu_table_16Snorm))
colnames(df_16S_norm) <- c("strain", "SampleID", "abundance")
df_16S_norm$type <- "16S_norm"

df_rel <- melt(as.matrix(otu_table_rel))
colnames(df_rel) <- c("strain", "SampleID", "abundance")
df_rel$type <- "rel"

df_16S_pITS_norm <- melt(as.matrix(otu_table_16S_ITSnorm))
colnames(df_16S_pITS_norm) <- c("strain", "SampleID", "abundance")
df_16S_pITS_norm$type <- "16S_pITS_norm"

#match sampleID with the corresponding WCS4174_deriv, compartment
df_rel$deriv <- design$WCS417r_deriv[match(df_rel$SampleID, design$SampleID)]
df_16S_norm$deriv <- design$WCS417r_deriv[match(df_16S_norm$SampleID, design$SampleID)]
df_16S_pITS_norm$deriv <- design$WCS417r_deriv[match(df_16S_pITS_norm$SampleID, design$SampleID)]
df_rel$compartment <- design$compartment[match(df_rel$SampleID, design$SampleID)]
df_16S_norm$compartment <- design$compartment[match(df_16S_norm$SampleID, design$SampleID)]
df_16S_pITS_norm$compartment <- design$compartment[match(df_16S_pITS_norm$SampleID, design$SampleID)]

# add combined column
df_rel$comb <- paste0(df_rel$type, "_", df_rel$deriv)
df_16S_norm$comb <- paste0(df_16S_norm$type, "_", df_16S_norm$deriv)
df_16S_pITS_norm$comb <- paste0(df_16S_pITS_norm$type, "_", df_16S_pITS_norm$deriv)


##subset samples
idx <- df_rel$deriv %in% c("wt", "pqqF", "cyoB", "wt+pqqF", 
                           "wt+pqqF+cyoB") &
  T
df_rel <- df_rel[idx, ]
df_rel$deriv<- factor(df_rel$deriv, levels= c("wt", "pqqF", "cyoB", "wt+pqqF",
                                              "wt+pqqF+cyoB"))

idx <- df_16S_norm$deriv %in% c("wt", "pqqF", "cyoB", "wt+pqqF",
                                "wt+pqqF+cyoB") &
  T
df_16S_norm <- df_16S_norm[idx, ]
df_16S_norm$deriv<- factor(df_16S_norm$deriv, levels= c("wt", "pqqF", "cyoB",
                                                        "wt+pqqF",
                                                        "wt+pqqF+cyoB"))

idx <- df_16S_pITS_norm$deriv %in% c("wt", "pqqF", "cyoB", "wt+pqqF",
                                     "wt+pqqF+cyoB") &
  T
df_16S_pITS_norm <- df_16S_pITS_norm[idx, ]
df_16S_pITS_norm$deriv<- factor(df_16S_pITS_norm$deriv, levels= c("wt", "pqqF", "cyoB",
                                                                  "wt+pqqF", 
                                                                  "wt+pqqF+cyoB"))

write.table(df_rel, "df_rel_abundance_A.txt")
write.table(df_16S_norm, "df_16S_norm_B.txt")
write.table(df_16S_pITS_norm, "df_plant_ITS_norm_C.txt")

#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    axis.text.x=element_text(colour="black", size=14, hjust = 1,vjust = 1,angle=45),
                    legend.position="top",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans", size=15))

nb.cols=16
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

##plot root 16S
ggplot(df_rel, aes(x=deriv, y=abundance, fill = strain)) +
  geom_bar(position="stack", stat="summary", width = 0.8) +
  scale_fill_manual(values = mycolors) +
  scale_x_discrete(c("wt", "pqqF", "cyoB", "wt+pqqF",
                     "wt+pqqF+cyoB"))+
  theme_bw() +
  main_theme +
  labs(y="abundance")+
  ggsave(paste(figures.dir, "stacked_bar_root_rel_noWCS_v2.pdf", sep=""), width=2.3, height=5)


ggplot(df_16S_norm, aes(x=deriv, y=abundance, fill = strain)) +
  geom_bar(position="stack", stat="summary", width = 0.8) +
  scale_fill_manual(values = mycolors) +
  scale_x_discrete(c("wt", "pqqF", "cyoB", "wt+pqqF",
                     "wt+pqqF+cyoB"))+
  theme_bw() +
  main_theme +
  labs(y="abundance")+
  ggsave(paste(figures.dir, "stacked_bar_root_16S_norm_noWCS_v2.pdf", sep=""), width=2.3, height=5)

ggplot(df_16S_pITS_norm, aes(x=deriv, y=abundance, fill = strain)) +
  geom_bar(position="stack", stat="summary", width = 0.8) +
  scale_fill_manual(values = mycolors) +
  scale_x_discrete(c("wt", "pqqF", "cyoB", "wt+pqqF",
                     "wt+pqqF+cyoB"))+
  theme_bw() +
  main_theme +
  labs(y="abundance")+
  ggsave(paste(figures.dir, "stacked_bar_root_16S_pITS_norm_noWCS_v2.pdf", sep=""), width=2.3, height=5)
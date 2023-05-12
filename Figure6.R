#############################################################################
############### Figure 6 ####################################################

# set working directory, change it accordingly
setwd("C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/")

# set directories
results.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/results/"
data.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/data/"
figures.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/figures/"

library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("tidyverse")
library("rstatix")

############## Figure 6 a - b ###############################################

otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# remove all strain except for 
idx <- rownames(otu_table) %in% c("Root101", "Root123D2", "Root1310", "Root131", "Root142",
                                    "Root265", "Root480", "Root61", "Root685", "Root68", "Root695",
                                    "Root77", "Root83", "Root935", "WCS358", "pBCC0023") 
otu_table <- otu_table[idx, ]
  
# spike 16S normalized otu tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC0023"
Spike <- otu_table[idx, ]
Spike <- unlist(Spike)
design <- cbind(design, Spike)
otu_table <- otu_table[!idx,  ]
otu_table_16Snorm <- sweep(otu_table, 2, Spike, `/`)
  
# rel abundance otu tables
otu_table_rel <- apply(otu_table, 2, function(x) x/sum(x)*100)
  
# dataframe with information plus RA
df_16S_norm <- melt(as.matrix(otu_table_16Snorm))
colnames(df_16S_norm) <- c("strain", "SampleID", "abundance")
df_16S_norm$type <- "16S_norm"
  
df_rel <- melt(as.matrix(otu_table_rel))
colnames(df_rel) <- c("strain", "SampleID", "abundance")
df_rel$type <- "rel"
  
#match sampleID with the corresponding WCS4174_deriv, compartment
df_rel$deriv <- design$WCS417r_deriv[match(df_rel$SampleID, design$SampleID)]
df_16S_norm$deriv <- design$WCS417r_deriv[match(df_16S_norm$SampleID, design$SampleID)]
df_rel$compartment <- design$compartment[match(df_rel$SampleID, design$SampleID)]
df_16S_norm$compartment <- design$compartment[match(df_16S_norm$SampleID, design$SampleID)]
  
# add combined column
df_rel$comb <- paste0(df_rel$compartment, "_", df_rel$deriv)
df_16S_norm$comb <- paste0(df_16S_norm$compartment, "_", df_16S_norm$deriv)
  
# only keep WCS358 abundances
idx <- df_rel$strain %in% c("WCS358") &
  T
df_rel <- df_rel[idx, ]
  
idx <- df_16S_norm$strain %in% c("WCS358") &
  T
df_16S_norm <- df_16S_norm[idx, ]
  
##subset samples
idx <- df_rel$deriv %in% c("wt", "pqqF", "cyoB") &
  df_rel$compartment %in% c("root", "matrix") &
  T
df_rel <- df_rel[idx, ]
df_rel$deriv<- factor(df_rel$deriv, levels= c("wt", "pqqF", "cyoB"))
  
idx <- df_16S_norm$deriv %in% c("wt", "pqqF", "cyoB") &
  df_16S_norm$compartment %in% c("root", "matrix") &
  T
df_16S_norm <- df_16S_norm[idx, ]
df_16S_norm$deriv<- factor(df_16S_norm$deriv, levels= c("wt", "pqqF", "cyoB"))
  
#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                      axis.line.y=element_line(color="black"),
                      axis.ticks=element_line(color="black"),
                      axis.text=element_text(colour="black", size=11),
                      axis.text.x=element_text(colour="black", size=11, hjust = 1,vjust = 1,angle=45),
                      #use none to remove legend, can choose top, right, left
                      legend.position="top",
                      legend.title=element_blank(),
                      legend.background=element_blank(),
                      legend.key=element_blank(),
                      text=element_text(family="sans", size=11))
  
#plot root 16S
ggplot(df_rel, aes(x=comb, y=abundance, fill = compartment)) +
    geom_bar(position = "dodge", stat="summary", width = 0.7) +
    scale_x_discrete(limits = c("matrix_wt", "matrix_pqqF", "matrix_cyoB"))+
    geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=1.0) +
    geom_jitter(aes(x=comb), position = position_dodge(width = 1), fill = "black", size = 0.7, alpha = 0.5) +
    theme_bw() +
    main_theme +
    labs(y="abundance")+
    ggsave(paste(figures.dir, "stacked_bar_rel_WCS358_matrix.pdf", sep=""), width=1.8, height=3.8)
  
ggplot(df_rel, aes(x=comb, y=abundance, fill = compartment)) +
    geom_bar(position = "dodge", stat="summary", width = 0.7) +
    scale_x_discrete(limits = c("root_wt", "root_pqqF", "root_cyoB"))+
    geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=1.0) +
    geom_jitter(aes(x=comb), position = position_dodge(width = 1), fill = "black", size = 0.7, alpha = 0.5) +
    theme_bw() +
    main_theme +
    labs(y="abundance")+
    ylim(0, 12) +
    ggsave(paste(figures.dir, "stacked_bar_rel_WCS358_root.pdf", sep=""), width=1.8, height=3.8)
  
  
ggplot(df_16S_norm, aes(x=comb, y=abundance, fill = compartment)) +
    geom_bar(position = "dodge", stat="summary", width = 0.7) +
    scale_x_discrete(limits = c("matrix_wt", "matrix_pqqF", "matrix_cyoB"))+
    geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=1.0) +
    geom_jitter(aes(x=comb), position = position_dodge(width = 1), fill = "black", size = 0.7, alpha = 0.5) +
    theme_bw() +
    main_theme +
    labs(y="abundance")+
    ggsave(paste(figures.dir, "stacked_bar_spike_norm_WCS358_matrix.pdf", sep=""), width=1.8, height=3.8)
  
ggplot(df_16S_norm, aes(x=comb, y=abundance, fill = compartment)) +
    geom_bar(position = "dodge", stat="summary", width = 0.7) +
    scale_x_discrete(limits = c("root_wt", "root_pqqF","root_cyoB"))+
    geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=1.0) +
    geom_jitter(aes(x=comb), position = position_dodge(width = 1), fill = "black", size = 0.7, alpha = 0.5) +
    theme_bw() +
    main_theme +
    labs(y="abundance")+
    ylim(0, 0.25) +
    ggsave(paste(figures.dir, "stacked_bar_spike_norm_WCS358_root.pdf", sep=""), width=1.8, height=3.8)

# stats #
df_rel %>%
  dunn_test(data =., abundance ~ comb) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  #filter(., p.adj<=0.05) %>%
  filter()->WCS358_stats

############## Figure 6 c - d ###############################################

otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# spike 16S normalized otu tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC0023"
Spike <- otu_table[idx, ]
Spike <- unlist(Spike)
design <- cbind(design, Spike)
otu_table <- otu_table[!idx,  ]
otu_table_norm <- sweep(otu_table, 2, Spike, `/`)

# remove all strain except for 
idx <- rownames(otu_table_norm) %in% c("WCS358", "WCS358_pBC92", "WC358_pqqF_pBC93", 
                                         "WC358_cyo0_pBC163") 
otu_table_BCs <- otu_table_norm[idx,  ]


# dataframe with information plus RA
df_BCs <- melt(as.matrix(otu_table_BCs))
colnames(df_BCs) <- c("strain", "SampleID", "RA")

#match sampleID with the corresponding WCS4174_deriv, compartment
df_BCs$compartment <- design$compartment[match(df_BCs$SampleID, design$SampleID)]
df_BCs$deriv <- design$WCS417r_deriv[match(df_BCs$SampleID, design$SampleID)]
df_BCs$strain_deriv <- paste0(df_BCs$strain, "_", df_BCs$deriv)

# error correct BCs
idx <- df_BCs$strain %in% c("WCS358_pBC92")
df_BC_wt <- df_BCs[idx,  ]
df_BC_wt$RA_norm <- df_BC_wt$RA / 5.1

idx <- df_BCs$strain %in% c("WC358_pqqF_pBC93")
df_BC_pqqF <- df_BCs[idx,  ]
df_BC_pqqF$RA_norm <- df_BC_pqqF$RA / 1.9

idx <- df_BCs$strain %in% c("WC358_cyo0_pBC163")
df_BC_cyoB <- df_BCs[idx,  ]
df_BC_cyoB$RA_norm <- df_BC_cyoB$RA / 1.8

df_BCs_cor <- rbind(df_BC_wt, df_BC_pqqF, df_BC_cyoB)

# root
idx <- df_BCs_cor$compartment %in% c("root") &
  df_BCs_cor$deriv %in% c("wt", "pqqF", "cyoB", "wt+pqqF", "wt+pqqF+cyoB") &
  T
BC_root <- df_BCs_cor[idx, ]

BC_root$deriv<- factor(BC_root$deriv, levels= c("wt", "pqqF", "cyoB", "wt+pqqF", 
                                                "wt+pqqF+cyoB"))
BC_root$strain<- factor(BC_root$strain, levels= c("WCS358_pBC92", "WC358_pqqF_pBC93",
                                                  "WC358_cyo0_pBC163"))

# matrix
idx <- df_BCs_cor$compartment %in% c("matrix") &
  df_BCs_cor$deriv %in% c("wt", "pqqF", "cyoB", "wt+pqqF",  "wt+pqqF+cyoB") &
  T
BC_matrix <- df_BCs_cor[idx, ]

BC_matrix$deriv<- factor(BC_matrix$deriv, levels= c("wt", "pqqF", "cyoB", "wt+pqqF", 
                                                    "wt+pqqF+cyoB"))

BC_matrix$strain<- factor(BC_matrix$strain, levels= c("WCS358_pBC92", "WC358_pqqF_pBC93",
                                                      "WC358_cyo0_pBC163"))

colors <- c(  "WCS358_pBC92" = "black",
              "WC358_pqqF_pBC93" = "#e07a5f",
              "WC358_cyo0_pBC163" = "#81b29a")


#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=11),
                    axis.text.x=element_text(colour="black", size=11, hjust = 1,vjust = 1,angle=45),
                    #use none to remove legend, can choose top, right, left
                    legend.position="top",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans", size=9))

# plot root
ggplot(BC_root, aes(x=deriv, y=RA_norm, fill = strain)) +
  geom_boxplot(outlier.shape = NA, size=.4, width=1, alpha=0.7) +
  geom_jitter(aes(color = strain),  position=position_jitterdodge(0.25), size=1, alpha=0.7) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  main_theme +
  labs(y="Normalized to spike")+
  ggtitle("Root_syncom") +
  ggsave(paste(figures.dir, "box_abs_BC_root_x_cond_v2.pdf", sep=""), width=4.4, height=4)

# plot matrix
ggplot(BC_matrix, aes(x=deriv, y=RA_norm, fill = strain)) +
  geom_boxplot(outlier.shape = NA, size=.4, width=1, alpha=0.7) +
  geom_jitter(aes(color = strain), position=position_jitterdodge(0.25), size=1, alpha=0.7) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  main_theme +
  labs(y="Normalized to spike")+
  ggtitle("Matrix_syncom") +
  ggsave(paste(figures.dir, "box_abs_BC_matrix_x_cond_v2.pdf", sep=""), width=4.4, height=4)

# stats #
BC_root %>%
  dunn_test(data =., RA ~ strain_deriv) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter(., p.adj<=0.05) %>%
  filter()->BC_root_stats

BC_matrix %>%
  wilcox_test(data =., RA ~ deriv) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter(., group1=="wt") %>% 
  filter()->BC_matrix_stats

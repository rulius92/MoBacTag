#############################################################################
############### Figure 4 ####################################################

### Scripts originally created by Ruben Garrido-Oter and modified by Ka-Wai Ma and Jana Ordon
# The scipt can be used to generate the PCA and CPcoA plots.

options(warn=-1)

# cleanup

rm(list=ls())

setwd("C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/")

library("scales")
library("grid")
library("vegan")
library("rstatix")
library("ggplot2")

# load plotting functions
source("C:/Users/Jana Ordon/Dropbox/R/MiSeq/plotting_functions.R")
source("C:/Users/Jana Ordon/Dropbox/R/MiSeq/cpcoa.func.R")

# directories
results.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/results/"
data.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/data/"
figures.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/subset_BC_vs_noBC/figures/"

# files
design.file <- paste(data.dir, "MoBacTag_V2_design.txt", sep="")
otu_table.file <- paste(results.dir, "merged_asv_table.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# combine columns in design file
design$deriv_compartment <- paste0(design$WCS417r_deriv,"_", design$compartment)


# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# only keep BC_wt vs no BC
idx <- design$WCS417r_deriv %in% c("wt", "wt_noBC") &
  TRUE

design <- design[idx, ]
otu_table <- otu_table[, idx]


# remove rows for BC of WCS358:pqqF and WCS358:cyoB mutants and 
# additional sequences, which were included in MiSeq run
idx <- rownames(otu_table) %in% c("Root418", "WC358_pqqF_pBC93", "WC358_cyo0_pBC163", 
                                  "Conta1_FromSpike", "pBCC069_Spike1_bc57","pBCC084_Spike2_bc72") 
otu_table <- otu_table[!idx,  ]

# spike normalized otu tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC0023"
pBCC0023 <- otu_table[idx, ]
pBCC0023 <- unlist(pBCC0023)
design <- cbind(design, pBCC0023)
otu_table <- otu_table[!idx,  ]
otu_table_norm <- sweep(otu_table, 2, pBCC0023, `/`)

################ Figure 4 a - c ######################################

# remove additional not needed counts for sequences
idx <- rownames(otu_table_norm) %in% c("WCS358_pBC92")
otu_table_norm_b <- otu_table_norm[!idx,  ]

idx <- rownames(otu_table_norm) %in% c("WCS358_pBC92", "WCS358")
otu_table_norm_c <- otu_table_norm[!idx,  ]


### compute beta diversity

colors <- c("wt_root" = "#3d405b", "wt_matrix" = "#3d405b", 
            "wt_noBC_root" ="#9EB3BD", "wt_noBC_matrix" ="#9EB3BD")
shapes <- c("a" = 17, "b" = 16, "c" = 15)


sqrt_transform <- T

########## PCoA Bray-Curtis: Figure 4 a

idx <- design$compartment %in% c("root", "matrix") &
  TRUE

d <- design[idx, ]
otu_table_norm <- otu_table_norm[, idx]


bray_curtis <- vegdist(t(otu_table_norm), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

# plot
ggplot(points, aes(x=x, y=y, color=deriv_compartment, shape = bio..Rep.)) +
  geom_point(alpha=.7, size=1.0) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values = shapes) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_abs_all.pdf", sep=""), width=2.4, height=3.2)


# stats 
# subset root
idx <- design$compartment %in% c("root") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)  

# subset rhizosphere
idx <- design$compartment %in% c("matrix") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

# stats
capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

########## PCoA Bray-Curtis: Figure 4 b

idx <- design$compartment %in% c("root", "matrix") &
  TRUE

d <- design[idx, ]
otu_table_norm_b <- otu_table_norm_b[, idx]


bray_curtis <- vegdist(t(otu_table_norm_b), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

# plot
ggplot(points, aes(x=x, y=y, color=deriv_compartment, shape = bio..Rep.)) +
  geom_point(alpha=.7, size=1.0) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values = shapes) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_abs_woBC.pdf", sep=""), width=2.4, height=3.2)


# stats 
# subset root
idx <- design$compartment %in% c("root") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)  

# subset rhizosphere
idx <- design$compartment %in% c("matrix") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

# stats
capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

########## PCoA Bray-Curtis: Figure 4 c

idx <- design$compartment %in% c("root", "matrix") &
  TRUE

d <- design[idx, ]
otu_table_norm_c <- otu_table_norm_c[, idx]


bray_curtis <- vegdist(t(otu_table_norm_c), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

# plot PCo 1 and 2

ggplot(points, aes(x=x, y=y, color=deriv_compartment, shape = bio..Rep.)) +
  geom_point(alpha=.7, size=1.0) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values = shapes) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_abs_woBC_woWCS358.pdf", sep=""), width=2.4, height=3.2)


# stats 
# subset root
idx <- design$compartment %in% c("root") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)  

# subset rhizosphere
idx <- design$compartment %in% c("matrix") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

# stats
capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

################ Figure 4 d - e ######################################

# remove all strain except for 
idx <- rownames(otu_table_norm) %in% c("Root101", "Root123D2", "Root1310", "Root131", "Root142",
                                       "Root265", "Root480", "Root61", "Root685", "Root68", "Root695",
                                       "Root77", "Root83", "Root935", "WCS358") 
otu_table_16S <- otu_table_norm[idx, ]


# dataframe with information plus RA
df_16S <- melt(as.matrix(otu_table_16S))
colnames(df_16S) <- c("strain", "SampleID", "RA")


#match sampleID with the corresponding WCS4174_deriv, compartment
df_16S$deriv <- design$WCS417r_deriv[match(df_16S$SampleID, design$SampleID)]
df_16S$compartment <- design$compartment[match(df_16S$SampleID, design$SampleID)]

##subset samples
#root
idx <- df_16S$compartment %in% c("root") &
  df_16S$deriv %in% c("wt", "wt_noBC") &
  T
syncom_root <- df_16S[idx, ]
syncom_root$deriv<- factor(syncom_root$deriv, levels= c("wt", "wt_noBC"))

#matrix
idx <- df_16S$compartment %in% c("matrix") &
  df_16S$deriv %in% c("wt", "wt_noBC") &
  T
syncom_matrix <- df_16S[idx, ]
syncom_matrix$deriv<- factor(syncom_matrix$deriv, levels= c("wt", "wt_noBC"))


colors <- c(  "wt_noBC" ="#9EB3BD",
              "wt" = "#3d405b")


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

##plot root 16S
ggplot(syncom_root, aes(x=strain, y=RA, fill = deriv)) +
  geom_boxplot(alpha=0.7, outlier.shape= NA, size=.4, width=.85) +
  geom_jitter(position=position_jitterdodge(0.3), size=.4, alpha=0.7) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  main_theme +
  ylim (0, 1.8) +
  labs(y="Normalized to spike")+
  ggtitle("Root_all_16S") +
  ggsave(paste(figures.dir, "box_abs_16S_root.pdf", sep=""), width=7, height=4)

## plot matrix Syncom
ggplot(syncom_matrix, aes(x=strain, y=RA, fill = deriv)) +
  geom_boxplot(alpha=.7, outlier.shape = NA, size=.4, width=.85) +
  geom_jitter(position=position_jitterdodge(0.3), size=.4, alpha=0.7) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  main_theme +
  ylim (0, 25) +
  labs(y="Normalized to spike")+
  ggtitle("Matrix_16S") +
  ggsave(paste(figures.dir, "box_abs_16S_matrix.pdf", sep=""), width=7, height=4)

syncom_root %>%
  group_by(strain) %>%
  dunn_test(data =., RA ~ deriv) %>%
  filter()->stats_root

syncom_matrix %>%
  group_by(strain) %>%
  dunn_test(data =., RA ~ deriv) %>%
  filter()->stats_matrix

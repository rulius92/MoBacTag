############## Figure 5 d - e ###############################################

# Script originally created by Ruben Garrido-Oter and modified by Ka-Wai Ma and Jana Ordon
# The script can be used to generate the PCA and CPcoA plots.

rm(list=ls())

# set working directory
setwd("C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/")

# load plotting functions
source("C:/Users/Jana Ordon/Dropbox/R/MiSeq/plotting_functions.R")
source("C:/Users/Jana Ordon/Dropbox/R/MiSeq/cpcoa.func.R")

# directories
results.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/results/"
data.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/data/"
figures.dir <- "C:/Users/Jana Ordon/Dropbox/R/MiSeq/MoBacTag_V2/figures/"

# files
design.file <- paste(data.dir, "MoBacTag_V2_design.txt", sep="")
otu_table.file <- paste(results.dir, "merged_asv_table.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# remove BCs, remove WCS358 --> relative differences in community
idx <- rownames(otu_table) %in% c("pBCC069_Spike1_bc57", "Conta1_FromSpike", 
                                  "pBCC084_Spike2_bc72", "WCS358_pBC92", "WC358_pqqF_pBC93", 
                                  "WC358_cyo0_pBC163", "WCS358")
depleted_strains <- otu_table[idx, ]
otu_table <- otu_table[!idx,  ]

# spike normalized otu tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC0023"
pBCC0023 <- otu_table[idx, ]
pBCC0023 <- unlist(pBCC0023)
design <- cbind(design, pBCC0023)
otu_table <- otu_table[!idx,  ]
otu_table_norm <- sweep(otu_table, 2, pBCC0023, `/`)


### CPCoA Bray-Curtis of samples conditioned by technical factors
colors <- c(  "wt" = "black",
              "pqqF" = "#e07a5f",
              "cyoB" = "#81b29a",
              "wt+pqqF" = "#9A381D",
              "wt+pqqF+cyoB" = "#ECB45B")


sqrt_transform <- T

# CPcoA for single root samples
idx <- design$compartment %in% c("root") &  
  design$WCS417r_deriv %in% c("wt", "cyoB", "pqqF") &
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")
capscale.gen <- capscale(bray_curtis ~ WCS417r_deriv + Condition(bio..Rep. * tech.rep.), data=d, add=F, sqrt.dist=sqrt_transform)
perm_anova.gen <- anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[, 1:2]
colnames(points_cpcoa) <- c("x", "y")
points_cpcoa <- cbind(points_cpcoa, d[match(rownames(points_cpcoa), d$SampleID), ])

# plot
ggplot(points_cpcoa, aes(x=x, y=y, color=WCS417r_deriv, shape=bio..Rep.)) +
  geom_point(alpha=.8, size=2) +
  scale_colour_manual(values=colors) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("CPCoA of root samples", paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="bottom")+
  ggsave(paste(figures.dir, "CPCoA_root_single_16S_woWCS358_abs.pdf", sep=""),  height=4, width=2.7)

### CPcoA for wt vs pqqF vs wt+pqqF
idx <- design$compartment %in% c("root") & 
  design$WCS417r_deriv %in% c("wt", "pqqF", "wt+pqqF", "wt+pqqF+cyoB") &
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")
capscale.gen <- capscale(bray_curtis ~ WCS417r_deriv + Condition(bio..Rep. * tech.rep.), data=d, add=F, sqrt.dist=sqrt_transform)
perm_anova.gen <- anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[, 1:2]
colnames(points_cpcoa) <- c("x", "y")
points_cpcoa <- cbind(points_cpcoa, d[match(rownames(points_cpcoa), d$SampleID), ])

##plot
ggplot(points_cpcoa, aes(x=x, y=y, color=WCS417r_deriv, shape=bio..Rep.)) +
  geom_point(alpha=.8, size=2.0) +
  scale_colour_manual(values=colors) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("CPCoA of root samples", paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="bottom")+
  ggsave(paste(figures.dir, "CPCoA_root_16S_abs_noWCS358_wt_pqqF_competition.pdf", sep=""),  height=4, width=4.2)

# stats 

# subset root
# change for each comparison
idx <- design$compartment %in% c("root") &  
  design$WCS417r_deriv %in% c("wt+pqqF", "wt+pqqF+cyoB") &
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

# stats
capscale.gen <- capscale((bray_curtis ~ WCS417r_deriv), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)  

############## Figure 5 f - g ###############################################

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
  df_16S$deriv %in% c("wt", "pqqF", "cyoB") &
  T
syncom_root_single <- df_16S[idx, ]

syncom_root_single$deriv<- factor(syncom_root_single$deriv, levels= c("wt", "pqqF", "cyoB"))

colors <- c(  "wt" = "#3d405b",
              "pqqF" = "#e07a5f",
              "cyoB" = "#81b29a")


#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    axis.text.x=element_text(colour="black", size=14, hjust = 1,vjust = 1,angle=45),
                    #use none to remove legend, can choose top, right, left
                    legend.position="top",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans", size=15))

########
#subset single data sets in high and low abundant strains

##root
idx <- syncom_root_single$strain %in% c("Root1310", "Root480", "Root61", "Root68", "Root935") &
  T
syncom_root_single_high <- syncom_root_single[idx, ]
syncom_root_single_low <- syncom_root_single[!idx, ]

ggplot(syncom_root_single_high, aes(x=strain, y=RA, fill = deriv)) +
  geom_boxplot(aes(fill = deriv), alpha=0.7, outlier.shape = NA, size=.4, width=.85) +
  geom_jitter(aes (color = deriv), shape=16, position=position_jitterdodge(0.25), size=.4, alpha=0.7) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  main_theme +
  labs(y="Normalized to spike")+
  ggtitle("Root_16S_single_high")+
  ggsave(paste(figures.dir, "box_abs_16S_single_high_root.pdf", sep=""), width=2.6, height=3.7)

ggplot(syncom_root_single_low, aes(x=strain, y=RA, fill = deriv)) +
  geom_boxplot(aes(fill = deriv), alpha=0.7, outlier.shape= NA, size=.4, width=.85) +
  geom_jitter(aes(color = deriv), shape= 16, position=position_jitterdodge(0.25), size=.4, alpha=0.7) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  main_theme +
  labs(y="Normalized to spike")+
  ggtitle("Root_16S_single_low")+
  ggsave(paste(figures.dir, "box_abs_16S_single_low_root.pdf", sep=""), width=4.5, height=3.7)

# stats #
# change according to comparisons
syncom_root_single %>%
  group_by(strain) %>%
  wilcox_test(data =., RA ~ deriv) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter(., group1=="wt") %>% 
  filter(., group2=="cyoB") %>%
  filter()->single_stats

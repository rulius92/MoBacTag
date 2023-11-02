###################################################################################
############### Extended Data Figure 3 A B ########################################



# originally by Julien Thouin
# thouin@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())


# libraries
library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("readxl")
library("multcompView")
library("tidyverse")
library("gridExtra")
library("lattice")
library("cowplot")
library('splitstackshape')
library('ggpubr')
library("ggsci")
library("reshape2")
library("data.table")
library("dplyr")
library("tidyr")
library("stringr")
library("RColorBrewer")
library("reshape")
library("wesanderson")

# folder localisation
path <- "/Users/thouin/Nextcloud/Post-Doc-Koln/Experiments/DNA_Bar_Code/StandardCurve_Primers/"


# process independently a given experiment (run)
run <- "Standard_Curve_Primers"

# directories

results.dir <- paste(path, "results/", sep="")
data.dir <- paste(path, "data/", sep="")
figures.dir <- paste(path, "figures/", sep="")


# files
design.file <- paste(data.dir, "mapping_Test-Primers.xlsx", sep="")
otu_table.file <- paste(data.dir, "Merged_ASV_table.xlsx", sep="")

# load data
design <- read_excel(design.file, sheet = "mapping")
otu_table <- read_excel(otu_table.file)

#RowName otu_table
names <-otu_table$Strain


# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]
rownames(otu_table) <- names


#re-arrangment of data
otu_table[c("pBCC084_Spike2_bc72"),] <- otu_table[c("pBCC084_Spike1_bc72"),] + otu_table[c("pBCC084_Spike2_bc72"),]
rownames(otu_table) <- names

otu_table[c("Root9"),] <- otu_table[c("cont_1_Root9"),] + otu_table[c("Root9"),]
rownames(otu_table) <- names

otu_table[c("F212"),] <- otu_table[c("Cont5_F212"),] + otu_table[c("F212"),]
rownames(otu_table) <- names

names -> otu_table$Strain
rownames(otu_table) <- names

#Remove contamination line

conta = c("cont_1_Root9","pBCC084_Spike1_bc72", "cont_2_pBCC069_bc57","cont_3_pBCC069_bc57","cont_4pBCC084_Spike1_bc72","cont_5_pBCC084_Spike1_bc72","Cont5_F212","Conta1_FromSpike")

idx <- rownames(otu_table) %in%  conta
otu_table <- otu_table[!idx,]

names <- otu_table$Strain

#Remove P046 because it contain plant reads while it's matrix sample, an error occured
idx <- design$SampleID %in% c("P046")
design <- design[!idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]
rownames(otu_table) <- names

# normalize otu tables
#Because of Excel import, I need to transform otu_table to numeric value
otu_table <- type.convert(otu_table, na.strings = "NA", as.is = FALSE, dec = ".", numerals = "allow.loss")

rownames(otu_table) <- names


#Use of Spike as reference
Spike <- colSums(otu_table[c("pBCC069_Spike1_bc57"),])
Spike <- t(Spike)
#library('splitstackshape') / Repeat row from Spike
Spike <- expandRows(Spike, count=nrow(otu_table), count.is.col=FALSE)
otu_table <- otu_table/Spike

#Merge otu_table and design
otu_table <- t(as.matrix(otu_table))

design <- cbind(design,otu_table)

df <- melt(design, id.vars = c('SampleID','compartment',"primers","concentration"),
                measure.vars = c("pBCC084_Spike2_bc72"))

colnames(df) <- c('SampleID','compartment',"primers","concentration","sequence","value")

#Log transform concentration
#df$concentration <- log10(df$concentration)

#Remove concentration > 1.5
#idx <- df$concentration < 1.5
#df <- df[idx, ]


##Theme
theme_RTN2 <- theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid=element_blank(),
        legend.key=element_blank(),
        legend.text.align=0,
        strip.text=element_text(face="bold"),
        axis.line=element_line(),

        plot.title = element_text(color="black", size=24),

        legend.title = element_text(colour="black", size=20),
        legend.text = element_text(colour="black", size=18),
        legend.key.size = unit(2, 'cm'),

        axis.line.x=element_line(),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.x = element_text(colour="black", size=17),
        #axis.title.x = element_blank(),
        
        axis.text.y = element_text(colour="black", size=18),
        axis.title.y = element_text(colour="black", size=20, vjust =1),
        axis.line.y=element_line(),
        
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.background=element_rect(fill="transparent", colour=NA),
        strip.placement="outside")

  lettersize = 12
  yaxis = "Sequence abundance (Spike normalized)"
  xaxis = "Concentration (ng)"

al=0.4
labelsize=9



#fungalITS - ITS
idx <- df$primers %in% c("fungalITS")
df2 <- df[idx, ]

#Plot nuage de points
b <-  ggplot(df2, aes(x=concentration, y=value, color=compartment)) + 
  geom_point() + 
  geom_smooth(method=lm, formula = y~x)+
  scale_color_brewer(palette="Dark2") + 
  stat_regline_equation(size= labelsize,
          aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  labs(title = "ITS - ITS1/ITS2", y=yaxis, x=xaxis)+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_RTN2 
    
ggsave(paste(figures.dir, "ITS.pdf", sep=""), b, width=10, height=10)


#pITS
idx <- df$primers %in% c("plantITS")
df2 <- df[idx, ]




#Barcode speficic primers
idx <- df$primers %in% c("BCspecific")
df2 <- df[idx, ]

#Plot nuage de points
d <-  ggplot(df2, aes(x=concentration, y=value, color=compartment)) + 
  geom_point() + 
  geom_smooth(method=lm, formula = y~x)+
  scale_color_brewer(palette="Dark2") + 
  stat_regline_equation(size= labelsize,
          aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  labs(title = "Barcode specific", y=yaxis, x=xaxis)+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_RTN2 
    
ggsave(paste(figures.dir, "BC.pdf", sep=""), d, width=10, height=10)


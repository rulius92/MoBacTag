#############################################################################
############### Figure 2 D ##################################################



# originally by Julien Thouin
# thouin@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# folder localisation
path <- "/Volumes/groups/dep_psl/grp_psl/rooters/Tn7_Barcodes/Barcode_Proof-of-concept/MiSeq_analysis_Jana/"


# load plotting functions
library("ggplot2")
library("vegan")
library("tidyverse")
library("reshape2")
library("ggpubr")
library('splitstackshape')
library("readxl")


# process independently a given experiment (run)
run <- "FinalFigure"

# directories
results.dir <- paste(path, "results/", sep="")
data.dir <- paste(path, "data/", sep="")
figures.dir <- paste(path, "figures/", sep="")


# files
design.file <- paste(path, "results/16SvsBC_WCS358.xlsx", sep="")
###########
#Data was manually extracted form ASV_table of the gnotopot experiment



# load data
design <- read_excel(design.file)


##Theme
theme_RTN2 <- theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid=element_blank(),
        legend.key=element_blank(),
        legend.text.align=0,
        strip.text=element_text(face="bold"),
        axis.line=element_line(),

        plot.title = element_text(color="black", size=24),

        legend.title = element_text(colour="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.key.size = unit(2, 'cm'),

        axis.line.x=element_line(),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.x = element_text(colour="black", size=17),
        #axis.title.x = element_blank(),
        
        axis.text.y = element_text(colour="black", size=16),
        axis.title.y = element_text(colour="black", size=18),
        axis.line.y=element_line(),
        
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.background=element_rect(fill="transparent", colour=NA),
        strip.placement="outside")

  
#colors <- brewer.pal(n = 11, name = "Paired")
labelsize=9

# WT #010101
# PQQf #E07A5F
# CYOb #82B29A




#Plot nuage de points
p1 <-  ggplot(design, aes(x=Endo16S, y=BC, colour=genotype)) + 
  geom_point() + 
  scale_colour_manual(values=c("#82B29A","#E07A5F","#010101"))+
 stat_smooth(method = "lm", intercept = 0, se=FALSE,formula = y~x+0,) + 
 # scale_color_brewer(palette="Dark2") + 
  stat_regline_equation(formula = y~x+0,) +
  stat_cor(label.x = 2.5)+

  
 # expand_limits(x = 0, y = 0)+
  labs(title = "16S vs Barcode sequences", y="Barcode sequences abundance (Spike normalized)", x="16S sequences abundance\n(Spike normalized)")+
  theme_RTN2 
    
ggsave(paste(path, run, "16SvsBC_WCS358-Origin.pdf", sep=""),width=10, height=10, p1)



#############################################################################
############### Figure 2 E ##################################################


# originally by Julien Thouin
# thouin@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# folder localisation
path <- "/Users/thouin/Desktop/Julien Relecture Figures/"


# load plotting functions
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

# process independently a given experiment (run)

run <- "FinalFigure"

# directories

results.dir <- paste(path, "results/", sep="")
data.dir <- paste(path, "data/", sep="")
figures.dir <- paste(path, "figures/", sep="")


# files
design.file <- paste(path, "MultipleBarcode effect on amplification_R.xlsx", sep="")
###########
#Data was manually extracted form ASV_table of gnotopot experiment with R13D strains


# load data
design <- read_excel(design.file)


##Theme
theme_RTN2 <- theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid=element_blank(),
        legend.key=element_blank(),
        legend.text.align=0,
        strip.text=element_text(face="bold"),
        axis.line=element_line(),

        plot.title = element_text(color="black", size=24),

        legend.title = element_text(colour="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.key.size = unit(2, 'cm'),

        axis.line.x=element_line(),
        axis.text.x = element_text(colour="black", size=18),
        axis.title.x = element_text(colour="black", size=17),
        #axis.title.x = element_blank(),
        
        axis.text.y = element_text(colour="black", size=16),
        axis.title.y = element_text(colour="black", size=18),
        axis.line.y=element_line(),
        
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.background=element_rect(fill="transparent", colour=NA),
        strip.placement="outside")

  
colors <- brewer.pal(n = 11, name = "Paired")

labelsize=9

#Plot nuage de points
p1 <-  ggplot(design, aes(x=R13D_16S, y=Tag_reads, color=tag)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, formula = y~x+0)+
  scale_color_brewer(palette="Dark2") + 
  stat_regline_equation(method=lm, se=FALSE, fullrange = TRUE ,formula = y~x+0,label.x = 5,size= labelsize,
          aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
 # expand_limits(x = 0, y = 0)+
  labs(title = "16S vs Barcode sequences", y="Barcode sequences abundance (Spike normalized)", x="16S sequences abundance\n(Spike normalized)")+
  theme_RTN2 
    
ggsave(paste(path, run, "16SvsBC.pdf", sep=""),width=10, height=10, p1)






#############################################################################
############### Figure 2 F ##################################################

# originally by Julien Thouin
# thouin@mpipz.mpg.de

options(warn=-1)

# cleanup
rm(list=ls())

# folder localisation
path <- "/Volumes/groups/dep_psl/grp_psl/rooters/Tn7_Barcodes/Barcode_Proof-of-concept/old/"


# load plotting functions
#source(paste(path, "scripts/plotting_functions.R", sep=""))
#source(paste(path, "scripts/plotting_parameters.R", sep=""))
#source(paste(path, "scripts/cpcoa.func.R", sep=""))

# load plotting functions
library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("readxl")
library("reshape2")
library("multcompView")
library("data.table")
library("tidyverse")
library("dplyr")
library("tidyr")
library("cowplot")
library("stringr")
library("gridExtra")
library("lattice")
library('splitstackshape')
library('ggpubr')
library("RColorBrewer")
library("reshape")
library("wesanderson")
library("reshape2")


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
        legend.key.size = unit(1, 'cm'),

        axis.line.x=element_line(),
        axis.text.x = element_text(colour="black", size=15, angle=45,hjust = 0.5, vjust = 0.5),
        #axis.title.x = element_text(colour="black", size=17),
        axis.title.x = element_blank(),
        
        axis.text.y = element_text(colour="black", size=15),
        axis.title.y = element_text(colour="black", size=17, vjust =1),
        axis.line.y=element_line(),
        
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.background=element_rect(fill="transparent", colour=NA),
        strip.placement="outside")

# process independently a given experiment (run)
run <- "Jan2023_4barcodes_gDNA_Corrected_Absolute_2SpikeSums"

# directories
results.dir <- paste(path, "results/", sep="")
data.dir <- paste(path, "data/", sep="")
figuresPCOA.dir <- paste(path, "figures/PCOA/", sep="")
figuresPlot.dir <- paste(path, "figures/BoxPlot/", sep="")
figures.dir <- paste(path, "figures/", sep="")

# files
design.file <- paste(data.dir, "MiSeq_Barcode_proof-of-concept.xlsx", sep="")
otu_table.file <- paste(data.dir, "BCProofOfConcept_merged_asv_table.xlsx", sep="")

# load data
design <- read_excel(design.file, sheet = "mapping")
otu_table <- read_excel(otu_table.file)
names <- otu_table$`#Strain`



# re-order data matrices
idx <- design$`#SampleID...1` %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$`#SampleID...1`, colnames(otu_table))
otu_table <- otu_table[, idx]


# normalize otu tables
#Because of Excel import, I need to transform otu_table to numeric value
otu_table <- type.convert(otu_table, na.strings = "NA", as.is = FALSE, dec = ".", numerals = "allow.loss")

#Rename line for each strain
rownames(otu_table) <- names

#Use of Spike as reference
#Each spike represent 0,0005 of DNA per PCR reaction
idx <- rownames(otu_table) %in% c("pBCC069_Spike1_bc57","pBCC084_Spike1_bc72")  #"pBCC069_Spike1_bc57" didn't work for all of them
Spike <- otu_table[idx,]
#Mean of the two spike
Spike <- colSums(Spike)
Spike <- t(Spike)
Spike <- data.frame(Spike)
#library('splitstackshape') / Repeat row from Spike
#number of colum otu_table
n <- nrow(otu_table)
Spike <- expandRows(Spike, count=n, count.is.col=FALSE)

otu_table2 <- otu_table/Spike


#Remove spike
idx <- rownames(otu_table2) %in% c("pBCC069_Spike1_bc57","pBCC084_Spike1_bc72")
otu_table2<- otu_table2[!idx,]

#Only sequence from gDNA barcoded
idx <- rownames(otu_table2) %in% c("R13D","R13Dmut685","R13mut420","R13Dmut695","R13Dmut140")
otu_table2<- otu_table2[idx,]

##subset samples of interest 
samples<- c("gDNA")
idx <- design$compartment %in% samples
design <- design[idx,]
otu_table2 <- otu_table2[, idx]

#Normalize data in %
#otu_table2 <- apply(otu_table2, 2, function(x) x/sum(x)*100)



#Bind otu_table and design
otu_table2 <- t(otu_table2)
otu_table2 <- data.frame(otu_table2)
design <- cbind(design, otu_table2)



#Mean of the 3 replicates
idx <- colnames(design) %in% c("SynCom","R13D","R13Dmut685","R13mut420","R13Dmut695","R13Dmut140")
design <- design[,idx]

#Correction BC vs 16S
design$R13Dmut685 <- design$R13Dmut685 / 1.6 #tag a
design$R13mut420 <- design$R13mut420 / 1.4  #tag b
design$R13Dmut695 <- design$R13Dmut695 / 1.4 #tag c
design$R13Dmut140 <- design$R13Dmut140 / 0.87 #tag d



design1<- aggregate(design[, 2:6], list(design$SynCom), mean)

design2<- aggregate(design[, 2:6], list(design$SynCom),  function(x) sd = sd(x))




#Convert Data for ggplot2
#library("reshape")
df <- melt(design1,id.vars = c("Group.1"),
                measure.vars = c("R13D","R13Dmut685","R13mut420","R13Dmut695","R13Dmut140"))


df2 <- melt(design2,id.vars = c("Group.1"),
                measure.vars = c("R13D","R13Dmut685","R13mut420","R13Dmut695","R13Dmut140"))

df3 <- melt(design,id.vars = c("SynCom"),
                measure.vars = c("R13D","R13Dmut685","R13mut420","R13Dmut695","R13Dmut140"))

colnames(df) <- c("barcode","sequence","absolue")
colnames(df2) <- c("barcode","sequence","sd")

colnames(df3) <- c("barcode","sequence","absolue")



df$sd <- df2$sd


#ReOrderData for plotting and Stat
#Plotting
df$barcode <- factor(df$barcode, levels = c("BCa","BCb","BCc","BCd","BCa, BCb","BCa, BCc","BCa, BCd","BCa, BCb, BCc","BCa, BCb, BCd","BCa, BCb, BCc, BCd"))


#Only more than one tag

idx <- df$barcode %in% c("BCa, BCb","BCa, BCc","BCa, BCd","BCa, BCb, BCc","BCa, BCb, BCd","BCa, BCb, BCc, BCd")
df <- df[idx,]

#Mean of the 3 replicates
idx <- df3$barcode %in% c("BCa, BCb","BCa, BCc","BCa, BCd","BCa, BCb, BCc","BCa, BCb, BCd","BCa, BCb, BCc, BCd")
df3 <- df3[idx,]

#############
####PLOT#####
#############
#Boxplot#
#############

#Settings

al=0.5

h=9
w=30


#Background
bblack=data.frame(x1=c(0.5,2.5,4.5), ymin=-Inf, ymax=Inf)
bblack$x2 <- bblack$x1+1

bgrey=data.frame(x1=bblack$x1+1,x2=bblack$x2+1 ,ymin=-Inf, ymax=Inf)

bgrey <- bgrey[-c(9),] 



p1<- ggplot(df, aes(x=barcode, y=absolue, fill=sequence))+
    geom_bar(stat = "identity", position=position_dodge())+
    geom_errorbar(aes(ymin=absolue-sd, ymax=absolue+sd), width=.2,
                 position=position_dodge(.9)) +
    geom_jitter(data=df3,inherit.aes=FALSE, aes(x=barcode, y=absolue, fill=sequence),
        position=position_jitterdodge(
                      jitter.width = 0.35,
                      jitter.height = 0,
                      dodge.width = 1,
                      seed = NA),
                        size=2)+
    scale_fill_manual(
        values=c("#FF0000","#00A08A","#F2AD00","#F98400","#5BBCD6"),
        name  ="Sequences",
        breaks=c("R13D","R13Dmut685","R13mut420","R13Dmut695","R13Dmut140"),
        labels=c("16S","Tag a","Tag b","Tag c","Tag d"))+

    geom_rect(data=bblack, inherit.aes = FALSE, aes(xmin=x1, xmax=x2, ymin=ymin, ymax=ymax), fill="black", alpha=0.1) +
    geom_rect(data=bgrey, inherit.aes = FALSE, aes(xmin=x1, xmax=x2, ymin=ymin, ymax=ymax), fill="grey", alpha=0.1) +

    geom_vline(xintercept = c(1.5,2.5,4.5), 
        linetype="dotted",color = "black", size=0.5) +
    geom_vline(xintercept = c(3.5,5.5), 
       color = "black", size=0.8) +
    
    scale_x_discrete(labels = c("BCa"="Tag a","BCb"="Tag b","BCc"="Tag c","BCd"="Tag d","BCa, BCb"="Tag a + b","BCa, BCc"="Tag a + c","BCa, BCd"="Tag a + d","BCa, BCb, BCc"="Tag a + b + c","BCa, BCb, BCd"="Tag a + b + d","BCa, BCb, BCc, BCd"="Tag a + b + c + d"))+
    labs(title = "Effect of the number of barcodes on the amplification rate", y="Absolute quantification\n(Normalized to spike)", x="")+
theme_RTN2 

ggsave(paste(figures.dir, run, "_BoxPlot.pdf", sep=""),width=h, height=h, p1)



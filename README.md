# Pipeline to reproduce the Figures from the MoBacTag publication (https://doi.org/10.1101/2023.04.20.537712 )

---------------------------

## Scripts of Ordon, Thouin, *et. al*, Simultaneous tracking of near-isogenic bacterial strains in synthetic Arabidopsis microbiota by chromosomally-integrated barcodes.

These scripts are made available to facilitate the reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper.

---------------------------

### Accession numbers and desing files:

Raw *16S* rRNA amplicon reads were deposited in the European Nucleotide Archive (ENA) under the accession number PRJEB61076.

The design files for the gnotobiotic SynCom experiment and the validation of the oligonucleotides:
- MoBacTag_V2_design.txt
- Validation-Primers_design.txt

---------------------------

### Scripts used for processing data, creating the figures and performing the statistical analysis reported in the manuscript:

Raw data corresponding to the SynCom experiment or the validation of the oligonucleotides were processed as described in the manuscript and using the script syncom_ASV-table.sh

The scripts to generate the plots from the figures are listed below:

Figure1: 
- Figure1C.R

Figure2: 
- Figure2A-D.R
- Figure2E-F.R

Figure4:
- Figure4.R

Figure 5: 
- Figure5A-C.R
- Figure5D-E.R

Figure6:
- Figure6.R

---------------------------

### Supplementary data:

Manually curated data are present as follow: 
- frequency_primer_per_position.txt
- WCS358_16SvsBC.txt
- MultipleBarcode.txt

---------------------------

## Authors

Jana Ordon (ordon@mpipz.mpg.de)

Julien Thouin (thouin@mpipz.mpg.de)

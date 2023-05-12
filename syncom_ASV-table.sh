### This pipeline is used for exact 
### sequence variants analysis of   
### microbiome data                 

### Originally by Pengfan Zhang
### Email: pzhang@mpipz.mpg.de

### Modified by Jana Ordon and Julien Thouin
### ordon@mpipz.mpg.de; thouin@mpipz.mpg.de

echo "All the sequencing data from a single \
lane should be stored in an exclusive \
directory and the file names must be \
changed to forward.fastq.gz, reverse.fastq.gz \
and barcode.fastq.gz."

# get path to scripts
scripts_dir=`dirname $0`

#########################
## activate QIIME2, etc
#########################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/etc/profile.d/conda.sh" ]; then
        . "/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/etc/profile.d/conda.sh"
    else
        export PATH="/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/biodata/dep_psl/grp_psl/Jana/Results/MiSeq/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate qiime2-2021.2


#################################
## load the config informations
#################################
# The absolute path of the mapping file
Mapping_file="/netscratch/dep_psl/grp_psl/JulienT/Rbec/V2-16S_ITS_TestPrimers/mapping_file.txt"

# The absolute path of the dereplicated reference sequences
Ref="/netscratch/dep_psl/grp_psl/JulienT/Rbec/V2-16S_ITS_TestPrimers/Ref_Seq_FG.fasta"

# The directory containing the raw data
DataPath="/netscratch/dep_psl/grp_psl/JulienT/Rbec/V2-16S_ITS_TestPrimers/data"

# Library ID
# This parameter only supports one input each time
l_list="TestPrimers_V2_16S-ITS"

# the maximum length and maximum percent of mismatches of overlapped regions when merging reads by Flash2
max_overlap=350
max_m=0.25

# number of subsampled reads for error probability calculation
depth=2000

# number of threads for parallel steps
threads=2

# the sequencing quality encoding system 
# The sequencing quality from our MiSeq platform in the basement is encoded by phred 33. If the data is from a HiSeq platform, this might be phred 64. To be sure of that, you need to check the quality score.
phred=33




CurrDir=`pwd`

## Import data
mkdir $CurrDir/Results
qiime tools import \
	--type EMPPairedEndSequences \
	--input-path $DataPath \
	--output-path $CurrDir/Results/${l_list}.qza


## Demultiplex
qiime demux emp-paired \
	--m-barcodes-file $Mapping_file \
	--m-barcodes-column BarcodeSequence \
	--p-rev-comp-mapping-barcodes \
	--i-seqs $CurrDir/Results/${l_list}.qza \
	--o-per-sample-sequences $CurrDir/Results/${l_list}_syncom_demux.qza \
	--o-error-correction-details $CurrDir/Results/${l_list}_syncom_demux-details.qza \
#	--p-no-golay-error-correction

#No golay correction in case of double barcoding

qiime demux summarize \
	--i-data $CurrDir/Results/${l_list}_syncom_demux.qza \
	--o-visualization $CurrDir/Results/${l_list}_syncom_demux.qzv

qiime tools export \
	--input-path $CurrDir/Results/${l_list}_syncom_demux.qza \
	--output-path $CurrDir/Results/${l_list}_syncom/

## Merge Reads
for F1 in $CurrDir/Results/${l_list}_syncom/*L001_R1_001.fastq.gz
do
	R1=${F1/R1_001.fastq.gz/R2_001.fastq.gz}
	random_num=`basename $F1 _L001_R1_001.fastq.gz| rev |awk -F '_' '{print $1}' | rev`
	samplename=`basename $F1 _${random_num}_L001_R1_001.fastq.gz`
	flash2 $F1 $R1 -M $max_overlap -o $samplename -x $max_m -d $CurrDir/Results/${l_list}_syncom/ -p $phred
done
rm $CurrDir/Results/${l_list}_syncom/*.hist
rm $CurrDir/Results/${l_list}_syncom/*.histogram
rm $CurrDir/Results/${l_list}_syncom/*notCombined*

## Filter reads with Ns
for seqs in $CurrDir/Results/${l_list}_syncom/*.extendedFrags.fastq
do
	name=`basename $seqs .extendedFrags.fastq`
	usearch -fastq_filter $seqs -fastqout $CurrDir/Results/${l_list}_syncom/$name.extendedFrags_filtered.fastq -fastq_maxns 0
done

awk 'BEGIN{seqs=""}{if(/^>/){if(seqs!=""){print seqs;seqs=""}; print $0} else{seqs=seqs""$0}}END{print seqs}' $Ref > $CurrDir/Results/${l_list}_new_ref.fasta
#Generate ref-based error correction bash scripts for each sample
mkdir $CurrDir/Results/Correction
for sc in $CurrDir/Results/${l_list}_syncom/*extendedFrags_filtered.fastq
do
	name=`basename $sc .extendedFrags_filtered.fastq`
	mkdir -p $CurrDir/Results/Correction/$name
	echo "/biodata/dep_psl/grp_rgo/pzhang/software/anaconda2/envs/qiime2-2019.4/bin/Rscript /netscratch/dep_psl/grp_rgo/pzhang/mock/scripts/Rbec.R $sc $CurrDir/Results/${l_list}_new_ref.fasta $threads $depth $phred $name" > $CurrDir/Results/Correction/bsub_correction_$name.sh
done

for i in $CurrDir/Results/Correction/*
do
	if [ -d $i ];then
		echo "$i/strain_table.txt" >> $CurrDir/Results/Correction/table.list
	fi	
done

grep '^>' $CurrDir/Results/${l_list}_new_ref.fasta|sed 's/>//' > $CurrDir/Results/Correction/strain.list

echo "Now please submit the reference-based-error-correction scripts for each sample in the folder \
$CurrDir/Results/Correction to a node seperately \
, so all the samples could run parallelly!
In order to do so, firstly you have to log into hpcws001 or other hpcwsXX nodes with the command 'ssh hpcws001';
then enter the directory with 'cd $CurrDir/Results/Correction' and submit the scripts simultaneously with the following command:
	'for i in $CurrDir/Results/Correction/bsub*.sh;do name=\`basename \$i .sh\`; bsub -n [number of thread] -M [memory] -o \$name.log -q [queue] \$i;done'.
Finally, use the following command to merge the count table form each sample:
	'perl /netscratch/dep_psl/grp_rgo/pzhang/tools/16S/syncom_pipeline/merge_table.pl $CurrDir/Results/Correction/strain.list $CurrDir/Results/Correction/table.list > merged_asv_table.txt'"



#!/bin/bash
###########
#simple script for Adapter triming (CutAdapt) and fastqc report of it
#within the folder  (DATASET_FOLDER_NAME) it takes all files with
# specified suffix (SUFFIX) 


DATASET_FOLDER_NAME=/storage/brno2/home/tarakaramji/Data_test
OUTPUT_FOLDER_FASTQC=/storage/brno2/home/tarakaramji/results/fastqc
SUFFIX="fastq.gz"  # raw file format
cd $SCRATCH # home folder of the raw files
mkdir -p $OUTPUT_FOLDER_FASTQC  # making directory to store the fastqc data from the raw reads
module load fastQC-0.11.5  #load module for FASTQC
module load python27-modules-gcc  # load module for Cutadapt
cp $DATASET_FOLDER_NAME/*$SUFFIX .   # copy files in the working directory

###CutAdapt variables##
##ADAPTER TRIMMING AND FILTERING
FILE_FORMAT=fastq
OUTPUT_FOLDER_CUTADAPT=/storage/brno2/home/tarakaramji/results/cutadapt/     # file to store all the cutadapt files others
OUTPUT_FOLDER_CUTADAPTgz=/storage/brno2/home/tarakaramji/results/cutadapt/trimmed/  # file to store the trimmed files separately for FASTQC processing
mkdir -p $OUTPUT_FOLDER_CUTADAPT  # making directory for the output of the cutadapt 
mkdir -p $OUTPUT_FOLDER_CUTADAPTgz # making directory for the trimmed files of the cutadapt
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

ERROR_RATE=0.10 #allowed error rate or adapters
MIN_OVERLAP=5 #minimal overlap of adapters; because we afterwards filter by number of mismatches AND WE ARE NOT
		# primary looking for isoforms we can go and trim by overlap of a single nucleotide
MAX=1 # Maximal number of adapters to remove
DISC_SHORT=15 #discard length of sequences
QT_THRESHOLD5=10 #threshold for quality trimming; we filter by number of mismatches so we need high quality reads
QT_THRESHOLD3=10 #threshold for quality filtering
MAX_LENGTH=30
# running the cutadapt tool for all the files
for sample in *${SUFFIX}
do
cutadapt -f $FILE_FORMAT -a $ADAPTER -A $ADAPTER2 -n 1 -q $QT_THRESHOLD5,$QT_THRESHOLD3 --trim-n --max-n=0 -e $ERROR_RATE -O $MIN_OVERLAP -m  $DISC_SHORT -M $MAX_LENGTH -o ${sample%.fastq*}.ad3trim.fq.gz --info-file=${sample%.fastq*}.ad3info --too-short-output=${sample%.fastq*}.ad3short.fq.gz --untrimmed-output=${sample%.fastq*}.ad3untrimmed.fq.gz $sample &>${sample%.fastq*}.cutadapt_output || exit 3
done
cp *.ad3trim.fq.gz $OUTPUT_FOLDER_CUTADAPTgz/     #copy trimmed files
mv *.cutadapt_output $OUTPUT_FOLDER_CUTADAPTgz/   # Adapter output files 
mv *.ad3short.fq.gz $OUTPUT_FOLDER_CUTADAPT/  
mv *.ad3untrimmed.fq.gz $OUTPUT_FOLDER_CUTADAPT/


##FASTQC check after Adapter trimming###
OUTPUT_FOLDER_FASTQC_TRIMMED=/storage/brno2/home/tarakaramji/results/fastqc_trimmed/
mkdir -p $OUTPUT_FOLDER_FASTQC_TRIMMED

fastqc -t 4 *.ad3trim.fq.gz || exit 4 #fastqc command line for the trimmed fastq files
mv *.zip $OUTPUT_FOLDER_FASTQC_TRIMMED/ || exit 5 #fastqc command line for all the folder with .gz extension
mv *.html $OUTPUT_FOLDER_FASTQC_TRIMMED/ || exit 6 #fastqc command line for all the folder with .gz extension

rm -rf $SCRATCH/* 


#CutAdapt parameters reference
# -f type of file (eg. fasta)
# -a 3`adapter seqeunce (-g for 5`)
# -n no of rounds
# -q phred score
# --trim-n trim the N after adapter removal 
# --max-n=0 dicards reads containing more than count N bases (either 0 or 1)
# -e error rate (10% or 0.1) by default for the adapter
# -O  overlap length (by default 3)
# -m throw away shorter than 10 reads or write to file --too-short-outputfile 
# -o output file 
# --info-file= details



#!/bin/bash
###########
#simple script for minion adapter search
#within the folder  (DATASET_FOLDER_NAME) it takes all files with
# specified suffix (SUFFIX) and looks for the adapters and compares
# them with adapter file (ADAPTERS)
#minion can also take .gz input files as itself

ADAPTERS=adapters_merge.txt
#DATASET_FOLDER_NAME=/home/jan/Projects/katka_mirna/2016/ibrutinib/mirnaNotrim/160421_NS500595_0024_AH253GBGXY
OUTPUT_FOLDER=Adapters
SUFFIX=.fastq

mkdir -p $OUTPUT_FOLDER

#cd $DATASET_FOLDER_NAME

for i in *${SUFFIX}
do
	minion search-adapter -i $i -show 3 -write-fasta $OUTPUT_FOLDER/${i%.*}.minion.fasta
	swan -r $ADAPTERS -q $OUTPUT_FOLDER/${i%.*}.minion.fasta > $OUTPUT_FOLDER/${i%.*}.minion.compare
done

#generate the file with sample name and adapters sequence 

for file in *.compare; do echo -n "$file"; sed '3!d' "$file"; done

#!/bin/bash

# 03-21-17 Tara

# Copy raw RNA seq data files into new directory
dest_dir=/bio/tgallagh/EE283/EE283_final/data/raw
for i in $(seq 21 29); do
cp /bio/tgallagh/RNAseq/data/raw_data/raw_reads/Raw_Reads/P$i\_READ1.fastq $dest_dir/P$i\_READ1.fastq

cp /bio/tgallagh/RNAseq/data/raw_data/raw_reads/Raw_Reads/P$i\_READ2.fastq $dest_dir/P$i\_READ2.fastq

done 

#Rename the files so that sample info is included
cd $dest_dir
tail -9  P1_RNAsamples.txt | while IFS=$'\t' read -r new old; do cp -i -- P$old\_READ2.fastq $new\_READ2.fastq; done
tail -9  P1_RNAsamples.txt | while IFS=$'\t' read -r new old; do cp -i -- P$old\_RE[tgallagh@compute-1-13 raw]$ tail -9  P1_RNAsamples.txt | while IFS=$'\t' read -r new old; do cp -i -- P$old\_READ1.fastq $new\_READ1.fastq
#Delete old files 
for i in $(seq 21 29); do rm -f  P$i\_READ* ; done


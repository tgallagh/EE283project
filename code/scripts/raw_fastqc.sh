#!/bin/bash
#$ -N P1_RNAseq_rawdata
#$ -q pub64,free64,bio
#$ -pe make 8 
#$ -R y
#$ -t 1-18

module load fastqc/0.11.5  

#### Script for P1 RNAseq fastqc

## Naming files for array job
cd /bio/tgallagh/EE283/EE283_final/data/raw 
#for x in *.fastq; do echo $x >> RNAseqlist.txt; done

#mkdir /bio/tgallagh/EE283/EE283_final/output/reports/fastqc_raw
input=$(head -n $SGE_TASK_ID RNAseqlist.txt | tail -n 1)
dir_dest=/bio/tgallagh/EE283/EE283_final/output/reports/fastqc_raw

## Generate fastqc reports of the raw reads
fastqc -o $dir_dest $input



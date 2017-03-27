#!/bin/bash
#$ -N P1_RNAseq_pipeline
#$ -q pub64,bio
#$ -m e

### Pipeline for subsampling bam, merging bam, and preparing for IGV
module load samtools/1.3
#mkdir /bio/tgallagh/EE283/EE283_final/data/processed/bowtie2/subsample
SUB=/bio/tgallagh/EE283/EE283_final/data/processed/bowtie2/subsample/
ALIGN=/bio/tgallagh/EE283/EE283_final/data/processed/bowtie2/

for x in $(cat /bio/tgallagh/EE283/EE283_final/data/raw/prefix.txt); do
samtools view -h $ALIGN$x\both.sorted.RG.bam| head -n 20000| samtools view -bS - > $SUB$x\both.sorted.RG.sub.bam
samtools view -h $ALIGN$x\singles.sorted.RG.bam| head -n 20000| samtools view -bS - > $SUB$x\singles.sorted.RG.sub.bam
done
module load bamtools/2.3.0
cd $SUB

bamtools merge -list listbamfinal -out $SUB\all.p1.merged.bam

samtools sort $SUB\all.p1.merged.bam -o $SUB\all.p1.merged.sorted.bam
samtools index $SUB\all.p1.merged.sorted.bam


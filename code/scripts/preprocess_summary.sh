#!/bin/sh

# Tara Script to summarize pre-processing of  RNAseq reads
#output file
module load samtools 
touch /bio/tgallagh/EE283/EE283_final/output/tables/P1_RNAseq_summary.txt
OUTPUT=/bio/tgallagh/EE283/EE283_final/output/tables/P1_RNAseq_summary.txt

echo "#RAW READS SUMMARY" > $OUTPUT
#Determine the number of reads per file
cd /bio/tgallagh/EE283/EE283_final/data/raw
PREFIX=/bio/tgallagh/EE283/EE283_final/data/raw/prefix.txt
for x in $(cat $PREFIX); do 

echo  Raw_Fwd_Read_$x\ $'\t' "$(grep ^"@" <$x\READ1.fastq \
      | wc -l)" >> $OUTPUT;
echo  Raw_Rev_Read_$x $'\t' "$(grep ^"@" <$x\READ1.fastq \
      | wc -l)" >> $OUTPUT
done

#Determine number of reads that pass trimmomatic QF
#change directory
cd /bio/tgallagh/EE283/EE283_final/data/processed/trimmed
echo "#TRIMMOMATIC OUTPUT SUMMARY" >> $OUTPUT

echo "#READS WHERE BOTH MATES PASSES QF" >> $OUTPUT
for x in $(cat $PREFIX); do 
echo  Trimm_Fwd_PairedRead_P$i\ $'\t' "$(grep ^"@" <$x\READ1_paired2.fastq \
      | wc -l)" >> $OUTPUT;
echo  Trimm_Rev_PairedRead_P$i\ $'\t' "$(grep ^"@" <$x\READ2_paired2.fastq \
      | wc -l)" >> $OUTPUT
done

echo "#READS WHERE 1 MATE PASSES QF" >> $OUTPUT

for x in $(cat $PREFIX); do 
echo  Trimm_Fwd_Read_$x\ $'\t' "$(grep ^"@" <$x\READ1_unpaired2.fastq \
      | wc -l)" >> $OUTPUT;
echo  Trimm_Rev_Read_$x\ $'\t' "$(grep ^"@" <$x\READ2_unpaired2.fastq \
      | wc -l)" >> $OUTPUT
done 

#Determine number of reads that are overlapping
echo "#OVERLAPPING PE READS" >> $OUTPUT
for x in $(cat $PREFIX); do
echo  Overlapping_reads_$x\ $'\t' "$(grep ^"@" <$x\READ1_READ2_merged.fastq.assembled.fastq\
| wc -l)" >> $OUTPUT;
done

#Determine number of reads with rRNA contaminants
echo "#FINAL CLEANED READS AFTER RRNA REMOVAL" >> $OUTPUT
TRIMMED=/bio/tgallagh/EE283/EE283_final/data/processed/trimmed/
for x in $(cat $PREFIX); do
cd $TRIMMED$x\READ1_READ2_cat_unpaired2.fasta.deconseq
echo  SE_Final_$x\ $'\t' "$(grep ^">" <*_clean.fa\
| wc -l)" >> $OUTPUT;
cd $TRIMMED$x\READ1_READ2_merged.fastq.unassembled.forward.fasta.deconseq
echo  Forward_Final_$x\ $'\t' "$(grep ^">" <*_clean.fa \
| wc -l)" >> $OUTPUT;
cd $TRIMMED$x\READ1_READ2_merged.fastq.unassembled.reverse.fasta.deconseq
echo  Reverse_Final_$x\ $'\t' "$(grep ^">" <*_clean.fa \
| wc -l)" >> $OUTPUT;
done

### Compare alignments
cd /bio/tgallagh/EE283/EE283_final/data/raw/bowtie2
echo "PERCENT ALIGNMENT OF RAW READS" >> $OUTPUT
for x in $(cat $PREFIX); do
samtools flagstat $x\_raw.sorted.bam >  $x\.stats.temp
echo $x\raw_alignment $(cat $x\.stats.temp | head -5 | tail -1) >> $OUTPUT
done 

cd /bio/tgallagh/EE283/EE283_final/data/processed/bowtie2
echo "PERCENT ALIGNMENT OF PROCESSED READS" >> $OUTPUT
for x in $(cat $PREFIX); do
samtools flagstat $x\.both.sorted.bam  >  $x\.stats.temp
echo $x\processed_alignment $(cat $x\.stats.temp | head -5 | tail -1) >> $OUTPUT
done






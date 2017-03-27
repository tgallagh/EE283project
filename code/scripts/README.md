#### Scripts for RNA seq pre-processing and analysis of P1 data
copy_rename.sh --- move raw data into new directory and rename so condition is included in header

raw_fastqc.sh --- fastqc reports of all 18 raw reads

preprocess_pipeline.sh --- job array of all 9 samples through quality filtering and alignment

merge_igv --- to subsample and merge bam files for viewing on IGV 


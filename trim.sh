#!/bin/bash
#FILENAME: Slater_HB_RNAseq
#SBATCH -A lenders
#SBATCH -n 1
#SBATCH -t 75:00:00

# load the modules
	module load bioinfo
	module load trimmomatic


k=$1
pidlist_sampe=""
endfq="_R1_001.fastq.gz"
endfq2="_R2_001.fastq.gz"
endTP="_1_TP.fastq.gz"
endTP2="_2_TP.fastq.gz"
endTU="_1_TU.fastq.gz"
endTU2="_2_TU.fastq.gz"
fq1=$k$endfq
fq2=$k$endfq2
TP1=$k$endTP
TP2=$k$endTP2
TU1=$k$endTU
TU2=$k$endTU2
trimmomatic PE -threads 9 -phred33 $fq1 $fq2 $TP1 $TU1 $TP2 $TU2 \
ILLUMINACLIP:/depot/bioinfo/apps/apps/trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
echo " "
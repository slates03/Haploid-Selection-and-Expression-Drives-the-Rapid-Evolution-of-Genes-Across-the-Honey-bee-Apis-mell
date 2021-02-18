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
endfq=".fastq.gz"
endTP="_TP.fastq.gz"

fq1=$k$endfq
TP1=$k$endTP



trimmomatic SE $fq1 $TP1  \
ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


echo " "
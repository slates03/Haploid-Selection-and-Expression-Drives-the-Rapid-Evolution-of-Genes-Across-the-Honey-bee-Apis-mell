#!/bin/sh
#FILENAME: kallisto DGE
#SBATCH -A lenders
#SBATCH -n 1
#SBATCH -t 75:00:00

module load bioinfo
module load kallisto

k=$1
pidlist_sampe=""
endfq="_1_TP.fq.gz"
endfq2="_2_TP.fq.gz"
endTP="output_"
fq1=$k$endfq
fq2=$k$endfq2
TP1=$endTP$k

kallisto quant -i EU_transcripts.idx -o $TP1 -b 100 $fq2 $fq1

echo " "
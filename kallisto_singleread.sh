#!/bin/sh
#FILENAME: kallisto DGE
#SBATCH -A lenders
#SBATCH -n 1
#SBATCH -t 75:00:00

module load bioinfo
module load kallisto

k=$1
pidlist_sampe=""
endfq="_TP.fastq.gz"
endTP="output_"
fq1=$k$endfq
TP1=$endTP$k

kallisto quant -i varroa.idx -o $TP1 --single -l 200 -s 20 $fq1


echo " "
# Haploid Selection Protocol
## RNAseq

The RNAseq library was prepared using the Illumina TruSeq stranded mRNA sample preparation kit according to the manufacturer’s guidelines. Poly-A containing RNA was purified from total RNA using poly-T oligo attached magnetic beads after which mRNA will be fragmented and reverse transcribed to first strand cDNA using random primers. The cDNA fragments will be ligated to adapters and purified cDNA libraries enriched with PCR. All sequencing was performed using Illumina HiSeq 2500 sequencing technology producing 150-bp length paired-end reads.

## Trimming
Build a 'PRJNAXXXXX_samps' file that contains the NCBI SRR IDs for your samples (and the ID's of the files you downloaded with the rpevious script). The script below will take those in and trim them with trimmomatic. It assumes you're using True-Seq adaptors.

```
filename=Slater_EU_RNAseq.txt 
exec 4<$filename
echo Start
while read -u4 k ; do
sbatch /scratch/snyder/s/slater20/Shells/trimmomatic.sh "$k"
```
done 

```
#!/bin/sh
# FILENAME: trimmer
#PBS -q lenders
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l naccesspolicy=singleuser

# load the modules
module load bioinfo
module load trimmomatic

# specify commands

cd $PBS_O_WORKDIR
pwd
date +"%d%B%Y%H:%M:%S"

echo " "

k=$1
pidlist_sampe=""
endfq="_1.fastq.gz"
endfq2="_2.fastq.gz"
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
```

## Kallisto
Build Index File

```
module load bioinfo
module load kallisto

kallisto index -i EU_transcripts.idx /depot/bharpur/data/ref_genomes/EDIL/Edil_OGS_v1.1_dna.fa
```

Build a 'PRJNAXXXXX_samps' file that contains the trimmed samples. The script below will take those in and estimate abundances with Kallisto. 

```
filename=Slater_EU_RNAseq.txt 
exec 4<$filename
echo Start
while read -u4 k ; do
sbatch /scratch/snyder/s/slater20/Shells/kallisto.sh "$k"
```

done

```
#!/bin/sh
	# FILENAME: bwa
	#PBS -q lenders
	#PBS -l nodes=1:ppn=1
	#PBS -l walltime=24:00:00
	#PBS -l naccesspolicy=singleuser
  
	module load bioinfo
	module load kallisto

	# specify commands

	cd $PBS_O_WORKDIR
	pwd
	date +"%d%B%Y%H:%M:%S"

	echo " "
	k=$1
	pidlist_sampe=""
	endfq="_R2_001_TAligned.sortedByCoord.out.fastq.gz"
	endfq2="_R1_001_TAligned.sortedByCoord.out.fastq.gz"
	endTP="output_"
	fq1=$k$endfq
	fq2=$k$endfq2
	TP1=$endTP$k

kallisto quant -i transcripts.idx -o $TP1 -b 100 $fq2 $fq1
```

## Differential Gene Expression (DEG) DESeq
###Requirements

```
#Install Bioconductor/tximport
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("tximport")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("SummarizedExperiment")
BiocManager::install("apeglm")
BiocManager::install("ImpulseDE2")
BiocManager::install("devtools")
BiocManager::install("biocLite")
BiocManager::install("polyester")

library(biomaRt)
library(tximport)
library(GEOquery)
library(DESeq2)
library(ImpulseDE2)
library(biocLite)
library(dplyr)
```
### Remove rRNA Function
```
rrna<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/RNAContaminationReads.csv")
```

```
########################################
##                                    ##
##           Filter Files             ##
##                                    ##
########################################
filter_files=function(rrna,Sample) {
  bees = as.character(unique(Sample$Name))
  for(bee in bees){
    test = Sample[which(Sample$Name==bee), ]
    test=test[,1]
    path<-"~/OneDrive/Manuscripts/Haploid Selection/All/"
    file<-paste(path,test,sep="")
    X<-read.table(file,header = TRUE)
    X<-data.frame(X)
    X<-X[!(X$target_id %in% rrna$target_id),]
    
    add<-"2_"
    final<-paste(add,test,sep="")
    final<-paste(path,final,sep="")
    write.table(X, file = final, row.names=FALSE, sep="\t")
    }
  }
```

### Merge Population Genetic Datasets from a single repository- Function

```
########################################
##                                    ##
##            Resource File           ##
##                                    ##
########################################

#Merge Resource Files
multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,by="target_id")}, datalist)}

AMEL<-multmerge("~/OneDrive/Manuscripts/Haploid Selection/PopGenResources")

```
### Run DESeq analysis in R- Function

```
#Parameters
Sample=read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/dev_samp_2.csv") #Sample File 
tx2gene=read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/tx2gene.csv") #Transcript File 

########################################
##                                    ##
##           RNASeq Kallisto          ##
##                                    ##
########################################

#RNASeq with Kallisto
kallisto = function(Sample,tx2gene,AMEL,Type) {
  #Sample Names
  files<-t(Sample[,1])
  files<-as.character(files)
  names(files) <- Sample[,3]
  
  setwd("~/OneDrive/Manuscripts/Haploid Selection/All")
  #Run Kallisto
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  #DESeq Analysis
  dds <- DESeqDataSetFromTximport(txi, Sample, design= Type)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)
}

dds<-kallisto(Sample,tx2gene,AMEL,Type=~Sex)
```
### Create Contrast Matrix (upregulated and all DEG)

```
########################################
##                                    ##
##           Contrast Matrix          ##
##                                    ##
########################################

#Upregulated
contrast_up = function(type,x,y) {
  lar_d <- data.frame(results(dds, contrast=c(type,x,y)))
  lar_d <- data.frame(lar_d[which(lar_d$padj<0.1),])
  lar_d<-lar_d[which(lar_d$log2FoldChange>0),]
  lar_d$target_id<-rownames(lar_d)
  lar_d<-merge(GB_to_NCBI,lar_d,by="target_id")
  lar_d<-lar_d[!duplicated(lar_d$GB), ]
}
  
#All DEGs
contrast = function(type,x,y) {
  lar_d <- data.frame(results(dds, contrast=c(type,x,y)))
  lar_d <- data.frame(lar_d[which(lar_d$padj<0.1),])
  lar_d$target_id<-rownames(lar_d)
  lar_d<-merge(GB_to_NCBI,lar_d,by="target_id")
  lar_d<-lar_d[!duplicated(lar_d$GB), ]
}


lar_dvq<-contrast_up("Sex","drone","queen")
```
### Full Analysis in R - Outline
```
# 1. Remove Contamination
# 2. Set up data
# 3. Larval RNASeq + Contrasts
# 4. Adult RNASeq  + Contrasts
# 5. Produce Gene Lists
# 6. Graphical Representation-Tissue
# 7. Statistics
# 8. Proportion DOS 
# 9. dNdS Alternative Dataset 
# 10.Funtional Analysis

```

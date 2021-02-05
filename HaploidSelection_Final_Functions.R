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

########################################
##                                    ##
##           Contrast Matrix          ##
##                                    ##
########################################

contrast_up = function(type,x,y) {
  lar_d <- data.frame(results(dds, contrast=c(type,x,y)))
  lar_d <- data.frame(lar_d[which(lar_d$padj<0.1),])
  lar_d<-lar_d[which(lar_d$log2FoldChange>0),]
  lar_d$target_id<-rownames(lar_d)
  lar_d<-merge(GB_to_NCBI,lar_d,by="target_id")
  lar_d<-lar_d[!duplicated(lar_d$GB), ]
}
  
contrast = function(type,x,y) {
  lar_d <- data.frame(results(dds, contrast=c(type,x,y)))
  lar_d <- data.frame(lar_d[which(lar_d$padj<0.1),])
  lar_d$target_id<-rownames(lar_d)
  lar_d<-merge(GB_to_NCBI,lar_d,by="target_id")
  lar_d<-lar_d[!duplicated(lar_d$GB), ]
}

########################################
##                                    ##
##              Graphics              ##
##                                    ##
########################################

graph_tissue=function(sum,xaxis,yaxis,min,max) {
ggplot(sum) +
  geom_bar(aes(x=sex, y=mean, fill=tissue), stat="identity", alpha=0.5, width=0.5,position=position_dodge(1)) +
  geom_errorbar(aes(x=sex,ymin=mean-se, ymax=mean+se,group=tissue), position=position_dodge(1),width=0.1, size=.7) +
  xlab(xaxis) +
  ylab(yaxis) +
  scale_colour_hue(name="Tissue Type", l=40) + 
  scale_x_discrete(labels = c("Control","Diploid","Haploid")) +
  scale_y_continuous(expand = c(0, 0),limits=c(min,max)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x   = element_text(size=10, color="black"),
        axis.title.y  = element_text(face="bold", size=12),
        axis.text.y   = element_text(size=10, color="black"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))

}


graph_all=function(sum,xaxis,yaxis,min,max) {
ggplot(sum) +
  geom_bar(aes(x=sex, y=mean), stat="identity", alpha=0.5,width=0.5, position=position_dodge(1)) +
  geom_errorbar(aes(x=sex,ymin=mean-se, ymax=mean+se), position=position_dodge(1),width=0.1, size=.7) +
  xlab(xaxis) +
  ylab(yaxis) +
  scale_colour_hue(name="Tissue Type", l=40) + 
  scale_x_discrete(labels = c("Control","Diploid","Haploid")) +
  scale_y_continuous(expand = c(0,0 ),limits=c(min,max)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x   = element_text(size=10, color="black"),
        axis.title.y  = element_text(face="bold", size=12),
        axis.text.y   = element_text(size=10, color="black"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black")) 
  
}


hist_sfs=function(data,xaxis,yaxis,bin) {
  ggplot(sel_s,aes(x = freq, fill = sex, color = sex))+ 
    geom_histogram(aes(y=..density..),position="dodge", alpha=0.5,bins=bin) +
    xlab(xaxis) +
    ylab(yaxis) +
    theme_bw() +
    theme(axis.title.x = element_text(face="bold", size=12),
          axis.text.x   = element_text(size=10, color="black"),
          axis.title.y  = element_text(face="bold", size=12),
          axis.text.y   = element_text(size=10, color="black"),
          legend.background = element_rect(size=0.5, linetype="solid",colour ="black")) 
  
}



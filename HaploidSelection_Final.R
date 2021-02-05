library(dplyr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(ggpubr)
############################################################################
#== OUTLINE
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
############################################################################
########################################
##                                    ##
##     1.  Remove Contamination       ##
##                                    ##
########################################
rrna<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/RNAContaminationReads.csv")
Sample<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/all_samp.csv")
XX<-filter_files(rrna,Sample)

########################################
##                                    ##
##      2.    Setup Data              ##
##                                    ##
########################################

#PopGen Resource File
AMEL<-multmerge("~/OneDrive/Manuscripts/Haploid Selection/PopGenResources")
GB_to_NCBI<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/PopGenResources/GB_to_NCBI.csv")

########################################
##                                    ##
##     3.    Larval RNAseq            ##
##                                    ##
########################################
#Larval Data
#lar_d = drone larvae
#lar_q = queen larvae
#Parameters
Sample=read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/dev_samp_2.csv") #Sample File 
tx2gene=read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/tx2gene.csv") #Transcript File 

#RNASeq with Kallisto
dds<-kallisto(Sample,tx2gene,AMEL,Type=~Sex)

#Contrast Matrix-Drone Larval
lar_dvq<-contrast_up("Sex","drone","queen")
lar_dvw<-contrast_up("Sex","drone","worker")
lar_d<-merge(lar_dvq,lar_dvw,by="GB")
colnames(lar_d)[2]<-"target_id"

#In order to Merge Datafile later
lar_d$DT<-0;lar_d$DB<-0;lar_d$QB<-0;lar_d$QO<-0;lar_d$sex<-"male";lar_d$tissue<-"male.larval"
lar_d<-dplyr::select(lar_d,target_id,GB,DT,DB,QB,QO,sex,tissue)
#Remove any duplicates
lar_d<-lar_d[!duplicated(lar_d$GB), ]

#Contrast Matrix-Queen Larval
lar_qvd<-contrast_up("Sex","queen","drone")
lar_qvw<-contrast_up("Sex","queen","worker")
lar_q<-merge(lar_qvd,lar_qvw,by="GB")
colnames(lar_q)[2]<-"target_id"

#In order to Merge Datafile later
lar_q$DT<-0;lar_q$DB<-0;lar_q$QB<-0;lar_q$QO<-0;lar_q$sex<-"female";lar_q$tissue<-"female.larval"
lar_q<-dplyr::select(lar_q,target_id,GB,DT,DB,QB,QO,sex,tissue)
#Remove any duplicates
lar_q<-lar_q[!duplicated(lar_q$GB), ]
lar<-rbind(lar_d,lar_q)
table(lar$sex)
########################################
##                                    ##
##        4.    Adult RNAseq          ##
##                                    ##
########################################
#Adult Data
#DT=Drone Testis
#DB=Drone Brain
#QO=Queen Ovaries
#QB=Queen Brain

#Parameters
Sample=read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/adult_samp.csv") #Sample File 
tx2gene=read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/tx2gene.csv") #Transcript File 

#RNASeq with Kallisto
dds<-kallisto(Sample,tx2gene,AMEL,Type=~Final)

#Contrast Matrix-Drone Adult
#Drone testis
DTvDB<-contrast_up("Final","DG","DB")
DTvQB<-contrast_up("Final","DG","QB")
DTvQO<-contrast_up("Final","DG","QG")
DT<-merge(DTvDB,DTvQB,by="GB")
DT<-merge(DT,DTvQO,by="GB")
#For Upset Plot
DT$DT<-1
DT<-dplyr::select(DT,GB,DT)

#Drone Brain
DBvDT<-contrast_up("Final","DB","DG")
DBvQB<-contrast_up("Final","DB","QB")
DBvQO<-contrast_up("Final","DB","QG")
DB<-merge(DBvDT,DBvQB,by="GB")
DB<-merge(DB,DBvQO,by="GB")
#For Upset Plot
DB$DB<-1
DB<-dplyr::select(DB,GB,DB)

#Queen Brain
QBvDT<-contrast_up("Final","QB","DG")
QBvDB<-contrast_up("Final","QB","DB")
QBvQO<-contrast_up("Final","QB","QG")
QB<-merge(QBvDT,QBvDB,by="GB")
QB<-merge(QB,QBvQO,by="GB")
#For Upset Plot
QB$QB<-1
QB<-dplyr::select(QB,GB,QB)

#Queen Ovaries
QOvDT<-contrast_up("Final","QG","DG")
QOvDB<-contrast_up("Final","QG","DB")
QOvQB<-contrast_up("Final","QG","QB")
QO<-merge(QOvDT,QOvDB,by="GB")
QO<-merge(QO,QOvQB,by="GB")
#For Upset Plot
QO$QO<-1
QO<-dplyr::select(QO,GB,QO)

#Final Gene Set
final<-merge(DT,DB,all=TRUE,by="GB")
final<-merge(final,QB,all=TRUE,by="GB")
final<-merge(final,QO,all=TRUE,by="GB")
final[is.na(final)] <- 0

#upset Plot
final<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/upset_final.csv")

library(UpSetR)
upset_final<-upset(final,order.by=c("freq"))
upset_final

########################################
##                                    ##
##       5.  Produce Gene Lists       ##
##                                    ##
########################################

#Identify and Extract Gene Lists
#DEGs:
hap_1 <- subset(final,final$DT==1 & final$DB==1 & final$QO==0 & final$QB==0)
hap_1$sex<-c("male")
hap_1$tissue<-c("male.both")
hap_1_1 <- subset(final,final$DT==1 & final$DB==1 & final$QO==0 & final$QB==0)
hap_1_1$sex<-c("male")
hap_1_1$tissue<-c("male.somatic")
hap_2 <- subset(final,final$DT==1 & final$DB==0 & final$QO==0 & final$QB==0)
hap_2$sex<-c("male")
hap_2$tissue<-c("male.gonadal")
hap_3<- subset(final,final$DT==0 & final$DB==1 & final$QO==0 & final$QB==0)
hap_3$sex<-c("male")
hap_3$tissue<-c("male.somatic")
dip_1<-subset(final,final$DT==0 & final$DB==0 & final$QO==1 & final$QB==1)
dip_1$sex<-c("female")
dip_1$tissue<-c("female.both")
dip_1_2$sex<-c("female")
dip_1_2$tissue<-c("female.somatic")
dip_2<-subset(final,final$DT==0 & final$DB==0 & final$QO==1 & final$QB==0)
dip_2$sex<-c("female")
dip_2$tissue<-c("female.gonadal")
dip_3<-subset(final,final$DT==0 & final$DB==0 & final$QO==0 & final$QB==1)
dip_3$sex<-c("female")
dip_3$tissue<-c("female.somatic")

########################
#Hap,Dip and Con gene list
########################
#Remove column to equal other
lar_d<-lar_d[,-1]
lar_q<-lar_q[,-1]

#Remove Duplicates
ma<-rbind(hap_1,hap_2,hap_3,lar_d)
ma<-ma %>% distinct(GB,.keep_all = TRUE)

fe<-rbind(dip_1,dip_2,dip_3,lar_q)
fe<-fe %>% distinct(GB,.keep_all = TRUE)

sel<-rbind(fe,ma)

#Removee NAs
sel<-rbind(ma,fe)
sel<-merge(sel,AMEL,by="GB",all=TRUE)

##########################
#Tissue Specific Gene List
#########################

sel_tis<-sel[!sel$tissue == "male.both", ]
sel_tis<-sel_tis[!sel_tis$tissue == "female.both", ]

write.csv(sel,"~/OneDrive/Manuscripts/Haploid Selection/Resources/haploid_final.csv")
sel<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/haploid_final.csv")
            
write.csv(sel_tis,"~/OneDrive/Manuscripts/Haploid Selection/Resources/haploid_tissue_final.csv")
sel_tis<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/haploid_tissue_final.csv")

#Final Gene list

list_all<-sel%>%group_by(sex)%>%dplyr::summarize(number=length(GB))
list_tissue<-sel%>%group_by(tissue)%>%dplyr::summarize(number=length(GB))

########################################
##                                    ##
##         6.  Graphics               ##
##                                    ##
########################################


#Graphics-Tissue Comparisons
sum = sel_tis[!is.na(sel_tis$piN), ] %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(piN),sd=sd(piN),N=length(piN),se=sd/sqrt(N),na.rm=TRUE)
piN=graph_tissue(sum,"geneType","piN",min=0,max=0.00045)
sum = sel_tis[!is.na(sel_tis$piS), ] %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(piS),sd=sd(piS),N=length(piS),se=sd/sqrt(N),na.rm=TRUE)
piS<-graph_tissue(sum,"geneType","piS",min=0,max=0.006)
sum = sel_tis[!is.na(sel_tis$piNpiS), ] %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(piNpiS),sd=sd(piNpiS),N=length(piNpiS),se=sd/sqrt(N))
piNpiS<-graph_tissue(sum,"geneType","piNpiS",min=0,max=0.15)
sum = sel_tis[!is.na(sel_tis$dos), ] %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(dos),sd=sd(dos),N=length(dos),se=sd/sqrt(N))
dos<-graph_tissue(sum,"geneType","DoS",min=-0.03,max=0.045)
sum = sel_tis[!is.na(sel_tis$rate), ] %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(rate),sd=sd(rate),N=length(rate),se=sd/sqrt(N))
rate<-graph_tissue(sum,"geneType","Recombination Rate",min=0,max=22)
sum = sel_tis[!is.na(sel_tis$GCperc), ] %>% group_by(sex,tissue) %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(GCperc),sd=sd(GCperc),N=length(GCperc),se=sd/sqrt(N))
GC<-graph_tissue(sum,"geneType","GC Percentage",min=0,max=0.45)
sum = sel_tis[!is.na(sel_tis$FstCvS), ] %>% group_by(sex,tissue) %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(FstCvS),sd=sd(FstCvS),N=length(FstCvS),se=sd/sqrt(N))
FST<-graph_tissue(sum,"geneType","FST MvsS",min=0,max=0.17)

diver<-ggarrange(piN,piS,ncol=2,nrow=1,labels=c("A","B"))
selection<-ggarrange(piNpiS,dos,ncol=2,nrow=1,labels=c("A","B"))
rate_GC<-ggarrange(GC,rate,ncol=2,nrow=1,labels=c("A","B"))

#Graphics-Non-tissue comparisons
sum = sel[!is.na(sel$piN), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(piN),sd=sd(piN),N=length(piN),se=sd/sqrt(N))
piN<-graph_all(sum,"geneType","piN",min=0,max=0.00045)
sum = sel[!is.na(sel$piS), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(piS),sd=sd(piS),N=length(piS),se=sd/sqrt(N))
piS<-graph_all(sum,"geneType","piS",min=0,max=0.006)
sum = sel[!is.na(sel$piNpiS), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(piNpiS),sd=sd(piNpiS),N=length(piNpiS),se=sd/sqrt(N))
piNpiS<-graph_all(sum,"geneType","piNpiS",min=0,max=0.15)
sum = sel[!is.na(sel$rate), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(rate),sd=sd(rate),N=length(rate),se=sd/sqrt(N))
rate<-graph_all(sum,"geneType","Recombination Rate",min=0,max=22)
sum = sel[!is.na(sel$GCperc), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(GCperc),sd=sd(GCperc),N=length(GCperc),se=sd/sqrt(N))
GC<-graph_all(sum,"geneType","GC Percentage",min=0,max=0.45)
sum = sel[!is.na(sel$dos), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(dos),sd=sd(dos),N=length(dos),se=sd/sqrt(N))
dos<-graph_all(sum,"geneType","DoS",min=-0.03,max=0.04)
sum = sel[!is.na(sel$FstCvS), ] %>% group_by(sex) %>% dplyr::summarize(mean=mean(FstCvS),sd=sd(FstCvS),N=length(FstCvS),se=sd/sqrt(N))
FST<-graph_all(sum,"geneType","sel$FstCvM",min=-0.0,max=.15)

diver<-ggarrange(piN,piS,ncol=2,nrow=1,labels=c("A","B"))
selection<-ggarrange(piNpiS,dos,ncol=2,nrow=1,labels=c("A","B"))
rate_GC<-ggarrange(GC,rate,ncol=2,nrow=1,labels=c("A","B"))

########################################
##                                    ##
##          7.   Statistics           ##
##                                    ##
########################################
#ACNOVA with GCperc as covariate
res.aov <- aov(rate~ GCperc + sex, data = larval)
summary(res.aov)

#ANOVA without GCperc as covariate
res.aov <- aov(rate~ sex, data = larval)
summary(res.aov)


########################################
##                                    ##
##         8. Proportion DOS          ##
##                                    ##
########################################

#Compare DOS between haploid and diploid selected genes using Chisquare 
dos<-sel[which(sel$dos>0 & sel$mk.res<0.05),]

#Put into a chisquare table
pos<-table(dos$sex)
all<-table(sel$sex)
list<-rbind(pos,all)
list<-list[,-1]
chisq.test(list)

#Compare DOS between haploid and diploid tissue-specific genes using Chisquare 
dos<-sel_tis[which(sel_tis$dos>0 & sel_tis$mk.res<0.05),]
pos<-table(dos$tissue)
all<-table(sel_tis$tissue)

#Tissue-Gonadal
pos_1<-pos[c("female.gonadal","male.gonadal")]
all_1<-all[c("female.gonadal","male.gonadal")]
list<-rbind(pos_1,all_1)
chisq.test(list)
(34*1292)/(40*671)
#results:1.63662

#Tissue-Somatic
pos_1<-pos[c("female.somatic","male.somatic")]
all_1<-all[c("female.somatic","male.somatic")]
list<-rbind(pos_1,all_1)
chisq.test(list)
(21*619)/(29*388)
#results: 1.155261

#Tissue-Larval
pos_1<-pos[c("female.larval","male.larval")]
all_1<-all[c("female.larval","male.larval")]
list<-rbind(pos_1,all_1)
chisq.test(list)
(93*1169)/(41*1651)
#results:1.606078

########################################
##                                    ##
##   9.  dNdS Alternative Dataset     ##
##                                    ##
########################################
#Get dnds file and merge with gene list
dnds<-read.csv("~/OneDrive/Manuscripts/Haploid Selection/Resources/dnds.csv")
dnds_sel<-merge(dnds,sel,by="GB")

### Non-tissue specific ####
#Summarize Werner dataset
sum = dnds_sel %>% group_by(sex) %>% group_by(sex) %>% dplyr::summarize(mean=mean(dN_dS_Wern),sd=sd(dN_dS_Wern),N=length(dN_dS_Wern),se=sd/sqrt(N))
wern_all<-graph_all(sum,"geneType","dNdS Wern",min=0,max=0.17)
#Summarize Kapheim dataset
sum = dnds_sel %>% group_by(sex) %>% group_by(sex) %>% dplyr::summarize(mean=mean(dN_dS_Kap),sd=sd(dN_dS_Kap),N=length(dN_dS_Kap),se=sd/sqrt(N))
Kap_all<-graph_all(sum,"geneType","dNdS Kap",min=0,max=0.17)

###   tissue specific   ####
dnds_sel_tis<-merge(dnds,sel_tis,by="GB")
#Summarize Werner dataset
sum = dnds_sel_tis %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(dN_dS_Wern),sd=sd(dN_dS_Wern),N=length(dN_dS_Wern),se=sd/sqrt(N))
wern_tis<-graph_tissue(sum,"geneType","dNdS Wern Tis",min=-0.0,max=.15)

#Summarize Kaphem dataset
sum = dnds_sel_tis %>% group_by(sex,tissue) %>% dplyr::summarize(mean=mean(dN_dS_Kap),sd=sd(dN_dS_Kap),N=length(dN_dS_Kap),se=sd/sqrt(N))
Kap_tis<-graph_tissue(sum,"geneType","dNdS Kap Tis",min=-0.0,max=.2)


########################################
##                                    ##
##         10.  Function              ##
##                                    ##
########################################

#Genes with known male reproductive function
drone<-read.csv("~/OneDrive/Projects/Protocols/Bioinformatics/Proteomics/AMEL_DroneSpecificGenes_2.csv")
drone<-drone%>%filter(drone$ID %in% c("spermatogenesis","seminal_fluid","seminal_vesicle","sperm_competition","sperm_motility","stored_sperm","sv_and_sperm"))
#Extract genes from gene list with known function
drone_1<-merge(drone,sel,by="GB")
drone_1<-drone_1[which(drone_1$sex=="male"),]
#Control Genes-Diploid selected genes
control<-sel[which(sel$sex=="female"),]
control$ID<-c("Diploid")
drone_1<-rbind.fill(drone_1,control)
#Statistics
comp <- subset(drone_1, ID=="Diploid" | ID=="sperm_competition")
spermato <- subset(drone_1, ID=="Diploid" | ID=="spermatogenesis")
fluid <- subset(drone_1, ID=="Diploid" | ID=="seminal_fluid")
vesi <- subset(drone_1, ID=="Diploid" | ID=="seminal_vesicle")
moti <- subset(drone_1, ID=="Diploid" | ID=="sperm_motility")
store <- subset(drone_1, ID=="Diploid" | ID=="stored_sperm")
sv_and_sperm<- subset(drone_1, ID=="Diploid" | ID=="sv_and_sperm")

#Graphical Representation of Function
drone_2 = drone_1[!is.na(drone_1$dos), ] %>% group_by(ID) %>% dplyr::summarize(mean=mean(dos),sd=sd(dos),N=length(dos),se=sd/sqrt(N))

ggplot(drone_1, aes(x=ID, y=dos)) + 
  geom_boxplot() + 
  geom_jitter(aes(color = ifelse(mk.res < 0.05, "red", "black")),size=2.5) +
  xlab("ID") +
  ylab("mean") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x   = element_text(size=10, color="black"),
        axis.title.y  = element_text(face="bold", size=12),
        axis.text.y   = element_text(size=10, color="black"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black")) 


##############################################
#Chi Square analysis of each functional group
##############################################

#Chi Square Analyis of Functional Genes over 0 dos
#Filter Genes with DOS over 1 and are significant (mk.res<0.05). Place into table
dos<-drone_1[which(drone_1$dos>0 & drone_1$mk.res<0.05),]
#Genes Under positive selection
dos<-table(dos$ID)
#All genes
all<-table(drone_1$ID)

#Summarize Table for analysis
dos<-dos%>%group_by(ID)%>%dplyr::summarize(number=length(GB))
dos<-data.frame(dos)

list_all<-drone_1%>%group_by(ID)%>%dplyr::summarize(number=length(GB))
list_all<-data.frame(list_all)

#Seminal Fluid
matrix<-matrix(c(92,3,2619,26),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 2.3258, df = 1, p-value = 0.1272

#Seminal Vesicle
matrix<-matrix(c(92,0,2619,17),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 0.0097632, df = 1, p-value = 0.9213


#Competition
matrix<-matrix(c(92,0,2619,7),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 5.7001e-26, df = 1, p-value = 1

#Motility
matrix<-matrix(c(92,1,2619,7),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 1.7416e-26, df = 1, p-value = 1

#Spermatogenesis
matrix<-matrix(c(92,1,2619,29),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 5.3854e-30, df = 1, p-value = 1

#Stored Sperm
matrix<-matrix(c(92,21,2619,93),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 60.482, df = 1, p-value = 7.424e-15

#sv and Sperm
matrix<-matrix(c(92,2,2619,10),ncol=2,byrow = TRUE)
colnames(matrix) <- c("cont","trt")
rownames(matrix) <- c("pos","no_pos")
matrix <- as.table(matrix)
chisq.test(matrix)
#X-squared = 2.9606, df = 1, p-value = 0.08532




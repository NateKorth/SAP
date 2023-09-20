##Library packages (install.packages() if you need to)
library(sommer)
library(ggplot2)
library(readxl)
library(data.table)
library(dplyr)
library(tidyverse)
library(bestNormalize)

#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/SAP/GWAS")

df1<-read.csv("SAP_S1_FamGenASVAlpha_RA_filtered2.csv")


#Traits<-as.data.frame(names(df1[100:258]))
#S1:
Traits<-as.data.frame(names(df1[255:258]))
#S2:
#Traits<-as.data.frame(names(df1[259:262]))

names(Traits)<-"Traits"

A<-read.csv("Centered_IBS_SAP.csv",check.names = FALSE)
A[1:5,1:5]
rownames(A)<-A$Line
A$Line<-NULL
A<-as.matrix(A)

#Prep data for mapping
#calculate BLUPS:

missing<-df1[-which(df1$PI %in% rownames(A)),]
missing1<-A[-which(rownames(A) %in% df1$PI),]

#Missing from genotype:
#PI533998
#PI534104
#PI574455
#PI656030
#PI656117
#PI659693

#Old SNP file:
df1<-subset(df1,PI!="PI_656081")

#NEW SNP file:
#df1<-subset(df1,PI!="PI_533998")
#df1<-subset(df1,PI!="PI_534104")
#df1<-subset(df1,PI!="PI_574455")
#df1<-subset(df1,PI!="PI_656030")
#df1<-subset(df1,PI!="PI_656117")
#df1<-subset(df1,PI!="PI_659693")

#make factors characters for random effx
df1$Batch<-as.character(df1$Batch)
df1$Plate<-as.character(df1$Plate)
df1$Column<-as.character(df1$Column)


#Calculate h2
fitH <-mmer(S1_PreIndT1~PC1+PC2+PC3, random=~vs(PI,Gu=A)+Batch+Batch:Plate+Plate:Row+Plate:Column, rcov=~units, data=df1) 
summary(fitH)
VC <- summary(fitH)$varcomp
VL <- VC[1,1]
VE<-VC[6,1]
h2<-(VL/(VL+VE))
h2List<-as.data.frame(cbind("Observed",h2))

names(h2List)<-c("Trait","h2")

#Heritability:
geth2<-function(genus){
  t <- formula(paste0(genus, '~PC1+PC2+PC3'))
  fit2<-mmer(t, random=~vs(PI,Gu=A)+Batch+Plate:Batch+Plate:Row+Plate:Column, rcov=~units, data= df1)
  VC <- summary(fit2)$varcomp
  VL <- VC[1,1]
  VE<-VC[6,1]
  h2<-(VL/(VL+VE))
  h2_df <- as.data.frame(cbind(genus,h2))
  names(h2_df)<-c("Trait","h2")
  return(h2_df)
}

#Make sure geth2 function above is in working order
for (i in Traits$Traits){
  h2List<-rbind(h2List,geth2(i))
}

#Write dataframe
#write.csv(h2List,"S2_SAP_h2_estimations.csv")

#If code crashes Read in h2 list:
#h2List<-read.csv("S2_SAP_h2_estimations.csv")

#remove low heritability taxa:
h2List2<-as.data.frame(h2List[-which(h2List$h2<=0.05  ),])
h2List02<-as.data.frame(h2List[which(h2List$h2<=0.05  ),])
Traits1<-as.data.frame(Traits[which(Traits$Traits %in% h2List2$Trait),])
names(Traits1)<-"Traits"


Traits1a<-as.data.frame(Traits1[1:50,])
names(Traits1a)<-"Traits"
Traits1b<-as.data.frame(Traits1[51:107,])
names(Traits1b)<-"Traits"

#Remove outliers
df2<-df1
X<-df2$Shannon
mean<-mean(X)
std<-sd(X)
X[X<(mean-5*std)]<-NA
X[X>(mean+5*std)]<-NA
df3<-as.data.frame(X)

rmOut<-function(OTU){
  X<-df2[[paste0(OTU)]]
  mean<-mean(X)
  std<-sd(X)
  X[X<(mean-5*std)]<-NA
  X[X>(mean+5*std)]<-NA
  return(X)
}

for (i in Traits$Traits){ 
  df3<-cbind(df3,rmOut(i))
}
summary(df3)
sum(is.na(df3))
df3$X<-NULL
names(df3)<-Traits$Traits

#rmOut<-function(OTU){
#  X<-OTU
#  mean<-median(X)
#  std<-sd(X)
#  X[X<(mean-2*std)]<-NA
#  X[X>(mean+2*std)]<-NA
#  return(X)
#}

replace_out <- function(column) {
  qnt <- quantile(column, probs=c(.25, .75),na.rm=TRUE)
  upper_whisker <- 2 * IQR(column,na.rm=TRUE)
  clean_data <- column
  clean_data[column < (qnt[1] - upper_whisker)] <- NA
  clean_data[column > (qnt[2] + upper_whisker)] <- NA
  clean_data
}


#Remove outliers within lines:
df3.1<-as.data.frame(cbind(df2$PI,df3))
names(df3.1)<-c("PI",names(df3))

df3.o <- df3.1 %>% group_by(PI) %>% mutate_if(is.numeric, replace_out)
sum(is.na(df3.1))
sum(is.na(df3.o))


#for (i in Traits$Traits){ 
#  df3<-cbind(df3,rmOut(i))
#}



#Export All traits for autoencoder analysis
#df4<-cbind(df1[14],df3.o)
#write.csv(df4,"MicrobiomeTraits_S1.csv")

df4<-cbind(df1[1:100],df3.o)

#Calculate BLUEs
fitB <-mmer(S1_PreIndT1~PI, random=~Batch+Batch:Plate+Plate:Row+Plate:Column, rcov=~units, data=df4) 
summary(fitB)
BLUEs2<-predict.mmer(object = fitB,classify = "PI")
BLUEs3<-BLUEs2$pvals[,2:3]
mean<-mean(BLUEs3$predicted.value)
std<-sd(BLUEs3$predicted.value)
max<-mean+4*std
min<-mean-4*std
BLUEs4<-BLUEs3
BLUEs4[BLUEs4<(mean-4*std)]<-NA
BLUEs4[BLUEs4>max]<-NA
histogram(BLUEs3$predicted.value)
histogram(BLUEs4$predicted.value)
histogram(df1$Faecalibacterium)
BLUEsdf<-as.data.frame(BLUEs4)
names(BLUEsdf)<-c("PI","observed")
BLUEsdf$PI<-BLUEs3$PI

#df4<-cbind(df1[1:100],df3)

getBLUE <- function(OTU){
  f <- formula(paste0(OTU, '~PI'))
  fit <- mmer(f , random=~Batch+Plate:Batch, rcov=~units, data=df4) 
  BLUEs2<-predict.mmer(object = fit,classify = "PI")
  BLUEs4<-BLUEs2$pvals[,3]
  mean<-mean(BLUEs4)
  std<-sd(BLUEs4)
  BLUEs4[BLUEs4<(mean-4*std)]<-NA
  BLUEs4[BLUEs4>(mean+4*std)]<-NA
  BLUEs5<-as.data.frame(BLUEs4)
  names(BLUEs5)<-paste0(OTU)
  return(BLUEs5)
  gc()
}

#+Plate:Row+Plate:Column

#Calculate BLUPs for heritable taxa:


for (i in Traits$Traits){ 
  BLUEsdf<-cbind(BLUEsdf,getBLUE(i))
}

#for (i in Traits1b$Traits){ 
#  BLUEsdf<-cbind(BLUEsdf,getBLUE(i))
#}


BLUEsdf2<-BLUEsdf
BLUEsdf2$PI<-NULL
BLUEsdf2$Line<-BLUEs3$PI

#Traits3<-Traits1$`Traits[which(Traits$Traits %in% names(h2List2), ]`


#names(BLUEsdf2)<-c("PI","observed",Traits$Traits)
BLUEsdf2$observed<-NULL

#BLUEsdf2a<-BLUEsdf2[1:36]
#BLUEsdf2b<-BLUEsdf2[37:72]
#BLUEsdf2c<-BLUEsdf2[73:103]

#list<-names
#BLUEsdf2b$Line<-BLUEsdf2a$Line
#BLUEsdf2c$Line<-BLUEsdf2a$Line

#BLUEsdf2a.2<-as.data.frame(BLUEsdf2a[-c(1)])
#BLUEsdf2a.3

#write files:
write.csv(BLUEsdf2,"S1_PreIndT_BLUEs.csv",quote = FALSE,row.names = TRUE)

#write.csv(BLUEsdf2a,"S1_SAP_Phenotypes2.1.csv",quote = FALSE,row.names = TRUE)
#write.csv(BLUEsdf2b,"S1_SAP_Phenotypes2.2.csv",quote = FALSE,row.names = TRUE)
#write.csv(BLUEsdf2c,"S1_SAP_Phenotypes2.3.csv",quote = FALSE,row.names = TRUE)

#Phenotypes1model: fitB <-mmer(Observed~PI, random=~Batch+Batch:Plate+Plate:Row+Plate:Column, rcov=~units, data=df1)

#Plot h2 by abundance
h2RA<-read.csv("S1_SAP_h2_estimations2.csv")
h2RA$logRA<-log10(h2RA$Ave_RA)
ggplot(h2RA,aes(x=logRA,y=h2))+geom_point(size=2)+geom_smooth(method="lm",size=2,color="firebrick3")+theme_bw()+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+ theme(legend.position="none",axis.line = element_line(size=1.5)) 


################################################################
################################################################

#VariancePartitioning:
#Wait for Tannin Measurments:





#######################################3
#Pull completed SNP data for plotting:
#Plot single Result
Results<-read.csv("./SNPs_Out/S2_Pheno1/Faecalibacterium_RMIP.csv")
Results1<-data.frame(do.call('rbind', strsplit(as.character(Results$Var1),'_',fixed=TRUE)))
Results2<-data.frame(do.call('rbind', strsplit(as.character(Results1$X1),'S',fixed=TRUE)))
gwasResults<-as.data.frame(cbind(as.numeric(Results2$X2),as.numeric(Results1$X2),as.numeric(Results$Freq)))
gwasResultsMeta<-as.data.frame(cbind(as.numeric(Results2$X2),as.numeric(Results1$X2),as.numeric(Results$Freq),"SeedAcidFiber_T","Seed"))
names(gwasResults)<-c("CHR","BP","FREQ")
names(gwasResultsMeta)<-c("CHR","BP","FREQ","Trait","Type")

sig<-15
data_cum <- gwasResults %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHR, bp_add)

gwasResults <- gwasResults %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

axis_set <- gwasResults %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))


ggplot(gwasResults, aes(x = bp_cum, y = FREQ, color = as.factor(CHR))) +
  geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") + geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_color_manual(values = rep(c("Firebrick","black"), unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  theme_minimal() + theme(legend.position = "none", panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,1,1,1),"cm"))

#Combine all results:
GWAS<-data.frame()

CombineGWAS <- function(trait,type){
  Results<-read.csv(paste0("./SNPs_Out/",type,"/",trait,"_RMIP.csv"))
  Results1<-data.frame(do.call('rbind', strsplit(as.character(Results$Var1),'_',fixed=TRUE)))
  Results2<-data.frame(do.call('rbind', strsplit(as.character(Results1$X1),'S',fixed=TRUE)))
  gwasResults<-as.data.frame(cbind(as.numeric(Results2$X2),as.numeric(Results1$X2),as.numeric(Results$Freq)))
  gwasResultsMeta<-as.data.frame(cbind(as.numeric(Results2$X2),as.numeric(Results1$X2),as.numeric(Results$Freq),paste0(trait),paste0(type)))
  names(gwasResultsMeta)<-c("CHR","BP","FREQ","Trait","Type")
  names(gwasResults)<-c("CHR","BP","FREQ")
  return(gwasResultsMeta)
}

Traits<-read.csv(file ="Traits.csv")
SeedTraits<-as.data.frame(Traits$Seed2[1:46])
BiochemicalTraits<-as.data.frame(Traits$Biochemical[1:43])
AgrinomicTraits<-as.data.frame(Traits$Agrinomic[1:89])
S1Traits<-as.data.frame(Traits$S1.2)

for (i in SeedTraits$Traits){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"Seed"))
}

for (i in BiochemicalTraits$Traits){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"biochemical"))
}
#for (i in AgrinomicTraits$Traits){ 
#  GWAS<-rbind(GWAS,CombineGWAS(i,"Agrinomic"))
#}
for (i in S1Traits$Traits){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"Subject1"))
}
#for (i in S2Traits$Traits){ 
#  GWAS<-rbind(GWAS,CombineGWAS(i,"Subject 2"))
#}


GWAS1<-as.data.frame(sapply(GWAS[1:3],as.numeric))
GWASM<-GWAS[4:5]

data_cum <- GWAS1 %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHR, bp_add)

GWAS1 <- GWAS1 %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

axis_set <- GWAS1 %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

GWAS2<-cbind(GWAS1,GWASM)
sig<-15

ggplot(GWAS2, aes(x = bp_cum, y = FREQ, color = Type)) +
  geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") + geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_color_manual(values = rep(c("Firebrick","DarkGreen","Purple","black"), unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,1,1,1),"cm"))





#TEST GWAS Loop:
ph<-as.data.frame(cbind(BLUEsdf$PI,BLUEsdf2a[,"Faecalibacterium"]))
set.seed(40)
for (i in 3:10){
  z <- sample(1:344,33)
  ph[,i] <- ph[,2]
  ph[z,i] <- NA
  rm(z)
}

length(phenolist)

############################################
#Plot various single traits
#Import TanninmMeasurments:
BC<-read.csv("phenoinput/SAP_BiochemicalPhenotypes.csv")
tan1<-BC[45:46]

dfP<-df4
dfT<-df4
dfP<-subset(df4, Pop!="NA")
dfT<-subset(df4, Tannins_D!="NA")
dfT$Tannins_D<-as.character(dfT$Tannins_D)
dfT$Tannins_E<-as.numeric(dfT$Tannins_E)

dfT1<-merge(dfP,tan1)

dfT$Polyphenols_E
r<-ggplot(data=dfT, aes(x=Tannins_D,y=Butyrate,fill=Tannins_D,na.rm=TRUE))+  geom_boxplot()+scale_fill_manual(values = c("Brown4","cornflowerblue","darkorchid","Red"))+
  ylab("Butyrate")+theme_classic()+theme(text=element_text(size=20,family = "sans"))+#ggtitle("Butyrate")+
  scale_color_manual(values=c("orange", "black","Tan","yellow","red","green4","gold2", "grey"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
r+stat_compare_means()

r<-ggplot(data=dfT1, aes(x=Tannins_D,y=Tannins_NK,fill=Tannins_D,na.rm=TRUE))+  geom_boxplot()+scale_fill_manual(values = c("firebrick3","cornflowerblue","forestgreen","darkorchid"))+
  ylab("Relative Abundance")+theme_classic()+theme(text=element_text(size=20,family = "sans"))+ggtitle("Butyrate")+
  scale_color_manual(values=c("red", "cornflowerblue","green4","violet","red","green4","gold2", "grey"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())#+
  #geom_jitter(aes(color=Tannins_D),height = 0,width=.1)
r+stat_compare_means()
dfT1$Escherichia.Shigella
d<-ggplot(data=dfT1,aes(x=Acetate,y=Butyrate))+scale_color_manual(values = c("firebrick2","cornflowerblue"))+geom_point(size=2)+geom_smooth(method = "lm",se=FALSE,size=1.5)+
  theme_classic()+theme(text=element_text(size=20,family = "sans"))+xlab("ASV3_Clostridium.sensu.stricto.1")+ylab("Butyrate")+
  stat_cor(method="spearman",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+ theme(legend.position="none",axis.line = element_line(size=1.5)) #+scale_x_continuous(breaks=c(.01,.02,.03))
d

#Make Correlation Plots:
library(corrplot)
#S1
#df5<-cbind(dfT1[132:179],dfT1[251:257])
#For S1: Genera to remove to cleanup plot
df5$X.Ruminococcus..torques.group<-NULL
df5$X.Eubacterium..coprostanoligenes.group<-NULL
df5$X.Eubacterium..eligens.group<-NULL
df5$Ruminococcaceae.UCG.014<-NULL
df5$Ruminococcaceae_uc<-NULL
df5$Parabacteroides<-NULL
df5$Butyricicoccus<-NULL
df5$Ruminococcaceae.NK4A214.group<-NULL
df5$Ruminococcaceae.UCG.003<-NULL

#S2
df5<-cbind(dfT1[129:178],dfT1[255:261])
df5$Desulfovibrionaceae_uc<-NULL
df5$Slackia<-NULL
df5$Lachnospiraceae.FCS020.group<-NULL
df5$Clostridiales.vadinBB60.<-NULL
df5$Ruminococcaceae.UCG.003<-NULL
df5$X.Ruminococcus..torques.group<-NULL
df5$X.Eubacterium..hallii.group<-NULL
df5$Fusicatenibacter<-NULL
df5$Ruminococcaceae.UCG.013<-NULL
df5$Ruminococcaceae.NK4A214.group<-NULL





df5.c<-cor(df5,method="spearman",use = "complete.obs")
corrplot(df5.c,method="color",type="lower",order="hclust",mar = c(0, 0, 0, 0),tl.col="black",tl.srt = 45,col=colorRampPalette(c("firebrick4", "firebrick3","firebrick1","white","royalblue2","royalblue3","royalblue4"))(200))









#Compile multiple files:
df1<-read.csv("SAP_S1_FamGenASVAlpha_RA_filtered3.csv")
df2<-read.csv("SAP_S2_FamGenASVAlpha_RA_filtered3.csv")
alelles<-read.csv("Genes/MEL6A_alleles_64_exon1.csv")
alelles1<-alelles[-c(2:6)]
alelles1<-alelles1[-c(328:849),]
alelles1<-subset(alelles1,CW1_Haplotype!="Het")
alelles1$Line<-sub("PI","PI_",alelles1$Line)

df1G<-df1[130:177]
df1G$Line<-df1$PI
df1GA<-merge(alelles1,df1G)
write.table(df1GA[-c(1)],"S1_genus4Lefse.txt",row.names=FALSE,sep="\t", quote = FALSE)

df1G_mean <- df1G %>%
  group_by(Line) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
df1GA<-merge(alelles1,df1G_mean)
write.table(df1GA[-c(1)],"S1_AVEgenus4Lefse.txt",row.names=FALSE,sep="\t", quote = FALSE)

df2G<-df2[126:175]
df2G$Line<-df2$PI
df2GA<-merge(alelles1,df2G)
write.table(df2GA[-c(1)],"S2_genus4Lefse.txt",row.names=FALSE,sep="\t", quote = FALSE)

df2G_mean <- df2G %>%
  group_by(Line) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
df2GA<-merge(alelles1,df2G_mean)
write.table(df2GA[-c(1)],"S2_AVEgenus4Lefse.txt",row.names=FALSE,sep="\t", quote = FALSE)









library(ggplot2)
library(tidyverse)
library(ggpubr)
library(sommer)
library(RColorBrewer)
library(gplots)
#library(PERMANOVA)

setwd("~/SAP/Validation/qPCR")

df<-read.csv("FaecaliResults2.csv")
#df<-read.csv("FaecaliResultsV5-6.csv")
#df<-read.csv("FaecaliResultsTEST.csv")
#df<-df[1:192,]

dfmeta<-df[1:9]
#dfmeta<-cbind(df[2],df[3:8])

df1 <- df %>%
  group_by(Sample) %>%
  summarise(CT_F = mean(CTValueF,na.remove=TRUE)) %>%
  dplyr::select(Sample,CT_F)

df2<- merge(dfmeta,df1)
df2$Well<-NULL

df3<-unique(df2)

#Faecali:
df3$LogCFU<-(df3$CT_F-33.308)/-3.3608
#Ecoli
#df3$LogCFU<-(df3$CT_F*-0.2916)+9.198
X<-df3$LogCFU
X[X<0]<-runif(1,0.01,0.036)
df3$LogCFU<-X
df3$Subject<-factor(df3$Subject,levels=c("S1","S2", "S3","S4","S5","S6","S7","S8","S9","S10","S11","S12"))
df3$Tannin<-factor(df3$Tannin,levels=c("T","NT"))
df3$Group2<-factor(df3$Group2,levels=c("Major_T","Major_NT","Minor_T","Minor_NT"))

dfS<-subset(df3,Treatment!="NA")
dfS<-subset(dfS,Treatment=="Sorghum")

dfMC<-subset(df3,Treatment!="NA")
dfMC<-subset(dfMC,Treatment!="Sorghum")
c<-list(c("A1","R"),c("A","A1"),c("A","A2"),c("A1","A2"))

c1<-list(c("T","NT"))
c2<-list(c("Major","Minor"))
c3<-list(c("Major_T","Minor_T"),c("Major_NT","Minor_NT"))

a<-ggplot(dfS,aes(x=Group2,y=LogCFU,fill=Tannin))+scale_fill_manual(values = c("gold","forestgreen"))+geom_boxplot()+facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
a+stat_compare_means()+ylab("Log10(CFU_Faecalibacterium)")

b<-ggplot(dfS,aes(x=Group2,y=LogCFU,fill=Pool))+scale_fill_manual(values = c("black","blue"))+geom_boxplot()+facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
b+stat_compare_means()+ylab("Log10(CFU_Faecalibacterium)")

c<-ggplot(dfS,aes(x=Allele,y=LogCFU,fill=Tannin))+scale_fill_manual(values = c("brown","blue"))+geom_boxplot()+
  facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=18,family = "sans"))
c+stat_compare_means(comparisons = c2)+ylab("Log10(CFU_Faecalibacterium)")

d<-ggplot(dfS,aes(x=Allele,y=LogCFU,fill=Tannin))+scale_fill_manual(values = c("brown","blue"))+geom_boxplot()+
  facet_wrap(~Subject)+theme_classic()+theme(text=element_text(size=18,family = "sans"))
d+stat_compare_means()+ylab("Log10(CFU_Faecalibacterium)")

####################
a<-ggplot(dfS,aes(x=Group2,y=LogCFU))+geom_point()+facet_wrap(~Subject)+ theme_classic()+
  theme(text=element_text(size=20,family = "sans"))+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
a+stat_compare_means(comparisons = c3)+ylab("Log10(CFU_Faecalibacterium)")

a<-ggplot(dfS,aes(x=Group2,y=LogCFU,color=Subject))+geom_point()+ theme_classic()+
  theme(text=element_text(size=20,family = "sans"))+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
a+stat_compare_means(comparisons = c3)+ylab("Log10(CFU_Faecalibacterium)")
######################

dfS_T<-subset(dfS,Tannin=="T")
dfS_NT<-subset(dfS,Tannin=="NT")

d<-ggplot(dfS_T,aes(x=Group2,y=LogCFU,fill=Pool))+scale_fill_manual(values = c("red","purple"))+geom_boxplot()+
  facet_wrap(~Subject)+theme_classic()+theme(text=element_text(size=18,family = "sans"))
d+stat_compare_means(comparisons = c2)+ylab("Log10(CFU_Faecalibacterium)")

a<-ggplot(dfS,aes(x=Group1,y=LogCFU,fill=Tannin))+scale_fill_manual(values = c("brown","blue"))+
  geom_boxplot()+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
a+stat_compare_means(comparisons = c2)+ylab("Log10(CFU_Faecalibacterium)")


fit<-mmer(LogCFU~1,random=~Subject+Tannin+Allele+Tannin:Allele,rcov=~units,data = dfS)
summary(fit)

fit<-lm(LogCFU~Subject+Tannin+Allele+Tannin:Allele+Subject:Tannin+Subject:Allele,data=dfS)
summary(fit)

fit0<-lm(LogCFU~Subject+Allele,data=dfS)
L0<-logLik(fit0)
fit1<-lm(LogCFU~Subject+Tannin+Allele,data=dfS)
L1<-logLik(fit1)
Tan<-1-pchisq(2*(L1-L0),df=14)

fit0<-lm(LogCFU~Subject+Tannin,data=dfS)
L0<-logLik(fit0)
fit1<-lm(LogCFU~Subject+Tannin+Allele,data=dfS)
L1<-logLik(fit1)
All<-1-pchisq(2*(L1-L0),df=14)

#PERMANOVA
effects<-factor(c(1,2,3,4))
levels(effects)<-c("Subject","Tannin","Allele","Interaction")
X<-dfS$LogCFU
X<-IniTransform(X)
D = DistContinuous(X)

per1<-PERMANOVA(D,dfS$Sample,Effects = Effects)

All<-aov(LogCFU~Subject+Tannin+Allele+Subject:Tannin+Subject:Allele+Tannin:Allele,data=dfS)
summary(All)
dfS1<-subset(dfS,Subject=="S2")
All<-aov(LogCFU~Tannin+Allele+Tannin:Allele,data=dfS1)
summary(All)

All<-aov(LogCFU~Tannin+Allele+Tannin:Allele,data=dfS)
summary(All)
##################
dfMC$Treatment<-factor(dfMC$Treatment,levels=c("Media","Tannin","FOS","Tan+FOS"))

a<-ggplot(dfMC,aes(x=Treatment,y=LogCFU))+scale_fill_manual(values = c("gold","forestgreen"))+
  geom_boxplot()+facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
a+stat_compare_means()+ylab("Log10(CFU_Ecoli)")













#df4<-subset(df3,Type!="K55")


a<-
b<-ggplot(df4,aes(x=Type,y=LogEcoliCFU,fill=Type))+scale_fill_manual(values = c("gold","forestgreen"))+geom_boxplot()+facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
b+stat_compare_means()+ylab("Log10(CFU_Ecoli)")

a<-ggplot(df3,aes(x=Type,y=LogRoseCFU,fill=Type))+scale_fill_manual(values = c("gold","brown","forestgreen"))+geom_point()+facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
a+stat_compare_means()+ylab("Log10(CFU_Roseburia)")

b<-ggplot(df3,aes(x=Type,y=LogEcoliCFU,fill=Type))+scale_fill_manual(values = c("gold","brown","forestgreen"))+geom_boxplot()+facet_wrap(~Subject)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))
b+stat_compare_means()+ylab("Log10(CFU_Ecoli)")



###### #log fold change heatmap:
dfS<-subset(dfS,LogCFU!="NA")

dfS_T<-subset(dfS,Tannin=="T")
dfS_NT<-subset(dfS,Tannin=="NT")
dfS_NT<-subset(dfS_NT,Sample!="S1_Minor_NTB3")
dfS_NT<-subset(dfS_NT,Sample!="S4_Major_NTA3")

dfS_T1 <- dfS_NT %>%
  group_by(Group2) %>%
  summarise(Log_A = mean(LogCFU,na.remove=TRUE)) %>%
  dplyr::select(Group2,Log_A)
dfS_T2<-as.data.frame(t(dfS_T1))


calculate_log_fold_change <- function(subject,df) {
  S1 <- subset(df, Subject == subject)
  dfS1 <- S1 %>%
    group_by(Group2) %>%
    summarise(CT_A = mean(CT_F,na.remove=TRUE)) %>%
    dplyr::select(Group2,CT_A)
  dfS2<-as.data.frame(t(dfS1))
  names(dfS2)<-c("Major","Minor")
  dfS2 <- dfS2[-c(1), ]
  dfS2 <- transform(dfS2, Major = as.numeric(Major), Minor = as.numeric(Minor))
  dfS2$LogFoldChange <- log10((dfS2$Minor + runif(1,.00009,.0001)) / (dfS2$Major + runif(1,.00009,.0001)))
  S_Log<-as.data.frame(dfS2$LogFoldChange)
  names(S_Log)<-paste0(subject)
  rownames(S_Log)<-rownames(dfS2)
  return(S_Log)
}

X<-calculate_log_fold_change("S1",dfS_T)
Subjects<-unique(dfS_T$Subject)
for (i in Subjects){ 
  X<-cbind(X,calculate_log_fold_change(i,dfS_T))
}
X2<-calculate_log_fold_change("S1",dfS_NT)
for (i in Subjects){ 
  X2<-cbind(X2,calculate_log_fold_change(i,dfS_NT))
}

Values<-rbind(X,X2)
rownames(Values)<-c("Tannin","Non-Tannin")
Values<-Values[-c(1)]

ValuesM<-as.matrix(Values)
max1<-round(max(ValuesM),4)
min1<-round(min(ValuesM),4)
med1<-max1+min1

col.order <- c("S6","S1", "S9","S5","S11","S7","S3","S8","S4","S12","S2","S10")
ValuesM<-ValuesM[ , col.order]

#svg("OverMajorFaecali.svg",width = 14,height=8)
heatmap(ValuesM,scale="none", Rowv = NA,col= colorRampPalette(brewer.pal(8, "Reds"))(25))
legend(x="bottomright", legend=c(min1, med1, max1), 
       fill=colorRampPalette(brewer.pal(8, "Reds"))(3))
#dev.off()

pal<-colorpanel(20,"forestgreen","white","purple")
pdf("MajorOverMinorFaecaliqPCR.pdf",width = 11,height=8)
heatmap.2(ValuesM,scale="none", Rowv = NA, Colv = NA, col = pal,trace = "none",margins = c(3,17))
dev.off()




#ecoliplot
#svg("MinorOverMajorecoli.svg",width = 14,height=8)
heatmap(ValuesM,scale="none", Rowv = NA,col= rev(colorRampPalette(brewer.pal(8, "Blues"))(25)))
legend(x="bottomright", legend=c(min1, med1, max1), 
       fill=rev(colorRampPalette(brewer.pal(8, "Blues"))(3)))
#dev.off()
d<-ggplot(dfS_NT,aes(x=Group2,y=LogCFU))+scale_fill_manual(values = c("red","purple"))+geom_boxplot()+
  facet_wrap(~Subject)+theme_classic()+theme(text=element_text(size=18,family = "sans"))
d+stat_compare_means()+ylab("Log10(CFU_Faecalibacterium)")

##
#Plot Faecali and E coli qPCR

dfF <- df %>%
  group_by(Sample) %>%
  summarise(CT_F = mean(CTValueF,na.remove=TRUE)) %>%
  dplyr::select(Sample,CT_F)
dfE <- df %>%
  group_by(Sample) %>%
  summarise(CT_E = mean(CTValueE,na.remove=TRUE)) %>%
  dplyr::select(Sample,CT_E)

df2<- merge(dfmeta,dfF)
df2<- merge(df2,dfE)
df2$Well<-NULL

df3<-unique(df2)

#Faecali:
df3$LogCFU_F<-(df3$CT_F-33.308)/-3.3608
#Ecoli
df3$LogCFU_E<-(df3$CT_E*-0.2916)+9.198

df4<-subset(df3,Treatment!="Sorghum")
df4<-subset(df4,Treatment!="NA")
df4<-subset(df4,Treatment!="Tannin")
df4<-subset(df4,Treatment!="Tan+FOS")
df4$Treatment<-factor(df4$Treatment,levels=c("Media","FOS"))

a<-ggplot(df4,aes(x=Treatment,y=LogCFU_F))+geom_boxplot()+facet_wrap(~Subject)+theme_classic()+theme(text=element_text(size=18,family = "sans"))+
  stat_compare_means()+ylab("Log10(CFU_Faecalibacterium)")+xlab(NULL)

b<-ggplot(df4,aes(x=Treatment,y=LogCFU_E))+geom_boxplot()+facet_wrap(~Subject)+theme_classic()+theme(text=element_text(size=18,family = "sans"))+
  stat_compare_means()+ylab("Log10(CFU_Ecoli)")+xlab(NULL)

ggarrange(a,b,labels =c("A","B"))

setwd("~/SAP/GWAS/genes")
library(dplyr)

Genes<-read.csv("SorghumGenes.csv")

Genes<-subset(Genes,Type=="gene")
Genes<-Genes[1:6]

GenesM1<-subset(Genes,Chr==1)
GenesM1<-subset(GenesM1,Stop>=54395657)
GenesM1<-subset(GenesM1,Start<=55681861)

GenesM1$MEL<-"MEL1"
GeneList<-GenesM1

GenesM2<-subset(Genes,Chr==2)
GenesM2<-subset(GenesM2,Stop>=6691792)
GenesM2<-subset(GenesM2,Start<=6785367)

GenesM2$MEL<-"MEL2"
GeneList<-rbind(GeneList,GenesM2)

GenesM2<-subset(Genes,Chr==2)
GenesM21<-subset(GenesM2,Stop>=(9668945-50000))
GenesM21<-subset(GenesM21,Start<=9668869+50000)

GenesM21$MEL<-"MEL2.1"
GeneList<-rbind(GeneList,GenesM21)

GenesM3<-subset(Genes,Chr==3)
GenesM3A<-subset(GenesM3,Stop>=7453735)
GenesM3A<-subset(GenesM3A,Start<=8984424)

GenesM3A$MEL<-"MEL3A"
GeneList<-rbind(GeneList,GenesM3A)

GenesM3B<-subset(GenesM3,Stop>=52867187-50000)
GenesM3B<-subset(GenesM3B,Start<=52868173+50000)

GenesM3B$MEL<-"MEL3B"
GeneList<-rbind(GeneList,GenesM3B)

GenesM3B1<-subset(GenesM3,Stop>=53552257)
GenesM3B1<-subset(GenesM3B1,Start<=55499313)

GenesM3B1$MEL<-"MEL3B.1"
GeneList<-rbind(GeneList,GenesM3B1)

GenesM4<-subset(Genes,Chr==4)
GenesM4A<-subset(GenesM4,Stop>=5883328-50000)
GenesM4A<-subset(GenesM4A,Start<=5893095+50000)

GenesM4A$MEL<-"MEL4A"
GeneList<-rbind(GeneList,GenesM4A)

GenesM4A1<-subset(GenesM4,Stop>=7197553-50000)
GenesM4A1<-subset(GenesM4A1,Start<=7246890+50000)

GenesM4A1$MEL<-"MEL4A.1"
GeneList<-rbind(GeneList,GenesM4A1)

GenesM4B<-subset(GenesM4,Stop>=61888295)
GenesM4B<-subset(GenesM4B,Start<=63000085)

GenesM4B$MEL<-"MEL4B"
GeneList<-rbind(GeneList,GenesM4B)

GenesM5<-subset(Genes,Chr==5)
GenesM5A<-subset(GenesM5,Stop>=7210037-50000)
GenesM5A<-subset(GenesM5A,Start<=7221872+50000)

GenesM5A$MEL<-"MEL5A"
GeneList<-rbind(GeneList,GenesM5A)

GenesM5A1<-subset(GenesM5,Stop>=12139539)
GenesM5A1<-subset(GenesM5A1,Start<=12356416)

GenesM5A1$MEL<-"MEL5A.1"
GeneList<-rbind(GeneList,GenesM5A1)

GenesM5B<-subset(GenesM5,Stop>=51418161-50000)
GenesM5B<-subset(GenesM5B,Start<=51419786+50000)

GenesM5B$MEL<-"MEL5B"
GeneList<-rbind(GeneList,GenesM5B)

GenesM6<-subset(Genes,Chr==6)
GenesM6A<-subset(GenesM6,Stop>=42667079)
GenesM6A<-subset(GenesM6A,Start<=43971287)

GenesM6A$MEL<-"MEL6A"
GeneList<-rbind(GeneList,GenesM6A)

GenesM6B<-subset(GenesM6,Stop>=51620829)
GenesM6B<-subset(GenesM6B,Start<=57307490)

GenesM6B$MEL<-"MEL6B"
GeneList<-rbind(GeneList,GenesM6B)

GenesM6C<-subset(GenesM6,Stop>=59335835)
GenesM6C<-subset(GenesM6C,Start<=60149021)

GenesM6C$MEL<-"MEL6C"
GeneList<-rbind(GeneList,GenesM6C)

GenesM7<-subset(Genes,Chr==7)
GenesM7<-subset(GenesM7,Stop>=56567874)
GenesM7<-subset(GenesM7,Start<=56989886)

GenesM7$MEL<-"MEL7"
GeneList<-rbind(GeneList,GenesM7)

GenesM80<-subset(Genes,Chr==8)
GenesM8<-subset(GenesM80,Stop>=4311803)
GenesM8<-subset(GenesM8,Start<=4463398)

GenesM8$MEL<-"MEL8"
GeneList<-rbind(GeneList,GenesM8)

GenesM81<-subset(GenesM80,Stop>=5732816-50000)
GenesM81<-subset(GenesM81,Start<=5741529+50000)

GenesM81$MEL<-"MEL8.1"
GeneList<-rbind(GeneList,GenesM81)

GenesM10<-subset(Genes,Chr==10)
GenesM10A<-subset(GenesM10,Stop>=11068459)
GenesM10A<-subset(GenesM10A,Start<=14614059)

GenesM10A$MEL<-"MEL10A"
GeneList<-rbind(GeneList,GenesM10A)

GenesM10A1<-subset(GenesM10,Stop>=15201036-50000)
GenesM10A1<-subset(GenesM10A1,Start<=15201036+50000)

GenesM10A1$MEL<-"MEL10A.1"
GeneList<-rbind(GeneList,GenesM10A1)

GenesM10B<-subset(GenesM10,Stop>=48735787-50000)
GenesM10B<-subset(GenesM10B,Start<=48798761+50000)

GenesM10B$MEL<-"MEL10B"
GeneList<-rbind(GeneList,GenesM10B)

GenesM10B1<-subset(GenesM10,Stop>=51339917-50000)
GenesM10B1<-subset(GenesM10B1,Start<=51343698+50000)

GenesM10B1$MEL<-"MEL10B.1"
GeneList<-rbind(GeneList,GenesM10B1)
GeneListIDs<-as.data.frame(GeneList$ID)
names(GeneListIDs)<-"PI"
#write.csv(GeneListIDs,"GenesInMEL.csv")

#Transcriptomics:
RNAseq<-read.csv("ALLGenes_Expression.csv")
Names<-read.table("Sbicolor_454_v3.1.1.synonym.txt")
names(Names)<-c("PI","Old")
Names<-Names[-c(1),]

Names$PI<-gsub("\\.\\d+$", "", Names$PI)


NamesU<-as.data.frame(unique(Names))
RNAseqN<-merge(RNAseq,NamesU,by.x = "gene.ID",by.y="Old")

GeneListTranscripts<-merge(GeneList,RNAseqN,by.x = "ID",by.y = "PI")

GeneListTranscripts$Seed<-pmax(GeneListTranscripts$Seed_5d_After_Pollination,GeneListTranscripts$Seed_10d_After_Pollination,GeneListTranscripts$Embryo_25d_After_Pollination,GeneListTranscripts$Endosperm_25d_After_Pollination)

GeneListTranscripts2<-subset(GeneListTranscripts,Seed>=1)
#GeneListTranscripts2$Seed<-NULL

GeneListTranscripts3<-subset(GeneListTranscripts2,MEL=="MEL10B")

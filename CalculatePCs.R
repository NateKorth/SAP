##Library packages (install.packages() if you need to)
library(sommer)
library(ggplot2)
library(readxl)
library(data.table)
library(dplyr)
library(tidyverse)
library(bestNormalize)
library(compositions)

#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/SAP/GWAS")

df1<-read.csv("SAP_S2_FamGenASVAlpha_RA_filtered2.csv")

df1<-subset(df1,PI!="PI_656081")

#make factors characters for random effx
df1$Batch<-as.character(df1$Batch)
df1$Plate<-as.character(df1$Plate)
df1$Column<-as.character(df1$Column)

Traits<-as.data.frame(names(df1[100:258]))
names(Traits)<-"Traits"

df2<-df1[177:252]
dfx<-as.data.frame(cbind(df1$PI,df1$Batch,df1$Plate))
names(dfx)<-c("PI","Batch","Plate")

df3<-as.data.frame(clr(df2))

df4<-cbind(dfx,df3)

a<-lm(formula = cbind(ASV1_Phascolarctobacterium, ASV2_Faecalibacterium)  ~PI,data=df4)

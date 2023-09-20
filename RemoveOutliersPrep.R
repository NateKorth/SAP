setwd("~/SAP/GWAS/")
library(sommer)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(bestNormalize)
library(ggpubr)

dfS1<-read.csv("SAP_S1_FamGenASVAlpha_RA_filtered2.csv")
TraitsS1<-as.data.frame(names(dfS1[101:254]))
names(TraitsS1)<-"Traits"

dfS2<-read.csv("SAP_S2_FamGenASVAlpha_RA_filtered2.csv")
TraitsS2<-as.data.frame(names(dfS2[100:258]))
names(TraitsS2)<-"Traits"

S1meta<-read.csv("../PlateLayouts_SampleSheets/SAP_S1_Metadata2.csv")
S2meta<-read.csv("../PlateLayouts_SampleSheets/SAP_S2_Metadata2.csv")

biochem<-read.csv("./phenoinput/SAP_BiochemicalPhenotypes.csv")
biochem2<-as.data.frame(cbind(biochem$Line,biochem$Tannins_NK))
names(biochem2)<-c("PI","Tannins_NK")

#remove outliers(SAME PROTOCOL AS BLUE PROCESSING)
rmOut<-function(OTU){
  X<-dfS1[[paste0(OTU)]]
  mean<-mean(X)
  std<-sd(X)
  X[X<(mean-5*std)]<-NA
  X[X>(mean+5*std)]<-NA
  return(X)
}

X<-dfS1$Shannon
mean<-mean(X)
std<-sd(X)
X[X<(mean-5*std)]<-NA
X[X>(mean+5*std)]<-NA
dfS1.1<-as.data.frame(X)
for (i in TraitsS1$Traits){ 
  dfS1.1<-cbind(dfS1.1,rmOut(i))
}
sum(is.na(dfS1.1))
dfS1.1$X<-NULL
names(dfS1.1)<-TraitsS1$Traits

rmOut<-function(OTU){
  X<-dfS2[[paste0(OTU)]]
  mean<-mean(X)
  std<-sd(X)
  X[X<(mean-5*std)]<-NA
  X[X>(mean+5*std)]<-NA
  return(X)
}
X<-dfS2$Shannon
mean<-mean(X)
std<-sd(X)
X[X<(mean-5*std)]<-NA
X[X>(mean+5*std)]<-NA
dfS2.1<-as.data.frame(X)
for (i in TraitsS2$Traits){ 
  dfS2.1<-cbind(dfS2.1,rmOut(i))
}
sum(is.na(dfS2.1))
dfS2.1$X<-NULL
names(dfS2.1)<-TraitsS2$Traits

#function2:
replace_out <- function(column) {
  qnt <- quantile(column, probs=c(.25, .75),na.rm=TRUE)
  upper_whisker <- 2.5 * IQR(column,na.rm=TRUE)
  clean_data <- column
  clean_data[column < (qnt[1] - upper_whisker)] <- NA
  clean_data[column > (qnt[2] + upper_whisker)] <- NA
  clean_data
}

#Remove outliers within lines:
dfS1.2<-as.data.frame(cbind(dfS1$PI,dfS1.1))
names(dfS1.2)<-c("PI",names(dfS1.1))
dfS1.3 <- dfS1.2 %>% group_by(PI) %>% mutate_if(is.numeric, replace_out)
sum(is.na(dfS1[101:254]))
sum(is.na(dfS1.2))
sum(is.na(dfS1.3))

dfS2.2<-as.data.frame(cbind(dfS2$PI,dfS2.1))
names(dfS2.2)<-c("PI",names(dfS2.1))
dfS2.3 <- dfS2.2 %>% group_by(PI) %>% mutate_if(is.numeric, replace_out)
sum(is.na(dfS2[100:258]))
sum(is.na(dfS2.2))
sum(is.na(dfS2.3))


dfS1.4<-cbind(dfS1[1],dfS1.3)
dfS2.4<-cbind(dfS2[2],dfS2.3)

##########################################################################

#Remove non-heritable taxa
h2S1List<-read.csv("S1_SAP_h2_estimations.csv")
h2S2List<-read.csv("S2_SAP_h2_estimations.csv")

h2S1List2<-as.data.frame(h2S1List[which(h2S1List$h2<=0.05  ),])
h2S2List2<-as.data.frame(h2S2List[which(h2S2List$h2<=0.05  ),])

dfS1.5<-as.data.frame(dfS1.4[,!(names(dfS1.4) %in% h2S1List2$Trait)])
dfS2.5<-as.data.frame(dfS2.4[,!(names(dfS2.4) %in% h2S2List2$Trait)])  


#########################################################################
#merge asvtable, metadata, and tannin measurments
dfS1.5$PI<-NULL
dfS2.5$PI<-NULL
dfS1.6<-merge(dfS1.5,S1meta,by.x="Name",by.y="ID_S1")
dfS1.6<-merge(dfS1.6,biochem2)
dfS2.6<-merge(dfS2.5,S2meta,by.x="Name",by.y="ID")
dfS2.6<-merge(dfS2.6,biochem2)

dfS1.6<-subset(dfS1.6,PI!="PI_656081")
dfS2.6<-subset(dfS2.6,PI!="PI_656081")

#Write file:
#write.csv(dfS1.6,"./AutoEncoder/MicrobiomeTraits_h2trim_S1.csv")
#write.csv(dfS2.6,"./AutoEncoder/MicrobiomeTraits_h2trim_S2.csv")

dfS1.6$Red<-as.character(dfS1.6$Red)
dfS1.6C<-subset(dfS1.6,Red!="NA")
dfS1.6CNT<-subset(dfS1.6C,Tannins_NK<=1)

dfS2.6$Red<-as.character(dfS2.6$Red)
dfS2.6C<-subset(dfS2.6,Red!="NA")
dfS2.6CNT<-subset(dfS2.6C,Tannins_NK<=1)
dfS2.6C$Tannins_D<-as.character(dfS2.6C$Tannins_D)

#Plot color by traits:
a<-ggplot(dfS1.6CNT,aes(x=Red,y=Isovaleric))+geom_dotplot(binaxis='y', stackdir='center',stackratio=.8,dotsize = 0.5, position=position_dodge(0.8))
a+stat_compare_means()


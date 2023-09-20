library(ggplot2)
library(bestNormalize)
setwd("~/SAP/GWAS")

df<-read.csv("./phenoinput/SAP_BiochemicalPhenotypes2.csv")
TanninsN_NK<-bestNormalize(df$Tannins_NK)
df$TanninsN_NK<-TanninsN_NK$x.t

hist(df$Tannins_E)
hist(df$TanninsN_NK)

df<-read.csv("./phenoinput/S2_SAP_PhenotypesF.csv")
hist(df$ProteinPCT)
hist(df$OilPCT)
hist(df$FiberPCT)
hist(df$AshPCT)
hist(df$StarchPCT)
#write.csv(df,"./phenoinput/SAP_BiochemicalPhenotypes.csv",quote = FALSE,row.names=FALSE)


hist(df2$Polyphenols_E)

df<-read.csv("./phenoinput/SAP_SeedPhenotypes.csv")
df2<-data.frame(df[1])

Traits<-as.data.frame(names(df[-c(1)]))
names(Traits)<-"Traits"

normal1 <- function(trait){
  N1<-bestNormalize(df[[trait]])
  return(N1$x.t)
  gc()
}

for (i in Traits$Traits){ 
  df2<-cbind(df2,normal1(i))
}

names(df2)<-names(df)

write.csv(df2,"./phenoinput/SAP_SeedPhenotypesN.csv",quote = FALSE,row.names=FALSE)

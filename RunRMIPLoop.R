library(rMVP)
phenotype1<-read.csv("../../input/SAP_CompiledBLUEs3.csv",head=TRUE)
genotype<-attach.big.matrix("/scratch/SAP.geno.desc")
map<-read.table("/scratch/SAP.geno.map", head = TRUE)
kinship<-attach.big.matrix("/scratch/SAP.kin.desc")
covariates_PC<-bigmemory::as.matrix(attach.big.matrix("/scratch/SAP.pc.desc"))

#remove columns not to be ran (i.e. the line info):
phenotype2 <- as.data.frame(phenotype1[-c(1)])
phenotype2 <- as.data.frame(phenotype2[c(367:370)])

#S1:
#phenotype2 <- as.data.frame(phenotype2[c(1:107)])

#S2
#phenotype2 <- as.data.frame(phenotype2[c(108:214)])

#Biochemical
#phenotype2 <- as.data.frame(phenotype2[c(215:246)])

#Seed
#phenotype2 <- as.data.frame(phenotype2[c(247:303)])

#PolyS1
#phenotype2 <- as.data.frame(phenotype2[c(304:327)])

#PolyS2
#phenotype2 <- as.data.frame(phenotype2[c(328:351)]

#Poly
#phenotype2 <- as.data.frame(phenotype2[c(352:355)])

#Agronomic
#phenotype2 <- as.data.frame(phenotype2[c(356:444)])


phenolist<-names(phenotype2)

#Function to perform 100 iterations of FarmCPU removing ~10% of the data to change amount alter z
thresh<-0.05/861521.36

RunRMIP<-function(x,column){
  ph<-as.data.frame(cbind(phenotype1$Line,x[,column]))
  
  set.seed(40)
  for (i in 3:103){
    z <- sample(1:328,32)
    ph[,i] <- ph[,2]
    ph[z,i] <- NA
    rm(z)
  }
  
  RMIP <- c()
  for(j in 3:ncol(ph)){
    imMVP <- MVP(
      phe=ph[, c(1, j)],
      geno=genotype,
      map=map,
      K=kinship,
      CV.FarmCPU=covariates_PC,
      priority="speed",
      maxLoop=10,
      ncpus=225,
      method.bin="FaST-LMM",
      method=c("FarmCPU"),
      file.output = F,
      p.threshold = (thresh))
    farm <- cbind(imMVP$map, imMVP$farmcpu.results)
    farm <- na.omit(farm[farm[,8]<thresh,])
    colnames(farm)[8] <- "pvalue"
    RMIP <- rbind(RMIP, farm)
    rm(farm, imMVP)
    gc()
  }
  table(RMIP$SNP)
  
  SAPSNPs<-table(RMIP$SNP)
  return(SAPSNPs)
}

for (i in 1:length(phenolist)){
  RMIPSNPs<-RunRMIP(phenotype2,i)
  if(nrow(RMIPSNPs)>1){
    write.csv(RMIPSNPs,paste0(phenolist[i],"_RMIP.csv"),row.names=FALSE)
  }
}

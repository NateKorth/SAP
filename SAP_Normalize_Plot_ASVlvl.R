library(metagenomeSeq)
library(phyloseq)
library(ape)
library(ggplot2)
library(reshape2)
library(vegan)
library(superheat)
#library(DESeq2)
library(ggfortify)
library(corrplot)
library(stringr)
library(sommer)
library(data.table)
library(ggpubr)
library(devtools)
library(tidyverse)
##Prepare data for correlation analysis at ASV level, remove FBBs!

setwd("~/SAP")

#tree=read_tree("exported_S1/SAP_S1_tree.nwk")
#tree=read_tree("exported_S2/SAP_S2_tree.nwk")
table.alpha<-read.delim("./exported_S1/SAP_S1_OTU_Table_silva_taxonomy.tsv")
table.alpha<-read.delim("./exported_S2/SAP_S2_OTU_Table_silva_taxonomy.tsv")
table.alpha[1:4,1:4]

#remove ID column and add it as column names
A.OTU<-table.alpha
OTU_IDs<-as.data.frame(table.alpha[1])
A.OTU<-table.alpha[-c(1)]
rownames(A.OTU)<-OTU_IDs[,1]

A.OTU[1:3,1:3]
A.OTU$Consensus.Lineage<-NULL

#trim low abundance taxa:
B.OTU<-A.OTU
B.OTU$Undetermined<-NULL
B.OTU$AveReads<-rowMeans(B.OTU)
B.OTU<-B.OTU[-which(B.OTU$AveReads<=.01),]
B.OTU$AveReads<-NULL

###################################################################
#read in Silva taxonomy | before, in excel seperate taxonomy into columns
taxa<-read.delim("./exported_S1/SAP_S1_silva_taxonomy_biom.tsv",stringsAsFactors = FALSE,sep="\t")
taxa<-read.delim("./exported_S2/SAP_S2_silva_taxonomy.tsv",stringsAsFactors = FALSE,sep="\t")
#taxa<-taxa[c(1:9)]
ASV.num<-as.data.frame(paste("ASV",seq(1,nrow(taxa), by=1),sep=""))
names(ASV.num)<-"ASV.num"
ASV.name<-as.data.frame(taxa$Genus)
names(ASV.name)<-"ASV.name"
taxa$ASV.IDs<-paste0(ASV.num$ASV.num,"_",ASV.name$ASV.name,sep="")
taxa$ASV.IDs

#taxaonomy<-read.delim("NAME")
#rownames(df)<-paste(seq(1,177),taxaonomy$GEUNS)
#df1<-df[1:5] 

taxa1<-taxa[which(taxa$OTUID %in% row.names(B.OTU)),]

ordB<- match(rownames(B.OTU), taxa1$OTUID)
taxa1 <- taxa1[ordB,]

OTUdata<-AnnotatedDataFrame(taxa1)
rownames(OTUdata)<-taxa1$OTUID

#normalize data
#load metadata
MetaS1<-load_phenoData("PlateLayouts_SampleSheets/SAP_S1_Metadata2.txt",sep="\t")
MetaS2<-load_phenoData("PlateLayouts_SampleSheets/SAP_S2_Metadata2.txt",sep="\t")
#name the columns in Meta data file
#colnames(MetaS2)<-colnames(MetaS1)

Meta<-MetaS1
Meta<-MetaS2

#Match names in meta data to OTU table
ordA<- match(colnames(B.OTU), rownames(Meta))
Ameta <- Meta[ordA,]
rownames(Ameta)<-names(B.OTU)
#insert subset into phenodataframe:
phenotypeDataA<-AnnotatedDataFrame(Ameta)

#Make MetagenomeSeq object
A.data<-newMRexperiment(B.OTU, phenoData = phenotypeDataA, featureData = OTUdata)

#trim that shit yo: taxa only present in 20% of taxa removed--Change this number based on how many samples you have!
A.datatrim<-filterData(A.data,present = 50)

#normalize by CSS:
pA<-cumNormStatFast(A.datatrim)
A.datatrim<-cumNorm(A.datatrim,p=pA)

#log=TRUE -log2 transform
#SAP_normal_log<-as.data.frame(MRcounts(A.datatrim,norm=TRUE,log=FALSE))
SAP_normal_log<-as.data.frame(MRcounts(A.datatrim,norm=FALSE,log=FALSE))
#UFMU_A_normal_log<-log10(UFMU_A_normal_log)
SAP_normal_log[1:5,1:5]
#trimmed relative abundance OTU table:
SAP.relabun.tr<-as.data.frame(MRcounts(A.datatrim,norm=TRUE,log=FALSE))
SAP.relabun.trim<-sweep(SAP.relabun.tr,2,colSums(SAP.relabun.tr),"/")

#####################################################################################

######################################
#prep data for plotting:
#choose read counts or relative abundance
##Use Relative abundance for plotting, and abesolute data for diversity
P.RA<-SAP_normal_log
#P.RA<-SAP.relabun.trim
taxa2<-taxa1[which(taxa1$OTUID %in% row.names(P.RA)),]


ASVIDs<-str_replace_all(taxa2$ASV.IDs, "`","")
rownames(P.RA)<-ASVIDs
P.RA1<-as.data.frame(t(P.RA))
BMeta<-Ameta[which(row.names(Ameta) %in% row.names(P.RA1)),]

names(BMeta)[names(BMeta)=="V6"]<-"SAMPLENUM"
P.RA2<-cbind(P.RA1,BMeta$SAMPLENUM)

P.RA2.1<-subset(P.RA2, `BMeta$SAMPLENUM`=="FBB0")

#All data with FBBs
P.RA2.1<-subset(P.RA2, `BMeta$SAMPLENUM`!="FBB0")
P.RA2.1<-subset(P.RA2.1, `BMeta$SAMPLENUM`!="FBB16")
P.RA2.1<-subset(P.RA2.1, `BMeta$SAMPLENUM`!="PCRNC")
P.RA2.1<-subset(P.RA2.1, `BMeta$SAMPLENUM`!="DNANC")
P.RA2.1<-subset(P.RA2.1, `BMeta$SAMPLENUM`!="EMPTY")
P.RA2.1$`BMeta$SAMPLENUM`<-NULL
P.RA2.2<-as.data.frame(t(P.RA2.1))

#make phyloseq object:
#OTU table:P.RA
#taxa
taxa3<-taxa2[2:8]
rownames(taxa3)<-taxa2$ASV.IDs

OTU<-otu_table(P.RA2.2,taxa_are_rows = TRUE)
TAX<-tax_table(as.matrix(taxa3))
physeq<-phyloseq(OTU,TAX)

CMeta<-Ameta[which(row.names(Ameta) %in% names(P.RA2.2)),]
sampledata<-sample_data(CMeta)

physeq1<-merge_phyloseq(physeq,sampledata)

physeq1
#Even sampling Depth
total<-median(sample_sums(physeq1))
standf<-function(x,t=total) round(t*(x/sum(x)))
physeq2<-transform_sample_counts(physeq1,standf)


#Trim OTUs that show up more than 5 times in half the samples:
wh0 <- genefilter_sample(physeq2, filterfun_sample(function(x) x > 5), A=.75*nsamples(physeq2))
physeq3<-prune_taxa(wh0,physeq2)
#export tables with relative abundance
physeq3RA<-transform_sample_counts(physeq3, function(x) x / sum(x) )

#ASV level table

ASV1<-as(otu_table(physeq3RA),"matrix")
ASV2<-as.data.frame(t(ASV1))

#family level table
phyFam<-tax_glom(physeq2,taxrank="Family")
wh1 <- genefilter_sample(phyFam, filterfun_sample(function(x) x > 5), A=0.75*nsamples(phyFam))
phyFam<-prune_taxa(wh1,phyFam)
phyFamRA<-transform_sample_counts(phyFam, function(x) x / sum(x) )

taxa_names(phyFamRA)<-tax_table(phyFamRA)[,5]
Fam1<-as(otu_table(phyFamRA),"matrix")
Fam2<-as.data.frame(t(Fam1))

#genus level table
phygen<-tax_glom(physeq2,taxrank="Genus")
wh1 <- genefilter_sample(phygen, filterfun_sample(function(x) x > 5), A=0.75*nsamples(phygen))
phygen<-prune_taxa(wh1,phygen)
phygenRA<-transform_sample_counts(phygen, function(x) x / sum(x) )



taxa_names(phygenRA)<-tax_table(phygenRA)[,6]
gen1<-as(otu_table(phygenRA),"matrix")
gen2<-as.data.frame(t(gen1))

#gen3<-as.data.frame(gen1)
#gen3$S1<-(rowSums(gen3)/4)
#gen3<-gen3[5]
#write.csv(gen3,"S2_FBB0.csv")
#####################
###Calculate average abundance of every genus for biplot:
dfAv<-as.data.frame(t(gen2))
dfAv1<-as.data.frame(rowMeans(dfAv))
write.csv(dfAv1,"S2_GenusAverage")


All<-cbind(CMeta[7],CMeta[52],Fam2,gen2,ASV2)
#write.csv(All,"SAP_S1_FamGenASV_RelativeAbundance_Strictfiltered.csv")
#Strict=x>5,A=.5
#Less=x>3,A=.2
#Barely=x>1,A=.1



#Plot Beta Diversity
physeq3.ord<-ordinate(physeq2,"PCoA","jaccard")

p1<-plot_ordination(physeq2,physeq3.ord,type="samples",color="Pop")

p1<-plot_ordination(physeq3,physeq3.ord,type="samples",color="Tannins_D")

p1<-plot_ordination(physeq3,physeq3.ord,type="samples",color="PericarpColor_D")

p1+stat_ellipse()+ geom_point(size=4)

S1_PC_Jac<-physeq3.ord$vectors

S1_PC1_Jac<-as.data.frame(S1_PC_Jac[,1])
S1_PC1_Jac<-cbind(S1_PC1_Jac,sampledata$SAMPLENUM)
names(S1_PC1_Jac)<-c("axis1","Line")
S1_PC1_Jac.meds<-S1_PC1_Jac %>% group_by(Line) %>% summarise_each(funs(mean(.,na.rm = TRUE)))

S1_PC1_Jac2<-merge(S1_PC1_Jac,S1_PC1_Jac.meds,by.x="Line",by.y="Line")
names(S1_PC1_Jac2)<-c("Line","axis1","ave1")

S1_PC2_Jac<-as.data.frame(S1_PC_Jac[,2])
S1_PC2_Jac<-cbind(S1_PC2_Jac,sampledata$SAMPLENUM)
names(S1_PC2_Jac)<-c("axis2","Line")
S1_PC2_Jac.meds<-S1_PC2_Jac %>% group_by(Line) %>% summarise_each(funs(mean(.,na.rm = TRUE)))

S1_PC2_Jac2<-merge(S1_PC2_Jac,S1_PC2_Jac.meds,by.x="Line",by.y="Line")
names(S1_PC2_Jac2)<-c("Line","axis2","ave2")
rownames(S1_PC1_Jac2)<-rownames(S1_PC1_Jac)
S1_PCs_Jac<-cbind(S1_PC1_Jac2,S1_PC2_Jac2[,2:3])

S1_PCs_Jac$dist<-sqrt((S1_PCs_Jac$ave1-S1_PCs_Jac$axis1)^2+(S1_PCs_Jac$ave2-S1_PCs_Jac$axis2)^2)

outliers<-S1_PCs_Jac[which(S1_PCs_Jac$dist>=0.2382),]
outliers<-S1_PCs_Jac[which(S1_PCs_Jac$dist>=0.21),]
#AlphaDiversity:

AlphaD<-estimate_richness(physeq)


#Generate New Table with outliers removed:
#Use All for traits for Autoencoder, use All1 to generate table for BLUPs and Mapping:
All<-cbind(CMeta[7],CMeta[52],Fam2,gen2,ASV2)
All1<-cbind(CMeta,AlphaD,Fam2,gen2,ASV2)
All2<-All[-which(rownames(All) %in% rownames(outliers)),]
All2<-All1[-which(rownames(All1) %in% rownames(outliers)),]
#write.csv(All2,"SAP_S2_FamGenASVAlpha_RA_filtered2.csv")
#EStrict=x>5,A=.75OutliersRemoved
#Strict=x>5,A=.5
#Less=x>3,A=.2
#Barely=x>1,A=.1

histogram(All2$Faecalibacterium)
histogram(BLUEsdf2$Faecalibacterium)


#Generate a stacked barchart
library(wesanderson)
FBB<-read.csv("FBB0.csv")
ggplot(FBB,aes(x=Subject,y=Value,fill=Genus))+ 
  geom_bar(position="fill", stat="identity")+ theme_classic()+
  theme(text=element_text(size=20,family = "sans"))+ 
  scale_fill_manual(values=wes_palette(n=57, name="IsleofDogs1",type=c("continuous")))+ 
  scale_y_continuous(expand = c(0,0)) +ylab(NULL)







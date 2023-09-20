setwd("~/SAP/GWAS/BigGWASOut//")
library(dplyr)
#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(rMVP)
library(ggplot2)
library(tidyverse)
#library(ggpubr)
library(circlize)
library(ComplexHeatmap)
#Make a list of each trait per/type
S1<-list.files("./S1/",pattern=".csv")
S1<-sub("*_RMIP.csv*","",S1)

S2<-list.files("./S2/",pattern=".csv")
S2<-sub("*_RMIP.csv*","",S2)

S1P<-list.files("./S1_Polymicrobial/",pattern=".csv")
S1P<-sub("*_RMIP.csv*","",S1P)

S2P<-list.files("./S2_Polymicrobial/",pattern=".csv")
S2P<-sub("*_RMIP.csv*","",S2P)

P<-list.files("./PolyMicrobial/",pattern=".csv")
P<-sub("*_RMIP.csv*","",P)

Biochemical<-list.files("./Biochemical/",pattern=".csv")
Biochemical<-sub("*_RMIP.csv*","",Biochemical)

Seed<-list.files("./Seed/",pattern=".csv")
Seed<-sub("*_RMIP.csv*","",Seed)

Agronomic<-list.files("./Agronomic/",pattern=".csv")
Agronomic<-sub("*_RMIP.csv*","",Agronomic)

#Pull just color data from biochemical:
#Color<-Biochemical[grep("Color",Biochemical)]

#Combine all results:
GWAS<-data.frame()

CombineGWAS <- function(trait,type,name){
  Results<-read.csv(paste0("./",type,"/",trait,"_RMIP.csv"))
  Results1<-data.frame(do.call('rbind', strsplit(as.character(Results$Var1),'_',fixed=TRUE)))
  Results2<-data.frame(do.call('rbind', strsplit(as.character(Results1$X1),'S',fixed=TRUE)))
  gwasResults<-as.data.frame(cbind(as.numeric(Results2$X2),as.numeric(Results1$X2),as.numeric(Results$Freq)))
  gwasResultsMeta<-as.data.frame(cbind(as.numeric(Results2$X2),as.numeric(Results1$X2),as.numeric(Results$Freq),paste0(trait),paste0(name)))
  names(gwasResultsMeta)<-c("CHR","BP","FREQ","Trait","Type")
  names(gwasResults)<-c("CHR","BP","FREQ")
  return(gwasResultsMeta)
}

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

for (i in S1){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"S1","Subject1"))
}

for (i in S2){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"S2","Subject2"))
}

for (i in S1P){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"S1_Polymicrobial","Subject1"))
}

for (i in S2P){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"S2_Polymicrobial","Subject2"))
}

#for (i in P){ 
#  GWAS<-rbind(GWAS,CombineGWAS(i,"PolyMicrobial","PrebioticIndex"))
#}


for (i in Agronomic){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"Agronomic","Agronomic"))
}

for (i in Seed){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"Seed","Seed"))
}

for (i in Biochemical){ 
  GWAS<-rbind(GWAS,CombineGWAS(i,"Biochemical","Biochemical"))
}


##Color only:
#for (i in Color){ 
#  GWAS<-rbind(GWAS,CombineGWAS(i,"Biochemical","Color"))
#}


GWAS["FREQ"][GWAS["FREQ"]=="101"]<-"100"

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
sig<-10

#Bin Chromosomes: 
binnum<-75

data_bins <- GWAS2 %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(chrsize = max_bp) %>%
  mutate(binsize=chrsize/binnum) %>%
  select(CHR, chrsize,binsize)  

GWAS3 <- GWAS2 %>%
  inner_join(data_bins, by = "CHR")


GWAS3$Bin<- ceiling(GWAS3$BP/GWAS3$binsize)  

#Remove snps with freq less than some value
GWAS4<-subset(GWAS3,FREQ>=10)
#write.csv(GWAS4,"GWASresults.csv")
A<-nrow(subset(GWAS4,Type=="Agronomic"))
B<-nrow(subset(GWAS4,Type=="Biochemical"))
S<-nrow(subset(GWAS4,Type=="Seed"))
S1<-subset(GWAS4,Type=="Subject1")
S2<-subset(GWAS4,Type=="Subject2")

GWAS4.1<-subset(GWAS4,Type!="Agronomic")
GWAS4.1<-subset(GWAS4.1,Type!="Biochemical")
GWAS4.1<-subset(GWAS4.1,Type!="Seed")
#Most common significant snp?
#SNP<-find_mode(GWAS4$BP)
#print(SNP)

#Or don't:
GWAS4a<-GWAS3
GWAS4.2<-subset(GWAS4a,Type!="Agronomic")
GWAS4.2<-subset(GWAS4.2,Type!="Biochemical")
GWAS4.2<-subset(GWAS4.2,Type!="Seed")



#Calculate number or significant SNPs in a givin bin:
bins2 <- GWAS4 %>%
  group_by(CHR) %>%
  count(Bin)

bins2$Bin2<-paste(bins2$CHR,bins2$Bin,sep="_")

data_cum2 <- bins2 %>%
  group_by(CHR) %>%
  summarise(max_bin = max(Bin)) %>%
  mutate(bin_add = lag(cumsum(max_bin), default = 0)) %>%
  select(CHR, bin_add)

bins3 <- bins2 %>%
  inner_join(data_cum2, by = "CHR") %>%
  mutate(bin_cum = Bin + bin_add)

axis_set2 <- bins3 %>%
  group_by(CHR) %>%
  summarize(center = mean(bin_cum))

GWAS4$Bin2<-paste(GWAS4$CHR,GWAS4$Bin,sep="_")

bins4<-bins3
bins4$CHR<-NULL
bins4$Bin<-NULL

GWAS5 <- GWAS4 %>%
  inner_join(bins4, by = "Bin2")

axis_set3 <- GWAS5 %>%
  group_by(CHR) %>%
  summarize(center = mean(bin_cum))

axis_set4 <- GWAS5 %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

GWAS5$n2<-1

#Count type and return as new column
binF<-read.csv("BinsF.csv")

GWAS4.1<-subset(GWAS4,Type=="Subject1")
binsX <- GWAS4.1 %>%
  group_by(CHR) %>%
  count(Bin2)
names(binsX)<-c("CHR","Bin2","S1")
binF1<-merge(binF,binsX[2:3],all=TRUE)

GWAS4.1<-subset(GWAS4,Type=="Subject2")
binsX <- GWAS4.1 %>%
  group_by(CHR) %>%
  count(Bin2)
names(binsX)<-c("CHR","Bin2","S2")
binF1<-merge(binF1,binsX[2:3],all=TRUE)

GWAS4.1<-subset(GWAS4,Type=="Agronomic")
binsX <- GWAS4.1 %>%
  group_by(CHR) %>%
  count(Bin2)
names(binsX)<-c("CHR","Bin2","Agrinomic")
binF1<-merge(binF1,binsX[2:3],all=TRUE)

GWAS4.1<-subset(GWAS4,Type=="Biochemical")
binsX <- GWAS4.1 %>%
  group_by(CHR) %>%
  count(Bin2)
names(binsX)<-c("CHR","Bin2","Biochemical")
binF1<-merge(binF1,binsX[2:3],all=TRUE)

GWAS4.1<-subset(GWAS4,Type=="Seed")
binsX <- GWAS4.1 %>%
  group_by(CHR) %>%
  count(Bin2)
names(binsX)<-c("CHR","Bin2","Seed")
binF1<-merge(binF1,binsX[2:3],all=TRUE)
binF1[is.na(binF1)]<-0
binF1<-binF1[order(binF1$bin_cum),]

binS<-as.matrix(binF1[6:7])
binA<-as.matrix(binF1[8:10])
binM<-as.matrix(binF1[5])

Genes<-read.csv("../Genes/SorghumGenes.csv")
Genes<-subset(Genes,Type=="gene")
Genes<-Genes[!grepl("super", Genes$Chr),]
data_bins <- Genes %>%
  group_by(Chr) %>%
  summarise(max_bp = max(Stop)) %>%
  mutate(chrsize = max_bp) %>%
  mutate(binsize=chrsize/binnum) %>%
  select(Chr, chrsize,binsize)  
Genes1 <- Genes[1:6] %>%
  inner_join(data_bins, by = "Chr")
Genes1$Bin<- ceiling(Genes1$Stop/Genes1$binsize) 
binsG <- Genes1 %>%
  group_by(Chr) %>%
  count(Bin)
binsG$Bin2<-paste(binsG$Chr,binsG$Bin,sep="_")

density<-read.csv("SNPsInBins.csv",header = TRUE)

binF0<-read.csv("BinsF.csv")
binD<-merge(binF0,density,all=TRUE)
binD<-as.matrix(binD[7])
binD[is.na(binD)]<-0

binF1<-read.csv("BinsF.csv")
binsG<-merge(binF1,binsG,all=TRUE)
binsG<-binsG[order(binsG$bin_cum),]
binG<-as.matrix(binsG[6])
binG[is.na(binG)]<-0

split2 = factor(binF1$CHR)
col_fun0 = colorRamp2(c(0, 1500,9500, 14687), c("white","yellow2", "red", "firebrick4"))
col_fun1 = colorRamp2(c(0,5,15), c("white", "cyan2","blue4"))
col_fun2 = colorRamp2(c(0,6,22), c("white", "green2","darkgreen"))
col_fun3 = colorRamp2(c(0,50,100,164), c("white", "orange","turquoise","forestgreen"))



circos.par(gap.after = c(2, 2, 2, 2,2,2,2,2,2, 8))

circos.heatmap(binM,col=col_fun1,split=split2,track.height=.03,show.sector.labels = FALSE,cluster=FALSE,bg.border="black",bg.lwd = 2)
circos.heatmap(binS,col=col_fun1,split=split2,track.height=.14,show.sector.labels = FALSE,cluster=FALSE)
circos.heatmap(binA,col=col_fun2,track.height=.2,cluster=FALSE)
circos.heatmap(binG,col=col_fun3,track.height=.08,cluster=FALSE)
circos.heatmap(binD,col=col_fun0,track.height=.08,cluster=FALSE)

circos.clear()

circos.par(gap.after = c(2, 2, 2, 2,2,2,2,2,2, 8))
circos.heatmap(dM,col=col_fun0,split=split1,track.height=.1,cluster=FALSE)
circos.clear()

#draw legends
lgn_S<-Legend(title="Microbiome Metrics",col_fun =col_fun1 )
lgn_A<-Legend(title="Agrinomic Metrics",col_fun =col_fun2 )
lgn_G<-Legend(title="Gene Density",col_fun =col_fun3)
lgn_M<-Legend(title="SNP Density",col_fun =col_fun0)
lgn_list<-packLegend(lgn_S,lgn_A,lgn_G,lgn_M)
draw(lgn_list)


#bin5 <- GWAS5 %>%
#  group_by(bin_cum) %>%
#  summarize(occurrence = n())

mat1 = rbind(cbind(matrix(rnorm(50*5, mean = 1), nr = 50), 
                   matrix(rnorm(50*5, mean = -1), nr = 50)),
             cbind(matrix(rnorm(50*5, mean = -1), nr = 50), 
                   matrix(rnorm(50*5, mean = 1), nr = 50))
)
rownames(mat1) = paste0("R", 1:100)
colnames(mat1) = paste0("C", 1:10)
mat1 = mat1[sample(100, 100), ] # randomly permute rows
split = sample(letters[1:5], 100, replace = TRUE)
split = factor(split, levels = letters[1:5])




#GWAS5$Trait<-factor(GWAS5$Trait,levels=c("Color_R","Color_G", "Color_B","Color_PC1","Color_PC2","Color_PC3","PericarpColor_D"))
GWAS5$CHR2<-as.character(GWAS5$CHR)
#Plot RMIP Color only:
ggplot(GWAS5, aes(x=bp_cum,y=FREQ))+ geom_point(aes(color=CHR2),size=2)+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  scale_color_manual(values = rep(c("grey","black","black","grey","black","grey","black","grey","black"),unique(length(axis_set3$CHR))))+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+ labs(fill='Trait Type')+geom_hline(yintercept=10, linetype="dashed",color="black")



ggplot(GWAS5, aes(x=bp_cum,y=FREQ,shape=CHR2))+ geom_point(aes(color=Type),size=2)+scale_shape_manual(values=c(16,17,17,16,17,16,17,16,17,16,17),guide="none")+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  scale_color_manual(values = rep(c("red1","ForestGreen","blue3","aquamarine2","brown4","orange","black"),unique(length(axis_set3$CHR))))+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+ labs(fill='Trait Type')+geom_hline(yintercept=10, linetype="dashed",color="black")

#Make a supplemental Plot:
GWASS1<-subset(GWAS5,Trait=="Fiber")

GWASS1<-subset(GWAS4a,Trait=="S1_Faecalibacterium")
GWASS1<-rbind(GWASS1,subset(GWAS4a,Trait=="S1_ASV25_Faecalibacterium"))
GWASS1<-rbind(GWASS1,subset(GWAS4a,Trait=="S1_PreIndBT"))
GWASS1<-rbind(GWASS1,subset(GWAS4a,Trait=="S1_PC2"))
GWASS1<-rbind(GWASS1,subset(GWAS4a,Trait=="S1_LV2"))

GWASS2<-subset(GWAS4a,Trait=="S2_Faecalibacterium")
GWASS2<-rbind(GWASS2,subset(GWAS4a,Trait=="S2_ASV2_Faecalibacterium"))
GWASS2<-rbind(GWASS2,subset(GWAS4a,Trait=="S2_PreIndBT"))
GWASS2<-rbind(GWASS2,subset(GWAS4a,Trait=="S2_LV10"))
GWASS2<-rbind(GWASS2,subset(GWAS4a,Trait=="S2_PC1"))

GWASBio<-subset(GWAS4a,Trait=="Tannins_NK")
GWASBio<-rbind(GWASBio,subset(GWAS4a,Trait=="StarchPCT"))
GWASBio<-rbind(GWASBio,subset(GWAS4a,Trait=="FiberPCT"))
GWASBio<-rbind(GWASBio,subset(GWAS4a,Trait=="Phenols"))
GWASBio<-rbind(GWASBio,subset(GWAS4a,Trait=="SeedWeight_J"))

a<-ggplot(GWASS1, aes(x=bp_cum,y=FREQ))+ geom_point(aes(color=Trait),size=2,alpha=.8)+scale_shape_manual(values=c(16,17,17,16,17,16,17,16,17,16,17),guide="none")+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 78)) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  scale_color_manual(values = rep(c("ForestGreen","turquoise","royalblue","purple","brown4","red1"),unique(length(axis_set3$CHR))))+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+ labs(fill='Trait Type')+geom_hline(yintercept=10, linetype="dashed",color="black")
  #geom_vline(xintercept=57743537,linetype="dashed",color="black")+geom_vline(xintercept=313785640,linetype="dashed",color="black")
b<-ggplot(GWASS2, aes(x=bp_cum,y=FREQ))+ geom_point(aes(color=Trait),size=2,alpha=.8)+scale_shape_manual(values=c(16,17,17,16,17,16,17,16,17,16,17),guide="none")+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 78)) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  scale_color_manual(values = rep(c("ForestGreen","turquoise","royalblue","purple","brown4","red1"),unique(length(axis_set3$CHR))))+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+ labs(fill='Trait Type')+geom_hline(yintercept=10, linetype="dashed",color="black")
c<-ggplot(GWASBio, aes(x=bp_cum,y=FREQ))+ geom_point(aes(color=Trait),size=2,alpha=.8)+scale_shape_manual(values=c(16,17,17,16,17,16,17,16,17,16,17),guide="none")+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 78)) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+
  scale_color_manual(values = rep(c("brown4","red1","chocolate","coral3","burlywood2"),unique(length(axis_set3$CHR))))+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5)
  )+theme(text=element_text(size=24,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+ labs(fill='Trait Type')+geom_hline(yintercept=10, linetype="dashed",color="black")

N<-ggarrange(a,b,c,ncol = 1,nrow = 3)
ggsave("SupManWleg.svg",device = "svg",path="../../Figures/",plot=N,dpi=300,width=480,height = 255,units="mm")

a<-a+theme(legend.position = "none")
b<-b+theme(legend.position = "none")
c<-c+theme(legend.position = "none")

N<-ggarrange(a,b,c,ncol = 1,nrow = 3)
ggsave("SupMan.svg",device = "svg",path="../../Figures/",plot=N,dpi=300,width=480,height = 255,units="mm")

GWAS9<-GWAS5

table(GWAS9$Type)

GWAS10<-rbind(GWAS9,subset(GWAS9,Type=='Subject1'))
GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Subject1'))
GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Subject1'))
GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Subject2'))
GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Subject2'))
#GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Subject1'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S1_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S1_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S1_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S1_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='S1_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S1_PCs'))

#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S2_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S2_PCs'))
#GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='S2_PCs'))

GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Biochemical')) 
GWAS10<-rbind(GWAS10,subset(GWAS9,Type=='Biochemical')) 
GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='Seed')) 
GWAS10<-rbind(GWAS10,subset(GWAS10,Type=='Agronomic')) 

table(GWAS10$Type)


GWAS10$Subject<-GWAS10$bp_cum
GWAS10$Subject[GWAS10$Type=="Agronomic"]<-NA
GWAS10$Subject[GWAS10$Type=="Biochemical"]<-NA
GWAS10$Subject[GWAS10$Type=="Seed"]<-NA
#GWAS10$Subject[GWAS10$Type=="PM_Subject1"]<-NA
#GWAS10$Subject[GWAS10$Type=="PM_Subject2"]<-NA

GWAS10$Agri<-GWAS10$bp_cum
GWAS10$Agri[GWAS10$Type=="Subject1"]<-NA
GWAS10$Agri[GWAS10$Type=="Subject2"]<-NA
#GWAS10$Agri[GWAS10$Type=="S1_PCs"]<-NA
#GWAS10$Agri[GWAS10$Type=="S2_PCs"]<-NA

#GWAS10$PM<-GWAS10$bp_cum
#GWAS10$PM[GWAS10$Type=="Subject1"]<-NA
#GWAS10$PM[GWAS10$Type=="Subject2"]<-NA
#GWAS10$PM[GWAS10$Type=="Agronomic"]<-NA
#GWAS10$PM[GWAS10$Type=="Biochemical"]<-NA
#GWAS10$PM[GWAS10$Type=="Seed"]<-NA

#Plot mirrored density plot
a<-ggplot(GWAS10, aes(x=x)) + geom_density(aes(x=Subject, y=..density..,fill=Type),color=NA,adjust=.013,position="stack")+ 
  geom_density(aes(x=Agri, y=-..density..,fill=Type),color=NA,adjust=.02,position="stack")+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + 
  scale_fill_manual(values = rep(c("chartreuse","green3","darkgreen","cyan2","blue3","violet","pink","darkgreen","red1","aquamarine2"), unique(length(axis_set3$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) + ylab("density") + xlab(NULL)+
  theme(panel.border = element_blank(),panel.grid = element_blank(),panel.background = element_blank(),
                          axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5),axis.text.y=element_blank()
  )+theme(text=element_text(size=16,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+ylab(NULL)


axis_setmin <- GWAS5 %>%
  group_by(CHR) %>%
  summarize(val = min(bp_cum))
axis_setmax <- GWAS5 %>%
  group_by(CHR) %>%
  summarize(val = max(bp_cum))

axis_set5<-rbind(axis_setmin,axis_setmax)
library(svglite)
b<-a+theme(axis.ticks = element_blank())
b
ggsave("Figure3.svg",device = "svg",path="../../Figures/",plot=b,dpi=300,width=180,height = 115,units="mm")


c<-a+scale_x_continuous(breaks=axis_setmax$val)
ggsave("Figure3.1.svg",device = "svg",path="../../Figures/",plot=c,dpi=300,width=180,height = 1115,units="mm")







##########################################################################################################3
#Insert LD info
thresh<-10
OutF<-GWAS4.2
#OutF<-GWAS5

LD<-read.table("../BigGWASOut/LD/MEL4A.1_LD.txt",header = TRUE)
SNP<-as.numeric(LD[1,2])
Site<-as.numeric(LD[1,3])

LD2<-LD[14:15]
names(LD2)<-c("R2","DPrime")


LDpos1<-as.data.frame(LD$Position1)
LDpos1<-as.data.frame(LDpos1[-(1:Site),])
names(LDpos1)<-"pos"
LDpos2<-as.data.frame(LD$Position2)
LDpos2<-as.data.frame(LDpos2[1:Site,])
names(LDpos2)<-"pos"
LDpos<-as.data.frame(rbind(LDpos2,LDpos1))
LD2$pos<-LDpos$pos

#insert missing row:
Missing1<-c(1,1,SNP)
LD2<-rbind(LD2,Missing1)

#Determine Region of interest:
LD3<-subset(LD2,R2>0.5)
LD3$DPrime<-NULL

LD4<-LD3[-nrow(LD3),]
Range<-paste0(LD4[1,2],"..",LD4[nrow(LD4),2])
print(Range)



OutF2<-merge(OutF,LD2,by.x='BP',by.y="pos")

#Remove any problem traits
#OutF2<-subset(OutF2,Trait!="Phenols")

OutF2<-OutF2[order(OutF2$R2),]
OutF2$BP<-OutF2$BP/1000000

OutF2<-subset(OutF2,CHR==6)

OutF3<-subset(OutF2,BP<=43.35)

OutF3<-subset(OutF3,BP>=43.2772)

ggplot(data=OutF3,aes(x = BP, y = FREQ, color=R2))+geom_point(size=4)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))+
  geom_hline(yintercept=thresh, linetype="dashed",color="black")+ylab("RMIP Frequency")+xlab("Position (MB)")+
  scale_color_gradient2(low = "black",mid="blue", high = "orange",midpoint=0.1) + scale_y_continuous(expand = c(0, 0),limits = c(0,35))


#write.csv(OutF3,"MEL6A_gwasresult_Zoomed.csv")

OutF3<-read.csv("MEL6A_gwasresult_Zoomed.csv")
OutF3<-subset(OutF3,Family!="Rikenellaceae")
OutF3<-subset(OutF3,Family!="Sutterellaceae")
OutF3<-subset(OutF3,Family!="Tannerellaceae")
OutF3<-subset(OutF3,Family!="Burkholderiaceae")
OutF3<-subset(OutF3,Family!="Christensenellaceae")
OutF3<-subset(OutF3,Family!="Odoribacteraceae")

#OutF3<-subset(OutF3,BP<=43.331000)

ggplot(data=OutF3,aes(x = BP, y = FREQ, color=Family))+geom_jitter(size=3,alpha=0.8,height=0,width=0)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))+
  geom_hline(yintercept=thresh, linetype="dashed",color="black")+ylab("RMIP Frequency")+xlab("Position (MB)")+
  scale_y_continuous(expand = c(0, 0),limits = c(0,35))+geom_vline(xintercept=43.277256)+#+scale_x_continuous(expand = c(0, 0),limits = c(43.32945,43.3305))
  geom_vline(xintercept=43.279871)+geom_vline(xintercept=43.304436)+geom_vline(xintercept=43.308204)+
  geom_vline(xintercept=43.320306)+geom_vline(xintercept=43.324194)+geom_vline(xintercept=43.344453)+
  geom_vline(xintercept=43.349422)


svg("MEL6ARMIP.svg",height = 8,width = 14)
ggplot(data=OutF3,aes(x = BP, y = FREQ, color=Family))+geom_jitter(size=3,alpha=0.75,height=.4,width=.0002)+ theme_classic()+theme(text=element_text(size=20,family = "sans"))+
    geom_hline(yintercept=thresh, linetype="dashed",color="black")+ylab("RMIP Frequency")+xlab("Position (MB)")+
    scale_y_continuous(expand = c(0, 0),limits = c(0,35))+scale_color_manual(values = c("#332288","#44AA99","#661100","#D55E00"))#+scale_x_continuous(expand = c(0, 0),limits = c(43.32945,43.3305))
dev.off()

PlotPrepELD<-function(Trait){
  Out1<-as.data.frame(OutF2$pos/1000000)
  names(Out1)<-"position"
  Out1$effect<-OutF2[[paste0(Trait,"_Effect")]]
  #Out1$effect<-abs(Out1$effect)
  Out1$p<--log10(OutF2[[paste0(Trait,"_V2.FarmCPU")]])
  Out1$R2<-OutF2$R2
  ggplot(data=Out1,aes(x=position,y=p,color=R2))+geom_point()+ theme_classic()+theme(text=element_text(size=20,family = "sans"))+
    geom_hline(yintercept=thresh, linetype="dashed",color="black")+ylab(bquote(-log[10](P)))+xlab("Position (MB)")+
    ggtitle(paste0(Trait))+ scale_y_continuous(expand = c(0, 0.15))+scale_color_gradient2(low = "grey",mid="blue", high = "orange",midpoint=0.2) 
}

#plots LD SiteByAll:
#for (i in list2){
#  print(PlotPrepELD(i))
#}

#write plot?
#for (i in list2){
#  ggsave(filename=paste0(i,".png"),plot=print(PlotPrepEM(i),scale=1))
#}


#Determine Region of interest:
LD3<-subset(LD2,R2>0.5)
LD3$DPrime<-NULL

LD4<-LD3[-nrow(LD3),]
Range<-paste0(LD4[1,2],"..",LD4[nrow(LD4),2])
print(Range)


#Plot LD matrix

LD<-read.table("LD_SNPsInMEL2/MEL10A_LD.txt",head=T)
#LD<-read.table("LD/MEL6A_e1_LD.txt",head=T)

#LD$Position1<-as.character(LD$Position1)
#LD$Position2<-as.character(LD$Position2)


ggplot(LD,aes(x=reorder(Position1,sort(as.numeric(Position1))),y=reorder(Position2,sort(as.numeric(Position2))),fill=R.2))+geom_tile(colour="black")+ theme_classic()+theme(text=element_text(size=12,family = "sans"))+
  scale_fill_gradientn(limits=c(0,1),colors=c("white","red","purple3","blue3"))+guides(fill=guide_colorbar(ticks = FALSE,title=expression(paste(R^2))))

ggplot(LD,aes(x=reorder(Position1,sort(as.numeric(Position1))),y=reorder(Position2,sort(as.numeric(Position2))),fill=R.2))+geom_tile(colour="black")+ theme_classic()+theme(text=element_text(size=12,family = "sans"))+
  scale_fill_gradientn(limits=c(0,1),colors=c("white","red","purple3","blue3"))+guides(fill=guide_colorbar(ticks = FALSE,title=expression(paste(R^2))))+
  theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.line = element_blank(),axis.text=element_blank(),legend.position = "none")

#Remove some snps:
#MEL1
LD2<-subset(LD,Position1>=54000000)
LD2<-subset(LD2,Position2>=54000000)
#MEL2
LD2<-subset(LD,Position1<=8709020)
LD2<-subset(LD2,Position2<=8709020)
LD2<-subset(LD2,Position1>=4678862)
LD2<-subset(LD2,Position2>=4678862)
#MEL3B
LD2<-subset(LD,Position1<=54709000)
LD2<-subset(LD2,Position2<=54709000)
#MEL4A
LD2<-subset(LD,Position1<=8600000)
LD2<-subset(LD2,Position2<=8600000)
#MEL4B
LD2<-subset(LD,Position1<=62890190)
LD2<-subset(LD2,Position2<=62890190)
LD2<-subset(LD2,Position1>=61942457)
LD2<-subset(LD2,Position2>=61942457)
#MEL5A
LD2<-subset(LD,Position1>=6800000)
LD2<-subset(LD2,Position2>=6800000)
#MEL5B
LD2<-subset(LD,Position1>=51350000)
LD2<-subset(LD2,Position2>=51350000)
#MEL6A
LD2<-subset(LD,Position1<=44073121)
#MEL6C
LD2<-subset(LD,Position1>=58260000)
LD2<-subset(LD2,Position2>=58260000)
#MEL7
LD2<-subset(LD,Position1>=56183450)
LD2<-subset(LD2,Position2>=56183450)
LD2<-subset(LD2,Position1<=57876870)
LD2<-subset(LD2,Position2<=57876870)
#MEL8
LD2<-subset(LD,Position1>=3800000)
LD2<-subset(LD2,Position2>=3800000)
#MEL10B
LD2<-subset(LD,Position1>=49520000)
LD2<-subset(LD2,Position2>=49520000)
LD2<-subset(LD2,Position1<=51960000)
LD2<-subset(LD2,Position2<=51960000)



a<-ggplot(LD2,aes(x=reorder(Position1,sort(as.numeric(Position1))),y=reorder(Position2,sort(as.numeric(Position2))),fill=R.2))+geom_tile(colour="black")+ theme_classic()+theme(text=element_text(size=12,family = "sans"))+
  scale_fill_gradientn(limits=c(0,1),colors=c("white","red","purple3","blue3"))+guides(fill=guide_colorbar(ticks = FALSE,title=expression(paste(R^2))))+
  theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.line = element_blank(),axis.text=element_blank(),legend.position = "none")
a

ggsave("MEL4ALD.png",plot=a,width=20,height=20,units="cm")
#Calculate local LD?
LD<-read.table("LD_SNPsInMEL2/MEL10B_LD.txt",head=T)
Max<-max(LD$R.2)
R<-mean(LD$R.2)
SNPs<-as.numeric(nrow(as.data.frame(unique(LD$Position1))))
Dist<-as.numeric((max(LD$Position1)-min(LD$Position1)))



#Allele Dist:
#Read in BLUE values:
BLUEs<-read.csv("../SAP_CompiledBLUEs3.csv")
Alleles<-read.csv("../Genes/MEL6A_alleles.csv")

Alleles<-Alleles[-c(1),]

df1<-merge(BLUEs,Alleles)
#df2<-subset(df1,MEL6A!="Het")
df2<-df1
df2.1<-subset(df2,PreIndBT!="#VALUE!")
df2.1<-subset(df2.1,PreIndBT!="NA")
df2.1$PreIndBT<-as.numeric(df2.1$PreIndBT)

a<-ggplot(df2,aes(x=MEL6A,y=S2_Coprococcus.3))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))#+geom_jitter(aes(color=Origin2),size=3)
a+stat_compare_means(label.x=1.4)

#Figure2 Pop structure of sorghum influences microbiome
ggplot(df2.1)+geom_density(aes(x=S2_PreIndBT),fill="blue2", alpha = 0.9,color=NA)+
  geom_density(aes(x=S1_PreIndBT),fill="cyan2", alpha = 0.8,color=NA)+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+xlab("Butyrate(mMol)")+
  coord_cartesian(ylim = c(0.476, 10))

ggplot(df2)+geom_histogram(aes(x=S2_Faecalibacterium),fill="blue2", alpha = 0.9,color=NA)+
  theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+xlab("S2_Faecalibacterium")

df2.1<-subset(df2,S1_PreIndBT!="#VALUE!")
df2.1<-subset(df2.1,S1_PreIndBT!="NA")
df2.1$S1_PreIndBT<-as.numeric(df2.1$S1_PreIndBT)

PRhist<-ggplot(df2.1)+geom_density(aes(x=S2_PreIndBT),fill="blue2", alpha = 0.9,color=NA)+
  geom_density(aes(x=S1_PreIndBT),fill="cyan2", alpha = 0.8,color=NA)+theme_classic()+
  theme(text=element_text(size=30,family = "sans"))+xlab("Prebiotic Index")+
  coord_cartesian(ylim = c(0.476, 10))

Buthist<-ggplot(df2.1)+geom_density(aes(x=S2_Butyrate),fill="blue2", alpha = 0.9,color=NA)+
  geom_density(aes(x=S1_Butyrate),fill="cyan2", alpha = 0.8,color=NA)+theme_classic()+
  theme(text=element_text(size=30,family = "sans"))+xlab("Butyrate(mMol)")+
  coord_cartesian(ylim = c(0.019, .4))

ggsave("PRhist.svg",plot=PRhist,width=19,height=15,units="cm")
ggsave("Buthist.svg",plot=Buthist,width=19,height=15,units="cm")



ggplot(df2.2,aes(x=ProteinPCT,y=S2_Bacteroides))+theme_classic()+geom_point()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="spearman",cor.coef.name="rho")+xlab("Protein Content")+ylab("Subject2 Bacteroides")

ggplot(df2.2,aes(x=SeedOil_T,y=S1_Akkermansia))+theme_classic()+geom_point()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+xlab("Oil Content")+ylab("Subject2 Akkermansia")

ggplot(df2.2,aes(x=SeedOil_T,y=S1_Akkermansia))+theme_classic()+geom_point()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="spearman",cor.coef.name="rho")+xlab("Oil Content")+ylab("Subject2 Akkermansia")

ggplot(df2.2,aes(x=SeedAcidFiber_T1,y=S2_PreIndBT))+theme_classic()+geom_point()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="spearman",cor.coef.name="rho")+xlab("Oil Content")+ylab("Subject2 Akkermansia")



df2$SeedAcidFiber_T1<-df2$SeedAcidFiber_T+3
df2$SeedOil_T<-df2$SeedOil_T+3
df2.2<-subset(df2,Tannins_NK<=20)
df2.2<-subset(df2.2,Tannins_NK>=0)

TanXS1F<-ggplot(df2.2,aes(x=Tannins_NK,y=S1_Faecalibacterium))+theme_classic()+geom_point()+
  theme(text=element_text(size=30,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  xlab("Tannin Content")+ylab("Subject1 Faeclaibacterium")


FibXS2PI<-ggplot(df2.2,aes(x=SeedAcidFiber_T1,y=S2_PreIndBT))+theme_classic()+geom_point()+
  theme(text=element_text(size=30,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  xlab("Acid Soluble Fiber")+ylab("Subject2 Prebiotic Index")

df2.3<-subset(df2,ProteinPCT>=8.50)

OilXS1Akk<-ggplot(df2.2,aes(x=SeedOil_T,y=S1_Akkermansia))+theme_classic()+geom_point()+
  theme(text=element_text(size=30,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  xlab("Oil Content")+ylab("Subject1 Akkermansia")


ggsave("TanXS1F.svg",plot=TanXS1F,width=18,height=16,units="cm")
ggsave("FibXS2PI.svg",plot=FibXS2PI,width=18,height=16,units="cm")
ggsave("OilXS1Akk.svg",plot=OilXS1Akk,width=18,height=16,units="cm")


ggplot(df2)+geom_violin(aes(x=K.Pop,y=S1_Butyrate),fill="cyan",color=NA,alpha=.7)+
  geom_violin(aes(x=K.Pop,y=S2_Butyrate),fill="blue2",color=NA,alpha=.6)+theme_classic()+
  theme(text=element_text(size=24,family = "sans"),axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  ylab("Prebiotic Index")+xlab("Population")

ggplot(df2.1,aes(x=K.Pop,y=PreIndBT))+geom_boxplot(fill="purple",color="black",alpha=1)+
  theme_classic()+
  theme(text=element_text(size=24,family = "sans"),axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  ylab("Prebiotic Index")+xlab("Population")+stat_compare_means(label.x=1.4)


#Look at only non-tannin lines 
#AND GENERATE LIST OF LINES IN POOLS
df4<-subset(df2.1,Tannins_NK<=1)
df4<-subset(df4,Tannins_NK>=-.4)
a<-ggplot(df4,aes(x=MEL6A,y=Tannins_NK))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))

a<-ggplot(df2,aes(x=MEL6A,y=S2_Faecalibacterium))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))


a+stat_compare_means(label.x=1.4)+geom_point()

set.seed(52)
Pool_NT_Minor<-subset(df4,MEL6A=="Alt")
Pool_NT_Minor<-subset(Pool_NT_Minor,PreIndBT<=mean(Pool_NT_Minor$PreIndBT))
Pool_NT_Minor1<-Pool_NT_Minor[sample(nrow(Pool_NT_Minor),10),]
Pool_NT_Minor1$Pool<-"NT_Minor_1"
Pool_NT_MinorA<-Pool_NT_Minor[-which(Pool_NT_Minor$Line %in% Pool_NT_Minor1$Line),]
Pool_NT_Minor2<-Pool_NT_MinorA[sample(nrow(Pool_NT_MinorA),10),]
Pool_NT_Minor2$Pool<-"NT_Minor_2"

Pool_NT_Major<-subset(df4,MEL6A=="Ref")
Pool_NT_Major<-subset(Pool_NT_Major,PreIndBT>=mean(Pool_NT_Major$PreIndBT))
Pool_NT_Major1<-Pool_NT_Major[sample(nrow(Pool_NT_Major),10),]
Pool_NT_Major1$Pool<-"NT_Major_1"
Pool_NT_MajorA<-Pool_NT_Major[-which(Pool_NT_Major$Line %in% Pool_NT_Major1$Line),]
Pool_NT_Major2<-Pool_NT_MajorA[sample(nrow(Pool_NT_MajorA),10),]
Pool_NT_Major2$Pool<-"NT_Major_2"
Pool_NT_MajorB<-Pool_NT_MajorA[-which(Pool_NT_MajorA$Line %in% Pool_NT_Major2$Line),]

Pool_NT<-rbind(Pool_NT_Minor1,Pool_NT_Minor2,Pool_NT_Major1,Pool_NT_Major2)

a<-ggplot(Pool_NT,aes(x=Pool,y=Tannins_NK))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)+geom_point()

#########################################
df4<-subset(df2.1,Tannins_NK<=30)
df4<-subset(df4,Tannins_NK>=1)
a<-ggplot(df4,aes(x=MEL6A,y=PreIndBT))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)+geom_point()

set.seed(50)
Pool_T_Minor<-subset(df4,MEL6A=="Alt")
Pool_T_Minor<-subset(Pool_T_Minor,PreIndBT<=mean(Pool_T_Minor$PreIndBT))
Pool_T_Minor1<-Pool_T_Minor[sample(nrow(Pool_T_Minor),6),]
Pool_T_Minor1$Pool<-"T_Minor_1"
Pool_T_Minor2<-Pool_T_Minor[-which(Pool_T_Minor$Line %in% Pool_T_Minor1$Line),]
Pool_T_Minor2$Pool<-"T_Minor_2"

Pool_T_Major<-subset(df4,MEL6A=="Ref")
Pool_T_Major<-subset(Pool_T_Major,PreIndBT>=mean(Pool_T_Major$PreIndBT))
Pool_T_Major1<-Pool_T_Major[sample(nrow(Pool_T_Major),10),]
Pool_T_Major1$Pool<-"T_Major_1"
Pool_T_MajorA<-Pool_T_Major[-which(Pool_T_Major$Line %in% Pool_T_Major1$Line),]
Pool_T_Major2<-Pool_T_MajorA[sample(nrow(Pool_T_MajorA),10),]
Pool_T_Major2$Pool<-"T_Major_2"

Pool_T<-rbind(Pool_T_Minor1,Pool_T_Minor2,Pool_T_Major1,Pool_T_Major2)

a<-ggplot(Pool_T,aes(x=Pool,y=Tannins_NK))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)+geom_point()


Lines_NT<-as.data.frame(cbind(Pool_NT$Line,Pool_NT$Pool))
Lines_T<-as.data.frame(Pool_T$Line,Pool_T$Pool)

write.csv(Lines_NT,"Non-tanninPool.csv")
write.csv(Lines_T,"TanninPool.csv")

Pool_A<-rbind(Pool_NT,Pool_T)
a<-ggplot(Pool_A,aes(x=Pool,y=Tannins_NK))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)+geom_point()

#Lines for Fructan Measurements:
a<-ggplot(df2.1,aes(x=MEL6A,y=PreIndBT))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)+geom_point()

set.seed(49)
Lines_Minor<-subset(df4,MEL6A=="Alt")
Lines_Minor<-subset(Lines_Minor,PreIndBT<=mean(Lines_Minor$PreIndBT))
Lines_Minor1<-Lines_Minor[sample(nrow(Lines_Minor),6),]

Lines_Major<-subset(df4,MEL6A=="Ref")
Lines_Major<-subset(Lines_Major,PreIndBT>=mean(Lines_Major$PreIndBT))
Lines_Major1<-Lines_Major[sample(nrow(Lines_Major),6),]

Lines<-rbind(Lines_Minor1,Lines_Major1)
a<-ggplot(Lines,aes(x=MEL6A,y=Tannins_NK))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)+geom_point()

Lines1<-as.data.frame(Lines$Line)

write.csv(Lines1,"Pool.csv")

a<-ggplot(df4,aes(x=MEL6A,y=S2_Escherichia.Shigella))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)

mean(df2$Tannins_NK,na.rm = TRUE)

#Plot Haplotypes:
Alleles<-read.csv("../Genes/MEL6A_alleles.csv")
Alleles<-Alleles[-c(1),]
df1<-merge(BLUEs,Alleles)

df1.1<-subset(df1,PreIndBT!="#VALUE!")
df1.1<-subset(df1.1,PreIndBT!="NA")
df1.1$PreIndBT<-as.numeric(df1.1$PreIndBT)

df1.2<-subset(df1.1,Haplotype!="NA")
df1.2<-subset(df1.2,Haplotype!="A7")
df1.2<-subset(df1.2,Haplotype!="A1")
df1.2<-subset(df1.2,Haplotype!="A13")
df1.2<-subset(df1.2,Haplotype!="A10")

c1<-list(c("A1","R"),c("A","A1"),c("A","A2"),c("A1","A2"))

a<-ggplot(df1.2,aes(x=Haplotype2,y=S1_Escherichia.Shigella))+geom_point(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(comparisons = c1)

#MEL6A_32.splice1
#MEL6A_32.splice2
#MEL6A_98.splice
#MEL6A_64.splice1-3
#MEL6A_30.FS1
df1.3<-subset(df1.1,MEL6A_64.Del!="Het")

a<-ggplot(df1.3,aes(x=MEL6A_64.Del,y=S2_PreIndBT))+geom_boxplot(binaxis='y', stackdir='center',stackratio=.9,dotsize = 0.2, position=position_dodge(0.8))+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))
a+stat_compare_means(label.x=1.4)

#Plot by exon1_haplotype
setwd("~/SAP/GWAS/BigGWASOut//")
BLUEs<-read.csv("../SAP_CompiledBLUEs3.csv")
Alleles<-read.csv("../Genes/MEL6A_alleles_64_exon1.csv")
Alleles<-Alleles[c(2:329),]


df1<-merge(BLUEs,Alleles)
df2<-subset(df1,CW1_Haplotype!="Het")

#Calculate log fold change for genus (for biplot)
df3<-cbind(df2[23:57],df2[125:158],df2[473])
CalcRatio<-function(Trait){
  meanMaj<-mean(df3[[Trait]][df3$CW1_Haplotype == "Major"],na.rm = TRUE)
  meanMin<-mean(df3[[Trait]][df3$CW1_Haplotype == "Minor"],na.rm = TRUE)
  Lratio<-as.data.frame(log2((meanMaj+ runif(1,.00005,.00009))/(meanMin+ runif(1,.00005,.00009))))
  rownames(Lratio)<-paste(Trait)
  return(Lratio)
}
X<-CalcRatio("S1_Phascolarctobacterium")
traits<-names(df3)
for (i in traits){
  X<-rbind(X,CalcRatio(i))
}
write.csv(X,"Out1.csv")


df2.1<-subset(df2,PreIndBT!="#VALUE!")
df2.1<-subset(df2.1,PreIndBT!="NA")
df2.1$PreIndBT<-as.numeric(df2.1$PreIndBT)

S1_F_Maj_mean<-mean(df2$S1_Faecalibacterium[df2$CW1_Haplotype == "Major"],na.rm = TRUE)
S1_F_Min_mean<-mean(df2$S1_Faecalibacterium[df2$CW1_Haplotype == "Minor"],na.rm = TRUE)
S1_F_Lratio<-log(S1_F_Maj_mean/S1_F_Min_mean)
df2$Starch
S1<-ggplot(df2,aes(x=CW1_Haplotype,y=StarchPCT))+geom_boxplot()+
  theme_classic()+theme(text=element_text(size=24,family = "sans"))+xlab(label=NULL)+ylab(label="% Starch")
#ggsave("Starch.svg",plot=S1,width=6,height=10,units="cm")

PlotRatio<-function(df,Trait){
  meanMaj<-mean(df2[[Trait]][df2$CW1_Haplotype == "Major"],na.rm = TRUE)
  meanMin<-mean(df2[[Trait]][df2$CW1_Haplotype == "Minor"],na.rm = TRUE)
  Lratio<-log(meanMaj/meanMin)
  #Effect<-log(meanMaj-meanMin)
  Effect<-log(abs(meanMaj-meanMin))*-1
  p<-ggplot(df2,aes(CW1_Haplotype,.data[[Trait]]))+geom_boxplot()+
    theme_classic()+ xlab(label=NULL)+
    ggtitle(paste0("Effect: ", round(Effect, 2)))+
    theme(text=element_text(size=24,family = "sans"))
  return(p)
}



S1F<-PlotRatio(df2,"S1_Faecalibacterium")
S2F<-PlotRatio(df2,"S2_Faecalibacterium")

S1E<-PlotRatio(df2,"S1_Escherichia.Shigella")
S2E<-PlotRatio(df2,"S2_Escherichia.Shigella")

ggsave("S1F.svg",plot=S1F,width=8,height=15,units="cm")
ggsave("S2F.svg",plot=S2F,width=8,height=15,units="cm")
ggsave("S1E.svg",plot=S1E,width=8,height=15,units="cm")
ggsave("S2E.svg",plot=S2E,width=8,height=15,units="cm")



PlotSig<-function(df,Trait){
  meanMaj<-mean(df2[[Trait]][df2$CW1_Haplotype == "Major"],na.rm = TRUE)
  meanMin<-mean(df2[[Trait]][df2$CW1_Haplotype == "Minor"],na.rm = TRUE)
  #Lratio<-log(abs(meanMaj/meanMin))
  #Effect<-log(meanMaj-meanMin)
  #Effect<-log(abs(meanMaj-meanMin))*-1
  Effect<-meanMin-meanMaj
  p<-ggplot(df2,aes(CW1_Haplotype,.data[[Trait]]))+geom_boxplot()+
    theme_classic()+ xlab(label=NULL)+
    ggtitle(paste0("Effect: ", round(Effect, 2)))+
    theme(text=element_text(size=24,family = "sans"))+stat_compare_means()
  return(p)
}

PlotSig(df2,"Starch")
PlotSig(df2,"StarchPCT")
PlotSig(df2,"Starch_T")

PlotSig(df2,"SeedAcidFiber_T")

PlotSig(df2,"ProteinPCT")
#group by alleles
library(dplyr)
dfgrouped<-df3%>%
  group_by(MEL6A_32.splice1, MEL6A_32.splice2,MEL6A_98.splice,MEL6A_64.splice1,MEL6A_64.splice2,MEL6A_64.splice3,MEL6A_30.FS1) %>%
  summarise(n=n())


#MLM
#Sorghum Population Structure:
PopPCs<-read.table("../Genes/Reseq_PCs/SAP_PCs.txt",header=TRUE)
PopPCs$Line<-sub("PI_","PI",PopPCs$Line)

df3<-merge(df2,PopPCs,by.x="Line",by.y="Line")

#Plot Popstructure by Prebiotic Index
df3.1<-subset(df3,PreIndBT!="#VALUE!")
df3.1$PreIndBT<-as.numeric(df3.1$PreIndBT)
ggplot(df3.1,aes(x=PC1,y=PC3,color=PreIndBT))+geom_point()+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+scale_color_gradient2(na.value="white",low="black",mid="blue4",high = "violet")


ggplot(df3,aes(x=PC1,y=PC2,color=S2_PreIndBT))+geom_point()+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+scale_color_gradient2(na.value="white",low="black",mid="blue4",high = "violet")

df3.1<-subset(df3,S1_PreIndBT!="#VALUE!")
df3.1$S1_PreIndBT<-as.numeric(df3.1$S1_PreIndBT)
ggplot(df3.1,aes(x=PC1,y=PC2,color=S1_PreIndBT))+geom_point()+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+scale_color_gradient2(na.value="white",low="black",mid="blue4",high = "violet")

#Tannins an popstructure?
#df3.1<-subset(df3,Tannins_NK<=35)
#ggplot(df3.1,aes(x=PC1,y=PC2,color=Tannins_NK))+geom_point()+theme_classic()+
#  theme(text=element_text(size=24,family = "sans"))+scale_color_gradient2(na.value="white",low="black",mid="blue4",high = "violet",midpoint = 2)

#Butyrate correlated with preindex
ggplot(df2.1,aes(x=S2_PreIndBT,y=S1_Butyrate))+geom_point()+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
#ggsave("S1ButXS1PI.svg",width = 15,height=15,units = "cm")

ggplot(df3.1,aes(x=S2_Faecalibacterium,y=Tannins_NK))+geom_point()+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))



df3.1<-subset(df3,Line!="PI656029")

ggplot(df3.1,aes(x=S2_PreIndBT,y=S2_Butyrate))+geom_point()+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+geom_smooth(method="lm",color="purple2", se = FALSE)+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
ggsave("S2ButXS2PI.svg",width = 15,height=15,units = "cm")

#fitH <-mmer(S2_Butyrate~MEL2B_SS1a, random=~PC1.y+PC2.y+PC3.y, rcov=~units, data=df1.3a) 
#summary(fitH)

fitG1<-glm(S1_Faecalibacterium~MEL6A+Tannins_NK+MEL6A*Tannins_NK,data=df3)
summary(fitG1)
fitG2<-glm(S1_Faecalibacterium~MEL6A+PC1+PC2+PC3,data=df3)
summary(fitG2)
p<-coef(summary(fitG2))[2,4]
fitG2<-glm(S1_Faecalibacterium~MEL6A+PC1+PC2+PC3+Tannins_NK,data=df3)
summary(fitG2)

fitG2<-glm(S2_Escherichia.Shigella~MEL6A+PC1+PC2+PC3+Tannins_NK+MEL6A*Tannins_NK,data=df3)
summary(fitG2)


#read SNPs:
SNPs<-read.csv("./BigMELsubsets/MEL6A.csv")
SNPs[1:5,1:5]
SNPs$Alt<-apply(SNPs,1,function(x) length(which(x=="1|1")))
SNPs1<-SNPs[SNPs$Alt >= 33,]
SNPs1$Alt<-NULL
SNPs1$Alt<-apply(SNPs1,1,function(x) length(which(x=="0|0")))
SNPs2<-SNPs1[SNPs1$Alt >= 33,]
SNPs2$Alt<-NULL
SNPs1<-SNPs2

SNPst<-as.data.frame(t(SNPs1[-c(1:4)]))
names(SNPst)<-SNPs1$SNP
SNPst[1:5,1:5]

SNPsList<-as.data.frame(names(SNPst))
names(SNPsList)<-"SNPs"


SNPst[SNPst=="0|0"]<-"Ref"
SNPst[SNPst=="0|1"]<-"Het"
SNPst[SNPst=="1|0"]<-"Het"
SNPst[SNPst=="1|1"]<-"Alt"

Metrics<-df1
#Metrics$Line<-gsub("_","",Metrics$Line)
SNPst$Line<-rownames(SNPst)
df<-merge(SNPst,Metrics)
df<-merge(df,PopPCs)

res1<-aggregate(df$S1_Butyrate,df[2],FUN=mean, na.rm = TRUE)
res2<-as.data.frame(t(res1[2]))
names(res2)<-res1$S06_41309945
ave<-mean(df$S1_Butyrate,na.rm=TRUE)
val<-(res2$Ref)/(res2$Alt)
maf<-((2*sum(df[2]=="Alt"))+sum(df[2]=="Het"))/656
Alt1<-subset(df,df[2]=="Alt")
Ref1<-subset(df,df[2]=="Ref")
#wil<-wilcox.test(Alt1$S2_Butyrate,Ref1$S2_Butyrate)
#p<-wil$p.value
fitG2<-glm(S1_Butyrate~S06_41309945+PC1+PC2+PC3,data=df)
p<-coef(summary(fitG2))[2,4]

outF<-as.data.frame(cbind(val,p,maf))

CalcEfx<-function(SNP,Trait){
  res<-aggregate(df[[paste0(Trait)]],list(df[[paste0(SNP)]]),FUN=mean, na.rm = TRUE)
  res3<-as.data.frame(t(res[2]))
  ave<-mean(df[[paste0(Trait)]],na.rm=TRUE)
  names(res3)<-res$Group.1
  #val<-log10((res3$Ref)/(res3$Alt))
  val<-(res3$Ref)/(res3$Alt)
  #val<-(7.12-res3$Ref)/(7.12-res3$Alt)
  #val<-abs(res3$Ref-res3$Alt)
  Alt2<-subset(df,df[[paste0(SNP)]]=="Alt")
  Ref2<-subset(df,df[[paste0(SNP)]]=="Ref")
  #fitG<-glm(paste0(Trait)~paste0(SNP)+PC1+PC2+PC3,data=df)
  #p<-coef(summary(fitG))[2,4]
  p="NA"
  maf1<-(2*sum(df[[paste0(SNP)]]=="Alt")+sum(df[[paste0(SNP)]]=="Het"))/656
  if(maf1<0.5){
    maf<-maf1
  } else {
    maf<-(1-maf1)
  }
  out<-as.data.frame(cbind(val,p,maf))
  return(out)
}


outF<-as.data.frame(cbind(val,p,maf))

for (i in SNPsList$SNPs){ 
  outF<-rbind(outF,CalcEfx(i,"S2_Butyrate"))
}


 outF<-outF[-1,]
#rownames(outF<-SNPsList$SNPs)
outF$pos<-SNPs1$POS
outF1<-subset(outF,maf>=0.2)

#ggplot(data=outF1,aes(x=pos,y=val,color=maf))+geom_point()+
#  theme_classic()+theme(text=element_text(size=24,family = "sans"))+geom_smooth(colour="red")+ geom_hline(yintercept=0, linetype="dashed",color="black")


#Pos and Neg:
outPos<-subset(outF1,outF1$val>0)
outNeg<-subset(outF1,outF1$val<0)

ggplot(data=outF1,aes(x=pos,y=val,color=maf))+geom_point()+ theme_classic()+theme(text=element_text(size=24,family = "sans"))+
  geom_hline(yintercept=1, linetype="dashed",color="black")

#Clustering?
df3<-subset(df3,S1_PreIndBT!="#VALUE!")

#Pop<-read.csv("./SAP.kinship.csv",header = TRUE)
#Pop<-as.matrix(df3[459:468])
#row.names(Pop)<-Pop$X
#Pop$X<-NULL
#SAPDist<-dist(Pop,"euclidean")

#SAPH<-hclust(d=dist(Pop))




df3$S1_PreIndBTR<-rank(df3$S1_PreIndBT)
df3$S2_PreIndBTR<-rank(df3$S2_PreIndBT)

SAPD<-cbind(df3$S1_PreIndBTR,df3$S2_PreIndBTR)
SAPD<-matrix(as.numeric(SAPD),ncol=2)
row.names(SAPD)<-df3$Line
colnames(SAPD)<-c("S1_PrebioticIndex","S2_PrebioticIndex")

SAPDist<-dist(as.data.frame(SAPD),"euclidean")
SAPH<-hclust(SAPDist,method = "complete")

svg("PreIndHeatUnclust.svg",width = 10,height=20)
heatmap(SAPD,scale="column", Colv=NA,Rowv = as.dendrogram(SAPH),labRow = "",labCol = "",col= colorRampPalette(brewer.pal(8, "Reds"))(25))
dev.off()

legend(x="bottomright", legend=c("min", "mid", "max"), 
       fill=colorRampPalette(brewer.pal(8, "Reds"))(3))



heatmap(SAPD,scale="column", Colv = NA,labRow = "",labCol = "",col= colorRampPalette(brewer.pal(8, "Blues"))(25))


###################################
#Plot effect manhattans:

PlotPrepMM<-function(Trait){
  Out1<-as.data.frame(OutM$pos/1000000)
  names(Out1)<-"pos"
  Out1$effect<-OutM[[paste0(Trait,"_Effect")]]
  Out1$p<--log10(OutM[[paste0(Trait,"_V2.MLM")]])
  #Out1$effect<-abs(Out1$effect)
  ggplot(data=Out1,aes(x=pos,y=effect,color=p))+geom_point()+ theme_classic()+theme(text=element_text(size=20,family = "sans"))+
    geom_hline(yintercept=0, linetype="dashed",color="black")+ylab("Effect")+xlab("Position (MB)")+ labs(color = bquote(-log[10](P)))+
    ggtitle(paste0(Trait))+ scale_y_continuous(expand = c(0, 0.15))+scale_color_gradient(low = "azure4", high = "firebrick2") 
}


setwd("MEL6A_effect/")

list1<-list.files("./",pattern=".csv")

OutF<-read.csv(list1[4])
OutF$p<-(log10(OutF$p))*-1

OutF$pos<-OutF$pos/1000000
OutF<-OutF[order(OutF$p),]

ggplot(data=OutF,aes(x=pos,y=val,color=p))+geom_point()+ theme_classic()+theme(text=element_text(size=20,family = "sans"))+
  geom_hline(yintercept=1, linetype="dashed",color="black")+ylab("Effect")+xlab("Position (MB)")+ labs(color = bquote(-log[10](P)))+
  scale_y_continuous(expand = c(0, 0.15))+scale_color_gradient(low = "azure4", high = "firebrick2") 


colnames(OutF)<-sub("*.csv*","",colnames(OutF))
list2<-sub(".csv*","",list1)

##
thresh<--log10(0.05/nrow(OutF))

#Save plots:
for (i in list2){
  ggsave(filename=paste0(i,".png"),plot=print(PlotPrepEF(i),scale=1))
}
for (i in list2){
  ggsave(filename=paste0(i,"_efx.png"),plot=print(PlotPrepMF(i),scale=1))
}


#Make a list of all fastq files
setwd("~/SAP/SAP_16S_Data/")
files<-as.data.frame(list.files(path = "S1_1_4/", pattern = NULL, all.files = TRUE,
                  full.names = FALSE, recursive = FALSE,
                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE))
files<-cbind(files,as.data.frame(list.files(path = "S1_5_8/", pattern = NULL, all.files = TRUE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)))
files<-cbind(files,list.files(path = "S1_9_12/", pattern = NULL, all.files = TRUE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE))
files<-cbind(files,list.files(path = "S2_1_4/", pattern = NULL, all.files = TRUE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE))
files<-cbind(files,list.files(path = "S2_5_8/", pattern = NULL, all.files = TRUE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE))
files<-cbind(files,list.files(path = "S2_9_12/", pattern = NULL, all.files = TRUE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE))

write.csv(files,"fastqs.csv")

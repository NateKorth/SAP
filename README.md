## Sorghum Diversity Panel Microbiome Project
Code from the GWAS analysis of microbiome phenotypes. Most code is in R, more annotations to come/by request (nate.korthATgmail.com)

All phenotype data referenced in this code availible in publication (LINK TO COME)
All sorghum genotype data are attached to this publication: https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.15853 and are available here: https://www.ebi.ac.uk/ena/browser/view/PRJEB51985
Publically available sorghum phenotypes used in this study were compiled by Mural et. al.:  https://pubmed.ncbi.nlm.nih.gov/34100945/ 

Most of the following code is done independently by subject (human microbiome)

## STEP 1 ASV table generation
Apply the script ASV_PickingNK.sh which uses the qiime2 pipeline to assign ASVs and taxonomy using the SILVA database (At this point, this is an old build of the database, unless you're trying to replicate this work I'd recommend a newer version)
The raw sequence data is availible at: NCBI SRA database as project accession PRJNA1012736
Metadata required for analysis is avilable as supplementary table 2 in publication (LINK TO COME)

##STEP 2 Normalize (CSS) ASV table, and identify and remove outliers
Apply script SAP_Normalize_Plot_ASVlvl.R which takes qiime output and generates a normalized counts table and relative abundance table for later use, this script also includes step one in the three step outlier detection (see methods)
Apply script RemoveOutliersPrep.R for outlier detection steps 2 and 3

## STEP3 Heritablity and BLUE calculation
Apply the scrip BLUPs_h2_SAP.R -> this scrip is designed to first calculate broad sense heritability and filter out taxa with low values for h2. Next it calculates BLUEs - Best linear unbiased estimators - for each taxa for use in GWAS
Both methods rely on linear mixed models in sommer. In heritability calculations the genotype is set as a random effect / for BLUE calculations it is set as a fixed effect.

#STEP 4 Prep Genotype Data for GWAS
Use VCF/BCF tools to filter vcf file (SNP set)
```
Generate LinesInPheno.csv in R:
head -1 SAP_imputed.hmp > hmpHeader.tsv
ml R/4.0
R

df1<-read.csv("../input/CompiledBLUEs.csv")

#Remove any lines without a PI identifier
df2<-subset(df1,PI!="NA")

#Import the list of of genotypes:
LinesInHmp<-read.delim("hmpHeader.tsv",header = FALSE)
LinesInHmp2<-as.data.frame(t(LinesInHmp[12:369]))
names(LinesInHmp2)<-"PI"

#While we're thinking about lines and their names...
#If you peak at the data sets you see the nominclature in the genotype file contains underscores while the phenotype file does not (Typical)
#Lets just add those underscores in the phenotype file:

df2$PI<-sub("PI","PI_",df2$PI)

#Remove lines from the phenotype file that aren't in the genotype:
df3<-df2[which(df2$PI %in% LinesInHmp2$PI),]
```

```
ml vcftools
vcftools --vcf SAP_imputed.vcf --out SAP_imputed_Filter1 --keep LinesInPheno.csv --recode

vcftools --vcf SAP_imputed_Filter1.recode.vcf --out SAP_imputed_Filter2 --maf 0.05 --recode

ml bcftools
bcftools filter SAP_imputed_Filter2.recode.vcf --exclude 'F_PASS(GT=="het") > 0.1' -o SAP_imputed_Filter3.vcf
```

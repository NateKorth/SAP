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

## STEP 2 Normalize (CSS) ASV table, and identify and remove outliers
Apply script SAP_Normalize_Plot_ASVlvl.R which takes qiime output and generates a normalized counts table and relative abundance table for later use, this script also includes step one in the three step outlier detection (see methods)
Apply script RemoveOutliersPrep.R for outlier detection steps 2 and 3

## STEP 3 Heritablity and BLUE calculation
Apply the scrip BLUPs_h2_SAP.R -> this scrip is designed to first calculate broad sense heritability and filter out taxa with low values for h2. Next it calculates BLUEs - Best linear unbiased estimators - for each taxa for use in GWAS
Both methods rely on linear mixed models in sommer. In heritability calculations the genotype is set as a random effect / for BLUE calculations it is set as a fixed effect.

## STEP 4 Prep Genotype Data for GWAS
Use R to get list of genotypes in phenotype and SNP files (use to trim both files)
```
#Generate LinesInPheno.csv in R:
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
Use VCF/BCF tools to filter vcf file (SNP set)
```
ml vcftools
vcftools --vcf SAP_imputed.vcf --out SAP_imputed_Filter1 --keep LinesInPheno.csv --recode

vcftools --vcf SAP_imputed_Filter1.recode.vcf --out SAP_imputed_Filter2 --maf 0.05 --recode

ml bcftools
bcftools filter SAP_imputed_Filter2.recode.vcf --exclude 'F_PASS(GT=="het") > 0.1' -o SAP_imputed_Filter3.vcf
```

## STEP 5 Conduct GWAS
For a general GWAS tutorial see https://github.com/NateKorth/GWASTutorial
Agrinomic phenotypes are available as supplemental table 3 in publication
BLUEs are available as supplemental dataset 1 in publication.

Using all phenotypes and trimmed SNP file, apply RunRMIPLoop.R to run RMIP GWAS for every trait in phenotype file

## STEP 6 Plot GWAS output
Here use PlotBigGWAS_LD.R to subset and annotate GWAS output, trim non significant hits, and generated many kinds of plots. Admittedly this code can be difficult to follow as some of the plots were not included in the publication and some lines of code were not completed. But within these 1200 lines of code are some very useful tools for plotting multiple GWAS results. As well as digging underneath a specific GWAS peak. If you're unable to get what you need from this code please let me know and I am happy to help 

For MEL dissection (applied within the previous script) a list of sorghum genes within each MEL was pulled using the ExtractGenes.R script

## STEP 7 Analyze qPCR data
Use Analyze_qPCR.R to convert qPCR CT values to log(CFU/mL) based on Faecalibacterium standard curve, do statistical analysis, and generate heatmap.
The qPCR output of Faecalibacterium including metadata for the MEL6A validation study in the publication are available as Supplemental Dataset 3


For autoencoder platform see: https://github.com/mtross2/autoencoder_hyperspec_ref
If you'd like to apply any of this code to your pipeline and are meeting difficulties, feel free to reach out - Since writing this code I've developed a few more elegant solutions to answer these and similar questions.

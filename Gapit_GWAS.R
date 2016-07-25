#################load phenotype file
drgphenos<-read.table("~/Google Drive/GWAS_lydia/drgphenos_NRCRI_GAPIT_traits.txt",stringsAsFactors=F, header = TRUE)
#################load the dosage file without reduction by maf############################################################
setwd("~/Google Drive/MultikernelV6/")
load("V6_NEXTGEN_Dosages_50115.Rdata")
table(rownames(snps)%in%drgphenos$taxa)
#########extract phenotypes with clone names in snp file
test<-drgphenos[!(drgphenos$taxa%in%rownames(snps)),]
#########extract snps with clone names in phenos file
merge<-snps[rownames(snps)%in%drgphenos$CLONE,]
merge[1:10,1:10]
M<-merge-1#A.mat requires values 1,0,-1 values
M[1:10,1:10]
M<-round(M)
#################I use the A.mat function of rrBLUP to impute if missing values if not skip this step
library(rrBLUP)
A1 <- A.mat(M,impute.method="mean",n.core=2,return.imputed=T,max.missing=NULL,min.MAF=NULL)
snps_imp<-A1$imputed
dim(snps_imp)
snps_imp[1:10,1:10]
snps<-snps_imp+1 #GAPIT requires 0,1,2 values
snps[1:10,1:10]
snps<-round(snps)###############round because imputing gives decimal values
save(snps,file="500_imp_Lydia.Rdata")

####################################################################
Now that formating is done we can do GWAS with GAPIT !
####################################################################

###############load phenotypes##########################################################################################
phenos_1<-read.table("~/Google Drive/GWAS_lydia/drgphenos_NRCRI_GAPIT_traits.txt",stringsAsFactors=F, header = TRUE)
phenos_1<-phenos_1[!duplicated(phenos_1$taxa),]# should not change but will do if some names were duplicated
###################################subset dosage file by accessions in pheno file#######################################
table(rownames(snps)%in%phenos_1$taxa)#check how many names match
gwassnps<-snps[rownames(snps)%in%phenos_1$taxa,]
gwassnps<-as.data.frame(gwassnps,stringsAsFactors=F)
class(gwassnps)
rownames(gwassnps)
gwassnps<-as.data.frame(gwassnps,stringsAsFactors=F)
RPed <- sapply(gwassnps, as.integer) 
RPed[1:10,1:10]
RPed<-as.data.frame(RPed,stringsAsFactors=F)
rownames(RPed)<-rownames(gwassnps)
RPed$taxa<-rownames(RPed)###############format of file to meet GAPIT requirements
RPed[1:10,1:10]
RPed<-RPed[,c(ncol(RPed),1:(ncol(RPed)-1))]
rownames(RPed) <- NULL
RPed[1:10,1:10]
####################################subset phenos based on marker file##################################################
phenos_1<-phenos_1[phenos_1$taxa%in%RPed$taxa,]
rownames(phenos_1) <- NULL
###################################map file of markers after imputation#################################################
map<-colnames(RPed)#################column names usually are i.e S1_3456.. this is the chr and position we can extract for the map
map<-as.matrix(map)
write.table(map,file="map_gwas_V6_impute_NRCRI_500.txt")#open in excel and format if needed
map<-read.table("~/Google Drive/GWAS_lydia/map_gwas_V6_impute_NRCRI_500.txt",header=T,stringsAsFactor=F)
map$Chromosome=as.numeric(map$Chromosome)
map$Position=as.numeric(map$Position)
########################################################################################################################
#run libraries  
  library('MASS') # required for ginv
  library(compiler) #required for cmpfun
  library("scatterplot3d")
  library(multtest)
  library(gplots)
  library(LDheatmap)
  library(genetics)
  library(compiler) #this library is already installed in R library("scatterplot3d")
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt") 
#rename files
  myY  <- phenos_1
  myGD <- RPed
  myGM <- map
#Run GAPIT check manual to see the different settings that apply
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    PCA.total=4,
    kinship.cluster = c("ward"),
    SNP.MAF = 0.05)
  
  
 

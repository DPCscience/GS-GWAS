load("~/Google Drive/CBSDrob-dunia/Dosage_uganda_1066.Rdata")
dups<-dups[!duplicated(rownames(dups)),]
snps<-dups
rownames(snps)<-gsub("Ug","UG",rownames(snps))
snps<-snps[,grep("S11",colnames(snps))]
snps<-round(snps)

###the SNP that will be the reference SNP is the one with the highest lod score
res<-read.table(file="~/Google Drive/CBSDrob-dunia/results_TP/TP_CBSD_3s.loco.mlma",header=TRUE)
gwasmrkRR<-data.frame(-(log10(res$p)),stringsAsFactors=FALSE) 
colnames(gwasmrkRR)<-"logP"
comp<-cbind(res,gwasmrkRR)
interval=comp[grep("S11",comp$SNP),] ###subset by chromosome
max(interval$logP, na.rm = TRUE)
grep("10.1805",interval$logP)
interval[1350:1352,]  #S11_23277578
interval2=interval[(interval$bp > 23077578 & interval$bp <= 23477578),] ###snps in a window of 200kb for ld calculation
snp_name<-interval2$SNP
snp_name<-as.character(snp_name)
snpsed<-snps[,colnames(snps)%in%snp_name]#subset snps for ld calculation
snps<-snpsed

#####convert snps to plink files for ld calculation
snps<-snps-1
FUN<-function(x){
  y<-x
  y<-round(y) # Round to 0,1,2
  y<-gsub("-1","G G",y,fixed=T)
  y<-gsub("0","A G",y,fixed=T)
  y<-gsub("1","A A",y,fixed=T)
  y[is.na(y)]<-"0 0"
  return(y)
}
RECODED<-apply(snps,2,FUN)
RECODED[1:10,1:10]
dim(RECODED)
PED<-data.frame(FID=1,IID=rownames(snps),PID=0,MID=0,SEX=0,PHENO=-9)
PED<-cbind(PED,RECODED)
PED[1:10,1:10]
dim(PED)
write.table(PED, file = "LD_CHR11.ped", sep = " ", quote = FALSE, col.names=FALSE, row.names=FALSE)

# PLINK: MAP
markers<-colnames(snps)
markers1<-gsub("S","",markers)
positions<-data.frame(matrix(as.numeric(unlist(strsplit(markers1,"_"))),byrow=T,ncol=2))
map<-data.frame(Chr=positions$X1,Marker=colnames(snps), Dist=0, Pos=positions$X2)
head(map)
write.table(map, file = "LD_CHR11.map", sep = " ", quote = FALSE, col.names=FALSE, row.names=FALSE)

#####run plink different options to run use the reference snp with the highest lod score in GWAS with the option --ld-snp
system("~/Desktop/plink_mac/./plink",intern=TRUE)
system("~/Desktop/plink_mac/./plink --file LD_CHR11 --make-bed --out LD_CHR11",intern=TRUE)
system("~/Desktop/plink_mac/./plink --bfile LD_CHR11 --r2 --ld-snp S11_23277578 --chr 11 --out LD_CHR11")
system("~/Desktop/plink_mac/./plink --bfile LD_CHR11  --r2 --ld-window-r2 0.0000000000000001 --ld-snp S11_23277578 --chr 11 --out LD_CHR11")

ld<-read.table(file="LD_CHR11.ld",header=TRUE)
#ld=ld[(ld$SNP_A == "S11_23277578" | ld$SNP_B =="S11_23277578"),] 
ld_res<-data.frame(ld$SNP_B,ld$R2) #subset the ld output
merged<-merge(interval2,ld_res,by.x="SNP",by.y="ld.SNP_B",all=TRUE) #merge with snps in the interval window
merged[is.na(merged)] <- 0 # plink only gives you an output with r2<0.2 so all the other values for snps can be replaced by 0


####plot
plot<-data.frame(ID=rep(0,52),color=rep(0,52))
plot[,1]<-merged$SNP
inter1<-as.character(merged[merged$ld.R2 < 0.2,"SNP"])
inter2<-as.character(merged[(merged$ld.R2 >= 0.2 & merged$ld.R2 < 0.6),"SNP"])
inter<-as.character(merged[(merged$ld.R2 >= 0.6),"SNP"])

plot[plot$ID%in%inter1,2]<-1
plot[plot$ID%in%inter2,2]<-4
plot[plot$ID%in%inter,2]<-6
plot[plot$ID == "S11_23277578",2]<-20
merge2<-merge(merged,plot,by.x="SNP",by.y="ID",all=TRUE)

library("ggplot2")
library("ggthemr")
ggthemr('solarized')

scatter<-ggplot(data = new, aes(x = bp, y = logP,col=as.factor(color)))   + geom_point(size = 3,alpha=1)+ scale_y_continuous(breaks=seq(0, 11, 0.5))+theme(axis.text=element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16, face="bold"))
scatter

####create your own palette
color = c("#00008B","#000000", "#0000FF","#FFB90F" )
ugly <- define_palette(
  swatch = swatch,
  gradient = c(lower = "#FFB90F", upper ="#000000" ),
  background ="#F0FFFF"
)
ggthemr(ugly)

#####traditional ld plot
library(trio)
ld.out2 <- getLD(snpsed,which = c("rSquare"))
plot(ld.out2,y = "rSquare",cex = .5,alpha = 300,ciLD = c(0.7, 0.98))
class(res)
res[is.na(res)] <- 0
class(res)
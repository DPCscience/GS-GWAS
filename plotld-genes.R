
load("snps.Rdata")##
snps<-(snps+1)
snps[1:5,1:20]
snps<-snps[,grep("S11",colnames(snps))]
snps<-round(snps)
dim(snps)
###the SNP that will be the reference SNP is the one with the highest lod score
res<-read.table(file="~/Desktop/candidate genes/tp_gwas_6map.results.loco.mlma",header=TRUE)
gwasmrkRR<-data.frame(-(log10(res$p)),stringsAsFactors=FALSE) 
colnames(gwasmrkRR)<-"logP"
comp<-cbind(res,gwasmrkRR)
interval=comp[grep("S11",comp$SNP),] ###subset by chromosome
max(interval$logP, na.rm = TRUE)##find maximum lod score
grep("9.5292",interval$logP)
interval[1396:1400,]  #S11_23277578   S11_24197620   24197620   bp + >test>  bp - window 500kb
interval2=interval[(interval$bp > 23197620 & interval$bp <= 25197620),] ###snps in a window of 200kb for ld calculation
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
write.table(PED, file = "LD_CHR11_6M_TP_GWAS.ped", sep = " ", quote = FALSE, col.names=FALSE, row.names=FALSE)

# PLINK: MAP
markers<-colnames(snps)
markers1<-gsub("S","",markers)
positions<-data.frame(matrix(as.numeric(unlist(strsplit(markers1,"_"))),byrow=T,ncol=2))
map<-data.frame(Chr=positions$X1,Marker=colnames(snps), Dist=0, Pos=positions$X2)
head(map)
write.table(map, file = "LD_CHR11_6M_TP_GWAS.map", sep = " ", quote = FALSE, col.names=FALSE, row.names=FALSE)

#####run plink different options to run use the reference snp with the highest lod score in GWAS with the option --ld-snp
system("~/Desktop/plink_mac/./plink",intern=TRUE)
system("~/Desktop/plink_mac/./plink --file LD_CHR11_6M_TP_GWAS --make-bed --out LD_CHR11_6M_TP_GWAS",intern=TRUE)
system("~/Desktop/plink_mac/./plink --bfile LD_CHR11_6M_TP_GWAS --r2 --ld-snp S11_24197620 --chr 11 --out LD_CHR11_6M_TP_GWAS")
#system("~/Desktop/plink_mac/./plink --bfile LD_CHR11  --r2 --ld-window-r2 0.0000000000000001 --ld-snp S11_23277578 --chr 11 --out LD_CHR11")

ld<-read.table(file="LD_CHR11_6M_TP_GWAS.ld",header=TRUE)
#ld=ld[(ld$SNP_A == "S11_23277578" | ld$SNP_B =="S11_23277578"),] 
ld_res<-data.frame(ld$SNP_B,ld$R2) #subset the ld output
merged<-merge(interval2,ld_res,by.x="SNP",by.y="ld.SNP_B",all=TRUE) #merge with snps in the interval window
merged[is.na(merged)] <- 0 # plink only gives you an output with r2<0.2 so all the other values for snps can be replaced by 0


####plot
plot<-data.frame(ID=rep(0,165),color=rep(0,165))
plot[,1]<-merged$SNP
inter1<-as.character(merged[merged$ld.R2 < 0.2,"SNP"])
inter2<-as.character(merged[(merged$ld.R2 >= 0.2 & merged$ld.R2 < 0.6),"SNP"])
inter<-as.character(merged[(merged$ld.R2 >= 0.6),"SNP"])

plot[plot$ID%in%inter1,2]<-"R2 < 0.2"
plot[plot$ID%in%inter2,2]<-"0.6>R2 >0.2"
plot[plot$ID%in%inter,2]<-"r2 >= 0.6"
plot[plot$ID == "S11_24197620",2]<-"ref SNP"
merge2<-merge(merged,plot,by.x="SNP",by.y="ID",all=TRUE)

library("ggplot2")
library("ggthemr")
ggthemr('solarized')

scatter<-ggplot(data = merge2, aes(x = bp, y = logP,col=as.factor(color)))   + geom_point(size = 3,alpha=1)+ 
  scale_y_continuous(breaks=seq(0, 11, 0.5))+theme_bw()+
  theme(axis.text.x=element_blank(), axis.title=element_text(size=14), plot.title=element_text(size=16, face="bold"))+
  theme(axis.title.x=element_blank(),axis.text.y = element_text(size=15, color="black"),axis.title.y=element_text(size=18,face="bold", color="black"))+
  xlim(23200000, 25250000)+
  scale_color_manual(values = c("#d1ceaa","#2b8e86",'#E69F00',"red"),labels = c(expression(r^"2"*"<0.2"), expression("0.6>"*r^"2"* ">0.2"),expression(r^"2"*">0.6"),"reference SNP hit")) +
  geom_hline(yintercept = 5.8,colour ="darkgrey", size =1.3)+ labs(y=(expression(-log[10]*"P value")))+
  theme(plot.margin = unit(c(1,2,0.05,3), "cm"))+
  theme(legend.background = element_rect(fill = "white"),legend.position = "top",legend.direction="horizontal",legend.text=element_text(size=14))+ 
  theme(legend.title=element_blank())

####combine ld plot with gene positions, use the arrowplot code
library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(scatter), ggplotGrob(p), size = "first"))


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
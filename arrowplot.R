## the ggplot package creates the graph
## geom_line is used to create each gene line
require(ggplot2)
library(directlabels)
gff=read.table("~/Google Drive/NaCCRI-GWASpop/test2.txt", header=T, stringsAsFactor=F)#label,name,start(gene start position)
p<-ggplot(gff, aes(x=start, y=as.factor(label))) + 
geom_point(size=0.0005) +
geom_line(arrow = arrow(length=unit(0.3,"cm"), ends="first", type = "closed"), size = 1)+
geom_dl(aes(label = name),method = list(dl.combine("first.points"), cex = 0.8))+
theme(axis.title.y=element_blank(),axis.title.x=element_text(size=16,face="bold"),axis.text.y = element_text(size=13))+
theme(aspect.ratio = 0.2)+ theme(legend.position="none")+
xlim(23200000, 25250000)+ labs(x = "Position on chromosome 11 (bp)")+theme(axis.title.x=element_text(margin=margin(20,0,0,0)))


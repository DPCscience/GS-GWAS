
rm(list = ls())
setwd("~/Google Drive/MultikernelV6/")
#######
load("IITAV6.Rdata")
snps <- snps[!duplicated(rownames(snps)),]
######reading files and matching to snp dosage files
drgphenos8=read.table("iitanewnamesDosage.txt", header=T, stringsAsFactor=F)
drgphenos2=read.table("namesdrgphenos.txt", header=T, stringsAsFactor=F)
drgphenos5<-merge(drgphenos8,drgphenos2,by="CLONE",ALL=FALSE)
drgphenos5<-drgphenos5[which(drgphenos5$CLONE%in%row.names(snps)),]


################# STARCH #########################################################
starch<-read.table(file="annotation/starch.snps",header=TRUE,stringsAsFactor=F)
starch<-starch[!duplicated(starch[,1]),]
starch<-as.data.frame(starch)
colnames(starch)<-"SNP"

root<-read.table(file="annotation/roots.snps",header=TRUE,stringsAsFactor=F)
root<-root[!duplicated(root[,1]),]
root<-as.data.frame(root,stringsAsFactor=F)
colnames(root)<-"SNP"


#######
snps_ind<-snps[rownames(snps) %in% drgphenos5$CLONE,]
snpsind<-rownames(snps_ind)
#######
starch<-read.table(file="annotation/starch.snps",header=TRUE,stringsAsFactor=F)
starch<-starch[!duplicated(starch[,1]),]
starch<-as.data.frame(starch,stringsAsFactor=F)
colnames(starch)<-"SNP"
#######
snps_col<-snps[,colnames(snps) %in% root$SNP]
snpscol<-colnames(snps_col)
#######
snpsed<-snps[snpsind,snpscol]
dim(snpsed)
########################################################
dataset_ed=merge(drgphenos5,snpsed,by.x="CLONE",by.y="row.names")
rownames(dataset_ed)<-dataset_ed$CLONE
dim(dataset_ed)
traits<-dataset_ed[,1:28]
predictors<-dataset_ed[,29:949]
dataset<-dataset_ed
#dataset<-cbind(traits,markers)
#trait=c("DM","logRTNO","logRTWT","logSHTWT")
#genoID="CLONE"
#nFolds=5
#nRepeats=25

############
library(foreach)
library(doParallel)
trait=c("DM","logRTNO","logRTWT","logSHTWT")
proctime<-proc.time()
cl<-makeCluster(2)
registerDoParallel(cl)
bio.RF.RTNO <- foreach(a=trait, bio=icount(), .inorder=TRUE) %dopar% {
  require(randomForest)
  traits<-c("DM","logRTNO","logRTWT","logSHTWT")
  crossval<-FoldCrossValidation.V1.RF(dataset,traits[bio],"CLONE",colnames(dataset_ed[,29:949]),5,25)
}

stopCluster(cl)
proc.time() - proctime
save(bio.RF.RTNO,file="starch_snps_logRTNO.Rdata")

################

DM_RT_RF<-FoldCrossValidation.V1.RF(dataset,trait, genoID,predictors,nFolds,nRepeats)
     

# Fold cross-validation for randomForest V1.0 --------------------------------------------------------------------
FoldCrossValidation.V1.RF <- function(dataset, trait, genoID, predictors, nFolds, nRepeats){
    require(randomForest)  
    data<-dataset  # [which(dataset[ ,genoID] %in% rownames(Klist[[1]])),]
    data<-data[!is.na(data[,trait]),]    
    # rownames(data)<-data$CLONE
    nInd <- dim(data)[1] 
    accuracies<-data.frame()
    blups<-data.frame()
    for (rep in 1:nRepeats){ 
        print(paste("Rep ",rep," of ",nRepeats,sep=""))
        folds <- sample(rep(1:nFolds, length.out=nInd))
        BLUPSthisRep<-data.frame()
        for (fold in 1:nFolds){
            print(paste("Fold ",fold," of ",nFolds,sep=""))
            indInFold <- which(folds == fold)
            indNotInFold <- which(folds != fold)
            ValidSet<-data[indInFold,genoID]
            TrainSet<-data[indNotInFold,genoID]
            out=randomForest(y = data[TrainSet,trait], x = data[TrainSet,predictors], xtest=data[ValidSet,predictors], data=data, ntree=500,mtry=300,importance=TRUE)
            IMPthisrep<-data.frame(IMP=out[["test"]]$importance)
            BLUPSthisFold<-data.frame(CLONE=names(out[["test"]]$predicted), PRED=out[["test"]]$predicted)
            BLUPSthisRep<-rbind(BLUPSthisRep,BLUPSthisFold)
        } 
        BLUPSthisRep[,"Rep"]<-rep
        # Calc accuracy after predicting each fold to complete the rep
        BLUPSthisRep<-merge(BLUPSthisRep,data[,c("CLONE",paste(trait,".ebv",sep=""))],by="CLONE")
        accuracy.thisrep<-cor(BLUPSthisRep$PRED,BLUPSthisRep[,paste(trait,".ebv",sep="")], use="complete.obs")
        AccuracyThisRep<-data.frame(Trait=trait,Rep=rep,Accuracy=accuracy.thisrep)
        accuracies<-rbind(accuracies,AccuracyThisRep) 
        blups<-rbind(blups,BLUPSthisRep)
    }
    return(list(accuracies=accuracies,blups=blups,importance=IMPthisrep))
}


imp<-round(importance(out), 2)

################ ROOT ############################################################
root<-read.table(file="annotation/roots.snps",header=TRUE,stringsAsFactor=F)
#root<-root[!duplicated(root[c("SNP")]),]
root<-root[!duplicated(root[,1]),]
root<-as.data.frame(root)
colnames(root)<-"SNP"
root<-root[!(root$SNP%in% starch$SNP), ]
root<-as.data.frame(root)
colnames(root)<-"SNP"

############## MRNA ###############################################################

mrna<-read.table(file="annotation/mirnas.snps",header=TRUE,stringsAsFactor=F)
mrna<-mrna[!duplicated(mrna[,1]),]
mrna<-as.data.frame(mrna)
colnames(mrna)<-"SNP"

tf<-read.table(file="annotation/TFs.snps",header=TRUE,stringsAsFactor=F)
tf<-tf[!duplicated(tf[,1]),]
tf<-as.data.frame(tf)
colnames(tf)<-"SNP"

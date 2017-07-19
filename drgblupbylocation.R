pheno2<-read.csv("/Users/duniapinodelcarpio//Downloads/Kigumba_GxE2015A_Final.csv",header = T,stringsAsFactors = F)
pheno1<-read.csv("/Users/duniapinodelcarpio//Downloads/Arua_GxE2015A_Final.csv",header = T,stringsAsFactors = F)

####Be aware of Clone column name, the model you specify and the name of the traits
library(foreach)
library(doParallel)
########### This is the source code, 
myrep<-function(traits,pheno){
  dfList <- list() 
  for(i in 1:length(traits))
  {
    require(lme4)
    model<-lmer(data=pheno,formula = get(traits[i]) ~  (1|Clones) ) #changes depending on the fixed and random effects
    deregress<-function(model, trait){##deregress code Uche Okeke
      BLUP <- ranef(model, condVar=TRUE)$Clones###CLONE column name
      PEV <- c(attr(BLUP, "postVar"))
      clone.var <- c(VarCorr(model)$Clones)###CLONE column name
      data<-as.data.frame(VarCorr(model))
      out <- BLUP/(1-(PEV/clone.var))#deregressed BLUPs
      r2 <- 1-(PEV/clone.var)
      h2 = clone.var/(sum(data$vcov))#change according to the random effects
      wt = (1-h2)/((0.1 + (1-r2)/r2)*h2)
      location=pheno$Location
      return(list(Trait=trait, drgBLUP=out, BLUP=BLUP, weights=wt,H2=h2,Location=location))
    }
    drg<-deregress(model,traits[i])#model is the lmer object
    dfList[[i]]<-drg
  }
  return(dfList)
}

########drgBLUPs output

listloc<-list(pheno1,pheno2)
traits<-c( "CBSDs","CBSDi","RootWt","HI", "DMC")
resList <- list() 
for(i in 1:length(listloc))
{
  traits<-c( "CBSDs","CBSDi","RootWt","HI", "DMC")
  res<-myrep(traits,listloc[[i]])
  resList[[i]]<-res
}

#########Heritability output
heritability2<- data.frame()
for(i in 1:length(resList)) {
  for (j in 1:length(traits)) {
    drg<-data.frame(H2 = (resList[[i]][[j]]$H2),trait=(resList[[i]][[j]]$Trait),Location=unique(resList[[i]][[j]]$Location))
    heritability2<-rbind(heritability2,drg)
  }
}

#########drgBLUPs output formatting
#i=number of location files you input
#j=number of traits
#########formatting
drgList <- list() 
for(i in 1:length(listloc)) {
  drgphenos<-data.frame(CLONE = unique(listloc[[i]]$Clones),stringsAsFactors=F)
  for (j in 1:length(traits)) {
    drg<-data.frame(CLONE = rownames(resList[[i]][[j]]$drgBLUP),stringsAsFactors=F)
    drg[,resList[[i]][[j]]$Trait]<-resList[[i]]$drgBLUP
    drg[,paste(resList[[i]][[j]]$Trait,"ebv",sep=".")]<-resList[[i]][[j]]$BLUP
    drg[,paste(resList[[i]][[j]]$Trait,"wt",sep=".")]<-resList[[i]][[j]]$weights
    #drg[,("Location")]<-resList[[i]][[i]]$Location
    drgphenos<-merge(drgphenos,drg,by="CLONE",all=T)
  }
  drgList [[i]]<-drgphenos
  for(i in 1:length(drgList)) {
    names(drgList)[i]<-unique(listloc[[i]]$Location)
  }
}


list2env(drgList, environment())

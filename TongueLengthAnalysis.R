library(MCMCglmm)
library(prevalence)
library(reshape)

data<-read.csv("Data.csv",header=T)
data$SpatialComposition<-as.factor(data$SpatialComposition)
reshapeddata<-reshape(data, direction="long", varying=list(names(data)[-c(1:4)]), times=names(data)[-c(1:4)])
colnames(reshapeddata)<-c("PoolID","Host","Count","SpatialComposition","Virus","Status","ID")

prevs<-vector("list",length(unique(reshapeddata$Virus))*length(unique(reshapeddata$Host)))

mats<-array(NA,dim=c(length(unique(reshapeddata$Host)),length(unique(reshapeddata$Virus)),1000))
rownames(mats)<-unique(reshapeddata$Host)
colnames(mats)<-unique(reshapeddata$Virus)

###generate posterior prevalences -- functionally a method of sampling from the likelihood of the data, and samples from likelihood in distance matrix
###when one pool prevalence as written does not work, code extracted from r package and made to run in the positive case, and analytical solution used in the negative case

k<-1
for (i in 1:length(unique(reshapeddata$Virus))) {
  for (j in 1:length(unique(reshapeddata$Host))) {
    workinglist<-vector("list",3)
    temp<-subset(reshapeddata,Virus==unique(reshapeddata$Virus)[i]&Host==unique(reshapeddata$Host)[j])
    if (nrow(temp)!=1) {
      workinglist[[1]]<-as.character(unique(reshapeddata$Host)[j])
      workinglist[[2]]<-as.character(unique(reshapeddata$Virus)[i])
      workinglist[[3]]<-do.call("c",lapply(truePrevPools(temp$Status,temp$Count,nchains=4)@mcmc$TP, function(x) {return(as.numeric(x))}))[seq(1,40000,by=10)]
      prevs[[k]]<-workinglist
      rm(workinglist)
    } else {
      if (temp$Status==0) {
        workinglist[[1]]<-as.character(unique(reshapeddata$Host)[j])
        workinglist[[2]]<-as.character(unique(reshapeddata$Virus)[i])
        workinglist[[3]]<-rbeta(4000,1,1+temp$Count)
        prevs[[k]]<-workinglist
        rm(workinglist)
      } else {
        data <- list(x = temp$Status, n = temp$Count, N = 1)
        mod<-jags.model(file = "~/Desktop/Pre-Glasgow/prevmodel.txt",
                   data = data,
                   n.chains = 4,
                   n.adapt = 1000)
        update(mod, n.iter = 10000, progress.bar = "none")
        samples <- coda.samples(mod, c("SE", "SP", "TP"), n.iter = 10000, thin = 10,
                                progress.bar = "none")
        workinglist[[1]]<-as.character(unique(reshapeddata$Host)[j])
        workinglist[[2]]<-as.character(unique(reshapeddata$Virus)[i])
        workinglist[[3]]<-c(as.numeric(samples[[1]][,3]),as.numeric(samples[[2]][,3]),as.numeric(samples[[3]][,3]),as.numeric(samples[[4]][,3]))
        prevs[[k]]<-workinglist
        rm(workinglist)
      }
    }
    mats[j,i,]<-sample(prevs[[j]][[3]],1000,replace=T)
    k<-k+1
  }
}

tonguelengths<-matrix(c(8.5,7.6,12.5,7.7,7.3,7.5,NA,8,6.4,7.5,8.5),nrow=length(row.names(mats)),ncol=1)
row.names(tonguelengths)<-row.names(mats)
tonguelengths<-tonguelengths[-7,]

tonguelengthdistances<-dist(tonguelengths)

correlationsraw<-rep(NA,1000)

for (i in 1:1000) {
  print(i)
  correlationsraw[i]<-mantel(tonguelengthdistances,dist(mats[-7,,i]),method="kendall")
}

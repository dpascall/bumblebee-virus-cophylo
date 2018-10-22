#!usr/bin/env Rscript

library(methods)
library(ape)
library(MCMCglmm)
library(rstan)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##read in trees

hosttree<-read.nexus("~/DoublePhyloStan/HostTreesReduced10.tree")
virustree<-read.nexus("~/DoublePhyloStan/VirusTreesReduced10.tree")

hosttreenumber<-length(hosttree)
virustreenumber<-length(virustree)

matrices<-hosttreenumber*virustreenumber

##read and reshape data

data<-read.csv("~/DoublePhyloStan/GeoNoHov.csv",header=T)
data$SpatialComposition<-as.factor(data$SpatialComposition)
reshapeddata<-reshape(data, direction="long", varying=list(names(data)[-c(1:4)]), times=names(data)[-c(1:4)])
colnames(reshapeddata)<-c("PoolID","Host","Count","SpatialComposition","Virus","Status","ID")
reshapeddata$VirusPhylogeny<-reshapeddata$Virus
reshapeddata$HostPhylogeny<-reshapeddata$Host
reshapeddata$HostConditional<-paste(reshapeddata$Host, reshapeddata$Virus, sep=".")
reshapeddata$VirusConditional<-paste(reshapeddata$Host, reshapeddata$Virus, sep=".")
reshapeddata$Interaction<-paste(reshapeddata$Host, reshapeddata$Virus, sep=".")

##generate matrices

virustrees<-array(NA,dim=c(matrices,length(unique(reshapeddata$Virus)),length(unique(reshapeddata$Virus))))
hosttrees<-array(NA,dim=c(matrices,length(unique(reshapeddata$Host)),length(unique(reshapeddata$Host))))

k<-1
for (i in 1:virustreenumber) {
  temp<-vcv(virustree[[i]],model="Brownian")
  for (j in 1:hosttreenumber) {
    virustrees[k,,]<-temp[order(rownames(temp)), order(rownames(temp))]
    k<-k+1
    print(k)
  }
}

k<-1
for (i in 1:virustreenumber) {
  for (j in 1:hosttreenumber) {
    temp<-vcv(hosttree[[j]],model="Brownian")
    hosttrees[k,,]<-temp[order(rownames(temp)), order(rownames(temp))]
    k<-k+1
    print(k)
  }
}

coevo<-array(NA,dim=c(matrices,length(unique(reshapeddata$Host))*length(unique(reshapeddata$Virus)),length(unique(reshapeddata$Host))*length(unique(reshapeddata$Virus))))
hostassem<-array(NA,dim=c(matrices,length(unique(reshapeddata$Host))*length(unique(reshapeddata$Virus)),length(unique(reshapeddata$Host))*length(unique(reshapeddata$Virus))))
virusassem<-array(NA,dim=c(matrices,length(unique(reshapeddata$Host))*length(unique(reshapeddata$Virus)),length(unique(reshapeddata$Host))*length(unique(reshapeddata$Virus))))

for (i in 1:matrices) {
  coevo[i,,]<-cov2cor(kronecker(virustrees[i,,], hosttrees[i,,]))
  hostassem[i,,]<-cov2cor(as.matrix(kronecker(hosttrees[i,,], Diagonal(nrow(virustrees[i,,])))))
  virusassem[i,,]<-cov2cor(as.matrix(kronecker(virustrees[i,,], Diagonal(nrow(hosttrees[i,,])))))
  virustrees[i,,]<-cov2cor(virustrees[i,,])
  hosttrees[i,,]<-cov2cor(hosttrees[i,,])
  print(i)
}

preddataset<-unique(reshapeddata[,c(2,5,12)])

dat<-list(datapoints=nrow(reshapeddata), virusspecies=length(unique(reshapeddata$Virus)),hostspecies=length(unique(reshapeddata$Host)),
          virushostcombinations=length(unique(reshapeddata$VirusConditional)),pools=length(unique(reshapeddata$PoolID)), 
          spatialcompositions=length(unique(reshapeddata$SpatialComposition)), matrices=matrices, y=reshapeddata$Status, 
          host=as.numeric(as.factor(reshapeddata$Host)),virus=as.numeric(as.factor(reshapeddata$Virus)),
          combination=as.numeric(as.factor(reshapeddata$Interaction)),pool=as.numeric(as.factor(reshapeddata$PoolID)),
          spatialcomposition=as.numeric(as.factor(reshapeddata$SpatialComposition)),samples=reshapeddata$Count,HostPhy=hosttrees,VirusPhy=virustrees,
          HostInter=hostassem,VirusInter=virusassem,CoevoInter=coevo, ID=c(1:nrow(reshapeddata)), 
          run_estimation=1, host_pred=as.numeric(as.factor(preddataset[,1])),virus_pred=as.numeric(as.factor(preddataset[,2])),
          combination_pred=as.numeric(as.factor(preddataset[,3])))

rm(data,hosttree,reshapeddata,temp,virustree,coevo,hostassem,hosttrees,i,virusassem,virustrees)

params<-c("sigma_virus", "sigma_virusphy", "sigma_host", "sigma_hostphy", "sigma_virusinter", "sigma_hostinter", "sigma_coevointer", "sigma_inter", "sigma_pool", 
          "sigma_spatial", "sigma_residual", "alpha", "ICC_virus", "ICC_virusphy", "ICC_host", "ICC_hostphy", "ICC_virusinter", "ICC_hostinter", "ICC_inter",
          "ICC_coevointer", "ICC_pool", "ICC_spatial", "ICC_residual", "ICC_nonphylogenetic", "ICC_phylogenetic", "denominator", "log_lik", "y_sim" , "y_pred", "lp__")

ptm <- proc.time()
fit <- stan("~/Exp1/DoublePhyloStan/CompleteCophylogeneticNC.stan", data = dat, iter = 3000, warmup = 2000, pars=params, thin = 1, chains=4, verbose = TRUE, init_r=0.01, refresh=1, control = list(adapt_delta = 0.9999, max_treedepth = 100))
proc.time() - ptm

save(fit, file="~/Exp1/DoublePhyloStan/Results/AllFull10.Rdata")

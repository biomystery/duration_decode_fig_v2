## load data 
rm(list = ls());
load(file='mRNA-Fit-avg-sp-v1c-pre.Rdata')#,list=c("genes.kdeg","genes.normCount",
source(file = "mRNA-Fit-avg-sp-v1c-funs.R")#rowMax, runModel.v1, runOptim.v1
library(gplots);library(RColorBrewer)
library(pheatmap)

for(N in 1:6){
  #N <- 2
  load(file = paste0('mRNA-Fit-avg-sp-v1c-run-N',N,'.Rdata'))#,list("v1.fit.MaxScale","lb","ub","initP"))
  ev <- new.env()
  load(file = paste0('mRNA-Fit-avg-sp-v1c-run-N',N,'_na.Rdata'),ev)
  # names  
  names(ev$v1.fit.MaxScale) <- lapply(ev$v1.fit.MaxScale,function(x) x$gene)
  names(v1.fit.MaxScale) <- lapply(v1.fit.MaxScale,function(x) x$gene)
  v1.fit.MaxScale[names(ev$v1.fit.MaxScale)] <- ev$v1.fit.MaxScale[names(ev$v1.fit.MaxScale)]
  tmp <- v1.fit.MaxScale[rownames(genes.normCount)]
  # model fit 
  v1.fit.1.r2 <- unlist(lapply(v1.fit.MaxScale, fun.calR2))
  for(i in 1:length(v1.fit.MaxScale)) v1.fit.MaxScale[[i]]$r2 <- v1.fit.1.r2[i]
  v1.fit.1.obj <- unlist(lapply(v1.fit.MaxScale, function(x) x$obj))  
  
  # plot the results 
  tmp.ord<- order(v1.fit.1.r2,decreasing = T)
  tmp.ord.2<- order(v1.fit.1.obj)
  v1.fit <- data.frame(
    r2 = v1.fit.1.r2,
    rmsd = v1.fit.1.obj
  )
  pdf(file=paste0("mRNA-Fit-avg-sp-v1c-run-N",N,".pdf"))
  par(mfrow=c(2,2))
  hist(v1.fit.1.r2)
  hist(v1.fit.1.obj)
  plot(v1.fit,main=paste("cor=",signif(cor(v1.fit.1.r2,v1.fit.1.obj,use = "na.or.complete"),2)))
  #identify(v1.fit$r2,v1.fit$rmsd,labels = rownames(v1.fit))
  abline(lm(v1.fit$rmsd~ v1.fit$r2),lwd=2,col=2) 
  scatter.smooth(x=v1.fit$r2,y=v1.fit$rmsd,lpars = list(col=2,lwd=4))
  
  #par(mfrow=c(3,3))
  #sapply(1:length(tmp.ord), function(x) fun.plotFit(v1.fit.MaxScale[[tmp.ord[x]]]))
  #par(mfrow=c(1,1))
  dev.off()
  save(file = paste0('mRNA-Fit-avg-sp-v1c-run-N',N,'_na.Rdata'),
       list=c("v1.fit.MaxScale","N"))
}











  




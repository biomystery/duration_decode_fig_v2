## load data 
rm(list = ls());
load(file='mRNA-Fit-avg-sp-v2-pre.Rdata')#,list=c("genes.kdeg","genes.normCount",
source(file = "mRNA-Fit-avg-sp-v2-funs.R")#rowMax, runModel.v1, runOptim.v1
library(gplots);library(RColorBrewer)
library(pheatmap)


for(N in 1:6){
  load(file = paste0('mRNA-Fit-avg-sp-v2-run-N',N,'.Rdata'))#,list("v1.fit.MaxScale","lb","ub","initP"))
  names(v1.fit.MaxScale) <- lapply(v1.fit.MaxScale,function(x) x$gene)
  eval(parse(text = paste0("r2.",N," <- unlist(lapply(v1.fit.MaxScale, fun.calR2))")))
}



r2.all <- data.frame(n1=r2.1,n2=r2.2,n3=r2.3,n4=r2.4,n5=r2.5,n6=r2.6) #,n4=r2.4,n5=r2.5,n6=r2.6
r2.final = data.frame(r2=rowMax(r2.all),
                      n= apply(r2.all,1, which.max),
                      delta_r2 = rowMax(r2.all) - r2.1)

r2.final$gene <- rownames(r2.final)
require(dplyr)
r2.final <- (arrange(r2.final,desc(r2)))
write.csv(file='results.csv',r2.final)


# load simulation data  ---------------------------------------------------

n.gene <- nrow(r2.final)
rownames(r2.final) <- r2.final$gene
getBestFit <-function(g=r2.final$gene[1]){
  g.res <- r2.final[g,]
  plot.tp_step <- 5 # 5mins step
  ev <- new.env()
  load(file=paste0('./mRNA-Fit-avg-sp-v2-run-N',g.res$n,'.Rdata'),ev)
  names(ev$v1.fit.MaxScale)<-sapply(1:n.gene,function(x) ev$v1.fit.MaxScale[[x]]$gene)
  fit.res <- ev$v1.fit.MaxScale[[g]]
  all.tps <- fit.res$bestFit.mt$LPS$time
  
  plot.tp_idx <- all.tps %in% seq(all.tps[1],all.tps[length(all.tps)],by=plot.tp_step)
  n.tps <- sum(plot.tp_idx)
  data.ts <- rbind(fit.res$bestFit.ctrl$LPS[plot.tp_idx,],
                   fit.res$bestFit.ctrl$TNF[plot.tp_idx,],
                   fit.res$bestFit.mt$LPS[plot.tp_idx,],
                   fit.res$bestFit.mt$TNF[plot.tp_idx,])
  data.ts$sti <- rep(c("lps","tnf","lps","tnf"),each=n.tps)
  data.ts$geno <- rep(c("wt","wt",'mt',"mt"),each=n.tps)
  data.ts$gene <- rep(g,n.tps*4)
  data.ts$mRNA[data.ts$geno =='wt']<- data.ts$mRNA[data.ts$geno=='wt']/fit.res$simScale.wt
  data.ts$mRNA[data.ts$geno =='mt']<- data.ts$mRNA[data.ts$geno=='mt']/fit.res$simScale.mt  
  rownames(data.ts) <- NULL
  data.ts
}


tmp<-sapply(1:n.gene, function(i) {print(i);getBestFit(rownames(r2.final)[i])},
            simplify = F);names(tmp) <- NULL ;tmp <- do.call('rbind',tmp)
tmp$time <- tmp$time/60

write.csv(file='bestFit.tc.csv',tmp,row.names = F)

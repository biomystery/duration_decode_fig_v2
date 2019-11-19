
rm(list = ls()); 
load(file='mRNA-Fit-avg-sp-v1c-pre.Rdata')#,list=c("genes.kdeg","genes.normCount",
source(file = "mRNA-Fit-avg-sp-v1c-funs.R")#rowMax, runModel.v1, runOptim.v1
# r2 compare  -------------------------------------------------------------
for(N in 1:6){
  load(file = paste0('mRNA-Fit-avg-sp-v1c-run-N',N,'.Rdata'))#,list("v1.fit.MaxScale","lb","ub","initP"))
  names(v1.fit.MaxScale) <- lapply(v1.fit.MaxScale,function(x) x$gene)
  eval(parse(text = paste0("r2.",N," <- unlist(lapply(v1.fit.MaxScale, fun.calR2))")))
  eval(parse(text = paste0("rmsd.",N," <- unlist(lapply(v1.fit.MaxScale, fun.calRMSD))")))
  
}

r2.all <- data.frame(n1=r2.1,n2=r2.2,n3=r2.3,n4=r2.4,n5=r2.5,n6=r2.6)
r2.all <- r2.all[complete.cases(r2.all),]
r2.final = data.frame(r2=rowMax(r2.all),
                      n= apply(r2.all,1, which.max))

r2.delta <- (r2.final$r2- r2.all$n1)
names(r2.delta) <- rownames(r2.final)
r2.final$delta <- r2.delta
r2.final$gene <- rownames(r2.final)
#require(dplyr)
#r2.final <- r2.final[order(r2.final[,1],decreasing = T),]


# parameters 
getBestFitPar <-function(g=r2.final$gene[1]){ # rmsd.final
  g.res <- r2.final[g,]
  ev <- new.env()
  load(file=paste0('./mRNA-Fit-avg-sp-v1c-run-N',g.res$n,'.Rdata'),ev)
  names(ev$v1.fit.MaxScale)<-sapply(1:length(ev$v1.fit.MaxScale),
                                    function(x) ev$v1.fit.MaxScale[[x]]$gene)
  df <- ev$v1.fit.MaxScale[[g]]$par
  names(df) <- c("kb","kt","Kd","k_deg","tau") 
  df
}
r2.final <- cbind(r2.final,t(sapply(rownames(r2.final), getBestFitPar))) 
write.csv(file='result.csv',r2.final)

# rmsd cmp ----------------------------------------------------------------

#for(N in 1:6){
#  load(file = paste0('mRNA-Fit-avg-sp-v1c-run-N',N,'.Rdata'))#,list("v1.fit.MaxScale","lb","ub","initP"))
#  names(v1.fit.MaxScale) <- lapply(v1.fit.MaxScale,function(x) x$gene)
#}

rmsd.all <- data.frame(n1=rmsd.1,n2=rmsd.2,n3=rmsd.3,n4=rmsd.4,n5=rmsd.5,n6=rmsd.6)
rmsd.final = data.frame(rmsd=rowMin(rmsd.all),
                      n= apply(rmsd.all,1, which.min))

rmsd.delta <- (rmsd.final$rmsd- rmsd.all$n1)
names(rmsd.delta) <- rownames(rmsd.final)
rmsd.final$delta <- rmsd.delta
rmsd.final$gene <- rownames(rmsd.final)
require(dplyr)
rmsd.final <- rmsd.final[order(rmsd.final[,1]),]
write.csv(file='result-RMSD.csv',arrange(rmsd.final,rmsd,rmsd.delta))



# load simulation data  ---------------------------------------------------
getBestFit <-function(g=r2.final$gene[1],r2.final=rmsd.final){
  g.res <- r2.final[g,]
  plot.tp_step <- 5 # 5mins step
  ev <- new.env()
  load(file=paste0('./mRNA-Fit-avg-sp-v1c-run-N',g.res$n,'.Rdata'),ev)
  names(ev$v1.fit.MaxScale)<-sapply(1:length(ev$v1.fit.MaxScale),function(x) ev$v1.fit.MaxScale[[x]]$gene)
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
  data.ts$mRNA[data.ts$geno == "wt"] <- data.ts$mRNA[data.ts$geno == "wt"]/fit.res$simScale.wt
  data.ts$mRNA[data.ts$geno == "mt"] <- data.ts$mRNA[data.ts$geno == "mt"]/fit.res$simScale.mt
  rownames(data.ts) <- NULL
  data.ts
}


tmp<-sapply(1:nrow(r2.final), function(i) {print(i);getBestFit(rownames(r2.final)[i])},
            simplify = F);names(tmp) <- NULL ;tmp <- do.call('rbind',tmp)
tmp$time <- tmp$time/60

write.csv(file='bestFit.tc.csv',tmp,row.names = F)



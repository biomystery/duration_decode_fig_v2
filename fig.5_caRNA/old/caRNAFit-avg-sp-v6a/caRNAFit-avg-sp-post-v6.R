rm(list=ls())
require(RColorBrewer)
require(pheatmap)

load(file = './caRNAFit-avg-sp-pre-v6.Rdata')
load(file = './caRNAFit-avg-sp-run-v6.Rdata')
source(file = 'caRNAFit-avg-sp-funs-v6.R')

nm <- 1:length(caFit)
for(i in 1:max(nm)) nm[i] <- caFit[[i]]$gene
names(caFit) <- nm

caFit.R2 <- sapply(1:length(caFit), calR2)
names(caFit.R2) <- nm

# parameters 
getBestFitPar <-function(g=names(caFit)[1]){ # rmsd.final
  df <- caFit[[g]]$par
  names(df) <- c("kt","k_deg","tau") 
  df
}

res <- data.frame(r2=caFit.R2,
                  gene=names(caFit.R2))
res <- cbind(res,t(sapply(names(caFit.R2), getBestFitPar)))
write.csv(file='caRNAFit-avg-sp-funs-v6-r2.csv',res,row.names = T)

# load simulation data  ---------------------------------------------------

n.gene <- length(caFit.R2)
getBestFit <-function(g=names(caFit.R2)[1]){
  g.res <- caFit.R2[g]
  plot.tp_step <- 5 # 5mins step
  fit.res <- caFit[[g]]
  all.tps <- fit.res$bestFit[[1]]$time
  plot.tp_idx <- all.tps %in% seq(all.tps[1],all.tps[length(all.tps)],by=plot.tp_step)
  n.tps <- sum(plot.tp_idx)
  data.ts <- rbind(fit.res$bestFit[[1]][plot.tp_idx,],
                   fit.res$bestFit[[2]][plot.tp_idx,],
                   fit.res$bestFit[[3]][plot.tp_idx,],
                   fit.res$bestFit[[4]][plot.tp_idx,])
  data.ts$sti <- rep(c("lps","tnf","lps","tnf"),each=n.tps)
  data.ts$geno <- rep(c("wt","wt",'mt',"mt"),each=n.tps)
  data.ts$gene <- rep(g,n.tps*4)
  data.ts$mRNA[data.ts$geno=="wt"]<- data.ts$mRNA[data.ts$geno=="wt"]/fit.res$simScale.wt
  data.ts$mRNA[data.ts$geno=="mt"]<- data.ts$mRNA[data.ts$geno=="mt"]/fit.res$simScale.mt
  rownames(data.ts) <- NULL
  data.ts
}


tmp<-sapply(1:n.gene, function(i) {print(i);getBestFit(names(caFit)[i])},
            simplify = F);names(tmp) <- NULL ;tmp <- do.call('rbind',tmp)
tmp$time <- tmp$time/60

write.csv(file='bestFit.tc.csv',tmp,row.names = F)


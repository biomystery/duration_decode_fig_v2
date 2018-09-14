rm(list=ls())
require(RColorBrewer)
require(pheatmap)
source(file = 'caRNAFit-avg-sp-funs-v6.R')
load(file = './caRNAFit-avg-sp-pre-v6.Rdata')
load(file = './caRNAFit-avg-sp-run-v6.Rdata')


nm <- 1:length(caFit)
for(i in 1:max(nm)) nm[i] <- caFit[[i]]$gene
names(caFit) <- nm

caFit.R2 <- sapply(1:length(caFit), calR2)
caFit.nRmsd <- sapply(1:length(caFit), function(x) caFit[[x]]$objective)
names(caFit.nRmsd)<- names(caFit.R2) <- nm
#caFit.R2 <-caFit.R2[order(caFit.R2,decreasing = T)]

genes.clust <- read.csv(file='../mRNA-Fit-avg-sp-v1e/mRNA.cluster.csv',stringsAsFactors = F)
pd <- data.frame(r2=caFit.R2,
                 nrmsd = caFit.nRmsd,
                 gene=names(caFit.R2),
                 stringsAsFactors = F)
all.equal(pd$gene,genes.clust$X)
pd <- data.frame(r2=caFit.R2,
                 nrmsd = caFit.nRmsd,
                 gene=names(caFit.R2),
                 cluster = genes.clust$cluster,
                 stringsAsFactors = F)

write.csv(file='caRNAFit-avg-sp-funs-v6-r2.csv',data.frame(r2=caFit.R2,
                                                           nrmsd = caFit.nRmsd,
                                                           gene=names(caFit.R2),
                                                           cluster=genes.clust$cluster),row.names = T)
# plot the results 
require(ggplot2)
require(dplyr)
pd %>% group_by(cluster) %>% summarise(ngenes = sum(nrmsd<=0.13))
ggplot(pd,aes(factor(cluster),nrmsd)) + geom_violin()+
  geom_hline(yintercept = .15) + geom_point(position = position_jitter(width = 0.2))
ggplot(pd,aes(factor(cluster),r2)) + geom_violin()+
  geom_hline(yintercept = .7) + geom_point(position = position_jitter(width = 0.2))


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


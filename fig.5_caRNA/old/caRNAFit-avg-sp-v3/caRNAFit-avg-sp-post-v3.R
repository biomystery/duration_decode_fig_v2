rm(list=ls())
require(RColorBrewer)
require(pheatmap)

load(file = './caRNAFit-avg-sp-pre-v2.Rdata')
load(file = './caRNAFit-avg-sp-run-v2.Rdata')
source(file = 'caRNAFit-avg-sp-funs-v2.R')

nm <- 1:length(caFit)
for(i in 1:max(nm)) nm[i] <- caFit[[i]]$gene
names(caFit) <- nm

caFit.R2 <- sapply(1:length(caFit), calR2)
names(caFit.R2) <- nm

ord <- order(caFit.R2,decreasing = T) 

pdf(file='caRNAFit-avg-sp-v2.pdf')
hist(caFit.R2,breaks = c(min(caFit.R2),seq(0,1,length.out = 11)),freq = T,xlim = c(0,1))
lapply(caFit[ord], plotFitInd)
dev.off()





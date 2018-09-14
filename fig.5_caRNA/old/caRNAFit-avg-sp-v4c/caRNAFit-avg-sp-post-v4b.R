rm(list=ls())
require(RColorBrewer)
require(pheatmap)

load(file = './caRNAFit-avg-sp-pre-v4b.Rdata')
load(file = './caRNAFit-avg-sp-run-v4b.Rdata')
source(file = 'caRNAFit-avg-sp-funs-v4b.R')


caFit.R2 <- sapply(1:length(caFit), calR2)
caFit.obj <- unlist(lapply(caFit, function(x) x$objective))
caFit.pars <- matrix(unlist(lapply(caFit, function(x) x$par)),ncol = 3,byrow = T)


names(caFit.R2) <- names(caFit)

ord <- order(caFit.R2,decreasing = T) 

write.csv(file='fit.result.csv',data.frame(r2=caFit.R2[ord],rmsd=caFit.obj[ord]),row.names = T)

pdf(file='caRNAFit-avg-sp-v4b.pdf')
hist(caFit.R2,breaks = c(min(caFit.R2),seq(0,1,length.out = 11)),freq = T,xlim = c(0,1))
lapply(caFit[ord], plotFitInd)
dev.off()



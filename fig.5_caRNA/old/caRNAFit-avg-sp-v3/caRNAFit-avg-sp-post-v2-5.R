rm(list=ls())
require(RColorBrewer)
require(pheatmap)

source(file = './caRNAFit-avg-sp-v2/caRNAFit-avg-sp-funs-v2.R')
load(file = './caRNAFit-avg-sp-v2/caRNAFit-avg-sp-pre-v2.Rdata')
files <- paste0("./caRNAFit-avg-sp-v",2:5,"/caRNAFit-avg-sp-run-v",2:5,'.Rdata')

for (i in 1:length(files)){
  ev <- new.env()
  load(file = files[i],ev)
  eval(expr = parse(text = paste0("caFit",i,".R2 <- sapply(1:length(ev$caFit), calR2)")))
}

# v2 - fit kt; v3 - delay + kt , v4 - delay + kt + kdeg , v5 - kt + kdeg
models <- c('kt','tau+kt','tau+kt+kdeg','kt+kdeg')
cols <- c('#999999','#E69F00',"dodgerblue","blueviolet","brown1","brown4")
plot(density(caFit1.R2),col=cols[1],main='R2 distribution',xlab = "R2",lwd=2)
lines



ord <- order(caFit.R2,decreasing = T) 

pdf(file='caRNAFit-avg-sp-v2.pdf')
hist(caFit.R2,breaks = c(min(caFit.R2),seq(0,1,length.out = 11)),freq = T,xlim = c(0,1))
lapply(caFit[ord], plotFitInd)
dev.off()





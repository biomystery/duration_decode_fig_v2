#setwd("~/Dropbox/Projects/DurationDecoding-code/Fig_code/Decode_figs/fig.1/")
rm(list=ls())

require(RColorBrewer)
require(pheatmap)

deltaAuc.all <- read.csv(file='../fig.1_wt/delta.avg.std.csv',header = T,row.names = 1,stringsAsFactors = F)
head(deltaAuc.all)

load("../fig.1_wt/kcluster_b1.Rdata")
pdf(file="spcificity_change2.pdf")


plot(deltaAuc.all$b1.TL.ctrl,deltaAuc.all$b1.TL.ko,pch=21,
     bg=ifelse(deltaAuc.all$b1.TL.ctrl<=-0.5 & deltaAuc.all$b1.TL.ko>-0.5 &deltaAuc.all$b1.TL.ko<0.5,
               'red','grey'),cex=1.5,
     xlim = c(-1.5,1.5),ylim = c(-1.5,1.5),
     xlab="TL specificity (ctrl.)",
     ylab="TL Specificity (ko.)")
abline(v=c(-.5,.5),lty=2)
abline(h=c(-.5,.5),lty=2)
abline(a=0,b=1,lwd=2,col="gray80")

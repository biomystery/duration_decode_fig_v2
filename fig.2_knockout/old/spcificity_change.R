#setwd("~/Dropbox/Projects/DurationDecoding-code/Fig_code/Decode_figs/fig.1/")
rm(list=ls())

require(RColorBrewer)
require(pheatmap)

deltaAuc.all <- read.csv(file='../fig.1_wt/delta.avg.std.csv',header = T,row.names = 1,stringsAsFactors = F)
head(deltaAuc.all)

load("../fig.1_wt/kcluster_b1.Rdata")
pdf(file="spcificity_change.pdf")

cols <- rainbow(8)
all.equal(k_ord$gene,rownames(deltaAuc.all))
plot(deltaAuc.all$b1.TL.ctrl,deltaAuc.all$b1.TL.ko,pch=21,bg=cols[k_ord$cluster],
     xlab="TL specificity (ctrl.)",
     ylab="TL Specificity (ko.)")
abline(v=c(-.5,.5),lty=2)
abline(h=c(-.5,.5),lty=2)
abline(a=0,b=1,lwd=2,col="gray80")
legend(0.6,0.6,legend = paste0('c',1:8),pch=16,col = cols,border = 'n')
#identify(deltaAuc.all$b1.TL.ctrl,deltaAuc.all$b1.TL.ko,labels = rownames(deltaAuc.all))

par(mfrow=c(3,3))
idx <- c(3  ,27 ,64 ,73,130 ,154)
genes.zscale <- read.csv(file='../fig.1_wt/zscale.csv',stringsAsFactors = F,
                         header = T,row.names = 1)
genes.select <- c(rownames(deltaAuc.all)[idx],"Fos","Ccl5")
all.equal(rownames(genes.zscale),rownames(deltaAuc.all))

fun.plotZscaleGene <- function(gname=genes.select[1]){
  pd <- genes.zscale[gname,c(1:14,22:35)]
  pd <- matrix(pd,nrow = 7)
  matplot(data.frame(time=c(0,.5,1,2,3,5,8)),pd,type='l',
          pch = 21,col=c(1,2,1,2),lwd=2,lty=c(1,1,2,2),
          xlab='Time(hr)',ylab = paste0(gname,'(z)'),
          main=paste0('sp.ctrl=',signif(deltaAuc.all[gname,'b1.TL.ctrl'],2),
                      ",sp.mt=",signif(deltaAuc.all[gname,'b1.TL.ko'],2)))
  matpoints(data.frame(time=c(0,.5,1,2,3,5,8)),pd,
          pch = 21,col=c(1,2,1,2),bg=c("grey",'lightcoral',"white","white"),lwd=2)
  }

sapply(genes.select, fun.plotZscaleGene)

dev.off()
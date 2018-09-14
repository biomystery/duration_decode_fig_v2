#setwd("~/Dropbox/Projects/DurationDecoding-code/Fig_code/Decode_figs/fig.1/")
rm(list=ls())

require(RColorBrewer)
require(pheatmap)

deltaAuc.all <- read.csv(file='../fig.1_wt/data/sp.fraction.maxscale.csv',header = T,row.names = 1,stringsAsFactors = F)
#deltaAuc.all <- log2(deltaAuc.all)
head(deltaAuc.all)

load("../fig.1_wt/data/kcluster_b1.Rdata")
#pdf(file="spcificity_change_max.pdf")

cols <- rainbow(8)
all.equal(k_ord$gene,rownames(deltaAuc.all))
plot(deltaAuc.all$b1.TL.ctrl,deltaAuc.all$b1.TL.ko,pch=21,bg=cols[k_ord$cluster],
     xlab="TL specificity (ctrl.)",
     ylab="TL Specificity (ko.)")
abline(v=c(-.5,.5),lty=2)
abline(h=c(-.5,.5),lty=2)
abline(a=0,b=1,lwd=2,col="gray80")
legend(-6,-1,legend = paste0('c',1:8),pch=16,col = cols,border = 'n')

#dev.off()
all.equal(k_ord$gene,rownames(deltaAuc.all))
pd <- data.frame(ctrl=deltaAuc.all$b1.TL.ctrl,ako=deltaAuc.all$b1.TL.ko)
write.csv(file='frac.sp.maxscale.csv',cbind(deltaAuc.all[,c("b1.TL.ctrl","b1.TL.ko")],cluster=k_ord$cluster))


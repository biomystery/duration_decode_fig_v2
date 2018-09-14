rm(list=ls())
require(RColorBrewer)
require(pheatmap)
sp.caRNA <- read.csv(file='./data/caRNA-zscale-sp.csv',header = T,
                     row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../fig.1_wt/delta.avg.std.csv",header = T,
                    row.names = 1,stringsAsFactors = F)


dim(sp.caRNA)
dim(sp.mRNA)

all.equal(rownames(sp.caRNA),rownames(sp.mRNA))

pdf(file='./fig/cmp-caRNA-mRNA-zsp.pdf')
par(mfrow=c(2,2))
pd <- data.frame(tlsp.mRNA.ctrl = sp.mRNA$b1.TL.ctrl,
                 tlsp.caRNA.ctrl=sp.caRNA$ctrl.tlsp,
                tlsp.mRNA.ko=sp.mRNA$b1.TL.ko,
                tlsp.caRNA.ko= sp.caRNA$ko.tlsp)
                
xylims <- c(min(pd),max(pd))                 
plot(pd[,1:2],pch=21,bg="grey",xlim=xylims,ylim=xylims,
     main=paste("cor =",signif(cor(pd[,1],y=pd[,2]),2)))
abline(a=0,b=1)

abline(v=c(-.5,.5))
abline(h=c(-.5,.5))

plot(pd[,3:4],pch=21,bg="grey",xlim=xylims,ylim=xylims,
     main=paste("cor =",signif(cor(pd[,3],y=pd[,4]),2)))
abline(a=0,b=1)
cor(pd[,3],y=pd[,4])
abline(v=c(-.5,.5))
abline(h=c(-.5,.5))




plot(pd[,c(2,4)],pch=21,bg=ifelse(pd[,1]<=-.5 & pd[,3]> -.5 & pd[,3]<.5,'coral1',"grey"),
     xlim=xylims,ylim=xylims,
     main=paste0('cor=',signif(cor(pd[,2],pd[,4]),2)))
abline(a=0,b=1)
abline(v=c(-.5,.5))
abline(h=c(-.5,.5))

plot(pd[,c(1,3)],pch=21,bg=ifelse(pd[,1]<=-.5 & pd[,3]> -.5 & pd[,3]<.5,'coral1',"grey"),
     xlim=xylims,ylim=xylims,
     main=paste0('cor=',signif(cor(pd[,1],pd[,3]),2)))
abline(a=0,b=1)
abline(v=c(-.5,.5))
abline(h=c(-.5,.5))

par(mfrow=c(1,1))
range(pd)
bks <- c(-1.5,-1,-0.5,0,.5,1,1.5)
cols <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                    "RdYlBu")))(length(bks)-1)
cols[4]<- cols[3]
pheatmap(pd, breaks = bks,color = cols,cutree_rows = 8, #kmeans_k = 8,
         scale = 'none',cluster_cols = F,show_rownames = F)

dev.off()
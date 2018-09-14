v1.fit <- read.csv(file = './data/v1-ctrl-goodness-of-fit.csv',header = T,row.names = 1)
tlsp.ctrl <- read.csv(file = '../fig.1_wt/delta.avg.std.csv',header = T,row.names = 1)
all.equal(rownames(v1.fit),rownames(tlsp.ctrl))
tlsp.ctrl <- tlsp.ctrl[rownames(v1.fit),]
colnames(tlsp.ctrl)
tlsp <- tlsp.ctrl[,c(1,4)]
pd <- cbind(v1.fit,tlsp)
colnames(pd)

range(pd$r2)
bks <- c(min(pd$r2),-1,-0.5,0,.5,.75,1)
bks <- c(min(pd$r2),-1,0,.75,1)
#cols<- terrain.colors(length(bks)-1)[cut(pd$r2,breaks = bks)]
cols<- colorRampPalette(colors = c("black","yellow1"))(length(bks)-1)

nf <- layout(matrix(c(1,2),byrow = T,nrow = 1),
       widths = c(5,2))
layout.show(nf)

#pdf(file = './fig/cmp-v1fit-TLsp-dy.pdf')
plot(pd$b1.TL.ctrl,pd$b1.TL.ko,pch=21,
     cex=1,bg=cols[cut(pd$r2,breaks = bks)],
     xlim = c(-1.5,1.5),ylim = c(-1.5,1.5),
     xlab="TL specificity (ctrl.)",
     ylab="TL Specificity (ko.)")
abline(v=c(-.5,.5),lty=2)
abline(h=c(-.5,.5),lty=2)
abline(a=0,b=1,lwd=2,col="gray80")

plotLegend <- function(cols,bks){
  require(RColorBrewer)
  #bks <- seq(round(min(rsums)),round(max(rsums))+1,length.out = 6)
  #cols<- colorRampPalette(c( "white", "blueviolet"))(5)
  #fnames<-'tmp.eps'
  #par(mar=c(0.2, 0.4, 0.4,0.1))
  barplot(rep(1,length(bks)-1),width=diff(bks),space = 0,border = NA,
          col = cols,axes = F,horiz = T,
          xaxs = "i", yaxs = "i",xlim = c(0,1),
          ylim=c(0,sum(diff(bks))),ylab="Fit (R2)"
          #xlim=c(0,1)
  )
  axis(2,labels = signif(bks[-1],2),tck = 0.2,at=cumsum(diff(bks)),las=2)
  box()
  
}
plotLegend(cols,bks)

#dev.off()


#

#dev.off()
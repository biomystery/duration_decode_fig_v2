v1.fit <- read.csv(file = './data/v1-ctrl-goodness-of-fit.csv',header = T,row.names = 1)
tlsp.ctrl <- read.csv(file = '../fig.1_wt/delta.avg.std.csv',header = T,row.names = 1)
all.equal(rownames(v1.fit),rownames(tlsp.ctrl))
tlsp.ctrl <- tlsp.ctrl[rownames(v1.fit),]
tlsp.ctrl <- tlsp.ctrl[,1]
pd <- cbind(v1.fit,tlsp.ctrl)


pdf(file = './fig/cmp-v1fit-TLsp.pdf')
par(mfrow=c(2,2))
fun.plotCmp <- function(x=pd$r2,y=pd$tlsp.ctrl,labx="r2",laby='tlsp.ctrl',
                        th =.75){
  lm.fit <- lm(y~x)
  plot(x,y,pch=21,bg=ifelse(x>=th,2,"grey"),
       xlab = labx,ylab = laby,
       main=paste('lm.r2=',signif(summary(lm.fit)$r.squared,2),
                  'cor=',signif(cor(x,y),2),
                  'n=',sum(x>th)))
  abline(v = th,col=2)
  abline(h=c(-.5,.5),lty=2)
  abline(lm.fit,lwd=2)
}

fun.plotCmp.2 <- function(x=pd$rmsd,y=pd$tlsp.ctrl,labx="rmsd",laby='tlsp.ctrl',
                        th =.15){
  lm.fit <- lm(y~x)
  plot(x,y,pch=21,bg=ifelse(x<=th,2,"grey"),
       xlab = labx,ylab = laby,
       main=paste('lm.r2=',signif(summary(lm.fit)$r.squared,2),
                  'cor=',signif(cor(x,y),2),
                  'n=',sum(x<=th)))
  abline(v = th,col=2)
  abline(h=c(-.5,.5),lty=2)
  abline(lm.fit,lwd=2)
}

fun.plotCmp()
fun.plotCmp.2()
dev.off()
# Check hf vs. sp ---------------------------------------------------------

genes.hf <- read.csv(file="v4-hf-final.csv",row.names = 1,stringsAsFactors = F)
#genes.sp <- read.csv(file="../fig.1_wt/delta.avg.std.csv",row.names = 1,stringsAsFactors = F)
genes.sp <- read.csv(file="../fig.1_wt/delta.avg.std.csv",row.names = 1,stringsAsFactors = F)
genes.sp <- genes.sp[rownames(genes.hf),]
colnames(genes.sp)

construct.raw <- read.csv(file="./hf_construct.csv",stringsAsFactors = F,header = T)
construct.raw[,2:ncol(construct.raw)] <-  2^construct.raw[,2:ncol(construct.raw)]

pd <- data.frame(hf_log10=genes.hf[,1],
                 lt_sp = genes.sp[,1])

pd.cor <- cor.test(x = pd$hf_log10,y = pd$lt_sp)

pdf(file='sp_ctl_lhf_cor.pdf')
par(mfrow=c(2,2))
plot(pd$hf_log10,(pd$lt_sp),cex=.75,pch=21,bg="grey",
     xlab = 'half-life(log10,mins)',
     ylab="TL specificity",
     main=paste0('cor=',signif(pd.cor$estimate,2),
                 ',p=',signif(pd.cor$p.value,2)))

#identify(x = pd$hf_log10,y = pd$lt_sp,labels = rownames(pd))
abline(lm(pd$lt_sp~pd$hf_log10),lwd=2,col=2)
boxplot(pd$lt_sp~ cut(pd$hf_log10,
                      c(1.0,1.5,2.0,2.5,3.0,3.5)),
                      #seq(min(pd$hf_log10),max(pd$hf_log10),length.out = 5)),
        outline = F,
        xlab = 'half-life(log10,mins)',
        ylab="TL specificity")


pd.tmp <- c(construct.raw[,2],construct.raw[,3])
pd.sp.ccl5 <- sum(construct.raw$ccl5.tnfp - construct.raw$ccl5.lpsp)/sd(pd.tmp)/4
pd.tmp <- c(construct.raw[,5],construct.raw[,6])
pd.sp.fos <- sum(construct.raw$fos.tnfp - construct.raw$fos.lpsp)/sd(pd.tmp)/4

matplot(construct.raw$Time,log2(construct.raw[,2:3]),type = "b",
        pch=21,lty=1,ylab = "log2 fold",xlab="Time (hr)",
        bg=c("grey",'lightcoral'),
        main=paste0("Ccl5,TL_sp=",signif(pd.sp.ccl5,2)))
legend(4,.4,legend = c("TNFp","LPSp"),lty=1,col=c(2,1))

matplot(construct.raw$Time,log2(construct.raw[,5:6]),type = "b",
        pch=21,lty=1,ylab = "log2 fold",xlab="Time (hr)",
        bg=c("grey",'lightcoral'),
        main=paste0("Fos,TL_sp=",signif(pd.sp.fos,2)))
par(mfrow=c(1,1))

dev.off()
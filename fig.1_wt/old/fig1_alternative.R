#setwd("~/Dropbox/Projects/DurationDecoding-code/Fig_code/Decode_figs/fig.1/")
rm(list=ls())

require(RColorBrewer)
require(pheatmap)
source('./auxilary_functions.R')
deltaAuc.all <- read.csv(file='./delta.avg.std.csv',header = T,row.names = 1,stringsAsFactors = F)
head(deltaAuc.all)

range(deltaAuc.all)
fun.getKord <- function(data=deltaAuc,k=10){
  set.seed(200)
  kclust <- kmeans(data,k,nstart = 1000,iter.max = 20000)
  data.frame(ord=unlist(sapply(1:k, function(x) which(kclust$cluster==x))),
                        clust = unlist(sapply(1:k, function(x) rep(x,sum(kclust$cluster==x)))))
}

auc.kord <- fun.getKord(data = deltaAuc.all[,1:6])


rainbow(n=10)
anno.row <- data.frame(clust= auc.kord$clust)
rownames(anno.row) <- rownames(deltaAuc.all[auc.kord$ord,])
anno.colors <- list(clust=rainbow(n=10))
pdf(file="fig1_alternative.pdf")
par(mfrow=c(2,2))


cols <- brewer.pal(12,"Paired")
boxplot(deltaAuc.all[,1:3],ylab='specificity',col=cols[c(2,4,6)])
abline(h = c(-.5,0,.5),lty=2,col='grey')
boxplot(deltaAuc.all[,c(1,4,2,5,3,6)],ylab='specificity',col=cols[c(2,1,4,3,6,5)],
        xaxt='n',outline = F)
axis(1,at=1:6,labels =colnames(deltaAuc.all[,c(1,4,2,5,3,6)]),las=2 )
abline(h = c(-.5,0,.5),lty=2,col='grey')


pd.tmp <- hist(deltaAuc.all[,1],breaks = c(min(deltaAuc.all[,1]),-.5,.5,max(deltaAuc.all[,1])),plot = F)$counts
names(pd.tmp) <-  c("LPS.sp",'none.sp',"TNF.sp")
barplot(pd.tmp,xlab = "TL.sp.ctrl",col  = c("dodgerblue","azure","yellow2"),
        ylim=c(0, max(pd.tmp)+15),ylab='# of genes')
box(); grid()

par(mfrow=c(1,1))
bks <- c(-1.5,-1,-.5,0,.5,1,1.5)
cols <- rev(colorRampPalette(c("yellow2","azure","dodgerblue"))(length(bks-1)))
cols <- rev(colorRampPalette(c("yellow2","white","dodgerblue"))(6))

cols[3] <- cols[4] <-"azure"

pheatmap(deltaAuc.all[auc.kord$ord,1:6],breaks = bks,color = cols,cellwidth = 24,cellheight = 2,
         scale = "none",show_rownames = F,cluster_cols = F,cluster_rows = F,
         gaps_row =  sapply(1:10,function(x) max(which(auc.kord$clust==x))),
         border_color = "black", gaps_col = 3,
         annotation_row = anno.row,annotation_colors = anno.colors,
         main="Specificity (std)",annotation_legend = F
)

dev.off()

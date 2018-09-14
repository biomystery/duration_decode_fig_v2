
# round 2: at 2016-10-24 ----------------------------------------------
# scaled to geno

rm(list = ls())
filename <- './data/maxscale.csv'
pd.rnaseq <- read.csv(file = filename,stringsAsFactors = F,
               row.names = 1)

#getSp moved to auxilary_functions 


### average metrics 
rsums <- sapply(0:11, function(i)
  rowSums(pd.rnaseq[,(i*7+1):((i+1)*7)])); 
sp.avg <- getSp(metrics = rsums)
write.csv(file='sp.avg.csv',sp.avg,row.names = T)

### peak metrics 
rsums <- sapply(0:11, function(i)
  apply(pd.rnaseq[,(i*7+1):((i+1)*7)], 1,max))
sp.peak <- getSp(metrics = rsums)


### write files 
pd <- cbind(sp.avg[,c("b1.LT.ctrl","b1.LT.ko")],
            sp.peak[,c("b1.LT.ctrl","b1.LT.ko")])
colnames(pd) <- c("sp.LT.avg.ctrl","sp.LT.avg.ko","sp.LT.peak.ctrl","sp.LT.peak.ctrl")
pd <- cbind(pd,sp.avg.ratio=pd[,1]/pd[,2],sp.peak.ratio=pd[,3]/pd[,4])
pd <- pd[,c(1,2,5,3,4,6)] ;rm(list = c("sp.avg",'sp.peak'))



### gene name convert & add cluster info 

load("./data/kcluster_b1.Rdata")
all.equal(k_ord$gene,rownames(pd)) ;all.equal(k_ord$gene,rownames(pd.rnaseq))

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")

tmp <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
      values = grep('ENSMUSG',k_ord$gene,value = T), mart = ensembl)
rownames(tmp) <- tmp$ensembl_gene_id

k_ord$gene[grep('ENSMUSG',k_ord$gene)] <- tmp[k_ord$gene[grep('ENSMUSG',k_ord$gene)],"mgi_symbol"]

rownames(pd) <- k_ord$gene; rownames(pd.rnaseq) <- k_ord$gene
#pd <- cbind(pd,cluster=k_ord$cluster)
write.csv(file='b1.sp.csv',pd,row.names = T)


### applying filtering 
th.frac <- 2^.5; th.ratio <- 2^0.5
lt.sp.avg.genes <- pd$sp.avg.ratio >= th.ratio # 1 in ctrl 
lt.sp.peak.genes <- pd$sp.peak.ratio >= th.ratio 
lt.sp.genes <- lt.sp.avg.genes | lt.sp.peak.genes
require(gplots)
venn(list(avg.sp.genes=which(lt.sp.avg.genes),
          peak.sp.genes=which(lt.sp.peak.genes)))

venn(data.frame(ctrl.lps=pd$sp.LT.avg.ctrl>=th.frac,
                 ctrl.tnf=pd$sp.LT.avg.ctrl<= 1/th.frac,
                 dy = lt.sp.avg.genes))


### clustering on sp 
require(pheatmap)
pd.sp <- pd ; rm(pd);
pd.sp.log2.cut <-pd.sp.log2 <- log2(pd.sp)
pd.sp.log2.cut[pd.sp.log2>2.5] <- 2.5; pd.sp.log2.cut[pd.sp.log2 < -2.5] <- -2.5

ncl <- function(mydata, k=min(20, nrow(mydata) - 1)) {
  wss <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))
  for(i in unique(c(2:10, round(seq(2, k, length.out=50))))) wss[i] <- kmeans(mydata, centers=i, iter.max=200, nstart=20)$tot.withinss
  plot(1:k, wss, xlab="Number of Clusters", ylab="Within groups sum of squares")
  return(wss)
}
wss <- ncl(pd.sp.log2, k=nrow(pd.rnaseq) - 1); abline(v=12,lwd=2,col=2)

k=12
set.seed(100)
kord <- kmeans(pd.sp.log2[,1:3],k,iter.max = 1000,nstart = 1000)
kord.ord.sp <- unlist(sapply(1:k, function(x) which(kord$cluster==x)))
kord.gap.sp <- unlist(sapply(1:k, function(x) which.max(kord$cluster==x)))

## hm. 1
require(RColorBrewer)
require(pheatmap)
cols <- colorRampPalette(c('azure',"gold"))(5)
cols.ratio <- c(rep("azure",5),cols)
bks.ratio <- c(-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5)
hm.1.ratio <- pheatmap(pd[kord.ord.sp,"sp.avg.ratio"],scale = 'none',cluster_cols = F,cluster_rows = F,
                 breaks = bks.ratio,color = cols.ratio,show_rownames = F,
                 show_colnames = F,cellwidth = 24) #score hm

cols <- colorRampPalette(c('azure',"green4"))(5)
cols.sp <- c(rep("azure",5),cols)
bks.sp <- c(-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5)
hm.1.sp <- pheatmap(pd.sp.log2.cut[kord.ord.sp,1:2],scale = 'none',cluster_cols = F,cluster_rows = F,
                 breaks = bks.sp,color = cols.sp,show_rownames = F,
                 show_colnames = F,cellwidth = 24) #score hm

all.equal(rownames(pd.rnaseq),rownames(pd.sp))

##  (plotly the scatter)
require(plotly)
bks <- bks.ratio
bks[length(bks)] <- max(pd.sp.log2$sp.avg.ratio)
pd3 <- cbind(pd.sp.log2,gene=rownames(pd.sp.log2),ratio=cut(pd.sp.log2$sp.avg.ratio,breaks = bks))
#colScale <- scale_color_manual(name='ratio',values = cols)
colScale <- scale_fill_manual(name='ratio',values = cols.ratio,breaks = levels(pd3$ratio))
p <- ggplot(aes(x=sp.LT.avg.ctrl,y=sp.LT.avg.ko),data=pd3)+
  geom_point(size=2,aes(text = paste("Gene:",gene),fill=ratio),shape=21,stroke=0.1,
             colour='black') + colScale + theme_bw()
(gg<- ggplotly(p))


## scatter 
bks.ratio.2 <- bks.ratio
bks.ratio.2[length(bks.ratio.2)] <- max(pd.sp.log2$sp.avg.ratio)
plot(pd.sp.log2$sp.LT.avg.ctrl,pd.sp.log2$sp.LT.avg.ko,
     xlim = range(pd.sp.log2[,1:2]),ylim=range(pd.sp.log2[,1:2]),
     pch=21,bg=cols.ratio[cut(pd.sp.log2$sp.avg.ratio,bks.ratio.2)],
     xlab='LP.sp.ctrl',ylab = 'LP.sp.ko')
abline(h=c(.5,-.5),lty=2); abline(v=c(.5,-.5),lty=2)
sapply(c(-0.5,0,0.5),function(x) abline(a=x,b=1,lty=2))
grid()



### RNA Seq
nbks.rnaseq <- 11
cols.rnaseq <- colorRampPalette(c("azure","firebrick1"))(nbks.rnaseq-1)
bks.rnaseq <- seq(0,1,length.out = nbks.rnaseq)
#anno.row <- data.frame(peak.gene = as.factor(lt.sp.peak.genes),
#                      avg.gene = as.factor(lt.sp.avg.genes))
#rownames(anno.row) <- rownames(pd.rnaseq)

hm.1.rnaseq <- pheatmap(pd.rnaseq[kord.ord.sp ,c(1:14,22:35)],
                 breaks = bks.rnaseq,
                 color = cols.rnaseq,
#                 annotation_row = anno.row,annotation_names_row=F,
                 show_rownames = F,
                 scale = 'none',cellwidth = 12,cellheight = 2,
                 cluster_cols = F,cluster_rows = F,show_colnames = F ,
                 gaps_col = c(7,14,21))

### filtering 
#idx <- anno.row$peak.gene == 'TRUE' | anno.row$avg.gene =="TRUE"
idx <- pd.sp.log2[kord.ord.sp,3] >= 0.5
names(idx) <- names(kord.ord.sp)

#anno.row <- anno.row[anno.row$peak.gene | anno.row$avg.gene,]


hm.3.ratio <- pheatmap(pd.sp.log2.cut[names(which(idx)),3],scale = 'none',cluster_cols = F,cluster_rows = F,
                 breaks = bks.ratio,color = cols.ratio,show_rownames = F,
                 #gaps_row = kord.gap.sp,
                 show_colnames = F,cellwidth = 24)

hm.3.sp <- pheatmap(pd.sp.log2.cut[names(which(idx)),1:2],scale = 'none',cluster_cols = F,cluster_rows = F,
                       breaks = bks.sp,color = cols.sp,show_rownames = F,
                       show_colnames = F,cellwidth = 24)


hm.3.rnaseq <- pheatmap(pd.rnaseq[names(which(idx)),c(1:14,22:35)],
                 breaks = bks.rnaseq,
                 color = cols.rnaseq,
#                 annotation_row = anno.row,annotation_names_row=F,
                 show_rownames = T,cutree_rows = 10,
                 scale = 'none',cellwidth = 12,
                 cluster_cols = F,cluster_rows = F,show_colnames = F ,
                 gaps_col = c(7,14,21))



write.table(file='sp.rnaseqFrac.csv',sep = ',',
          cbind(pd[kord.ord[names(which(idx))],],pd.rnaseq[kord.ord[names(which(idx))],c(1:14,22:35)]))

### final - pick up the 

write.csv(file = './data/rnaseq.avgDyRatio.csv',pd.rnaseq[kord.ord[names(which(idx))],c(1:14,22:35)],col.names = T)

pd.rnaseq <- read.csv(file = './')
pd.rnaseq <- pd.rnaseq[kord.ord[names(which(idx))],c(1:14,22:35)]
pd.sp.log2.avg <- pd.sp.log2[kord.ord[names(which(idx))],1:3]
pd.rnaseq <- pd.rnaseq[kord.ord[names(which(idx))],]
all.equal(rownames(pd.sp.log2.avg),rownames(pd.rnaseq))
ncl <- function(mydata, k=min(20, nrow(mydata) - 1)) {
  wss <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))
  bss <- wss; bss[1] <-0 
  for(i in unique(c(2:10, round(seq(2, k, length.out=50))))) {
    tmp <- kmeans(mydata, centers=i, iter.max=200, nstart=20)
    wss[i] <- tmp$tot.withinss
    bss[i] <- tmp$betweenss
    }
  plot(1:k, wss, xlab="Number of Clusters", ylab="Sum of squares")
  points(1:k,bss,col=2)
  return(list(wss,bss))
}

pd.rnaseq.cut <- pd.rnaseq[names(which(idx)),c(1:14,22:35)]

wss <- ncl(pd.rnaseq.cut, k=nrow(pd.rnaseq.cut) - 1); 
abline(v=12,lwd=2,col=2)

set.seed(100)
k=12
kord <- kmeans(pd.rnaseq.cut, centers=k, iter.max=200, nstart=20)
kord.ord.rnaseq <- unlist(sapply(1:k, function(x) which(kord$cluster==x)))
kord.gap.rnaseq <- unlist(sapply(1:k, function(x) which.max(kord$cluster==x)))


hm.4 <- pheatmap(pd.rnaseq.cut[kord.ord.rnaseq,],
                 breaks = bks.rnaseq,
                 color = cols.rnaseq,
                 gaps_row = kord.gap.rnaseq,
                 legend = F,
                 show_rownames = T,
                 scale = 'none',cellwidth = 12,
                 cluster_cols = F,cluster_rows = F,show_colnames = F ,
                 gaps_col = c(7,14,21))


pd.sp.log2.cut.2 <- pd.sp.log2.cut[names(which(idx)),]
all.equal(rownames(pd.rnaseq.cut),rownames(pd.sp.log2.cut.2))
hm.4.ratio <-  pheatmap(pd.sp.log2.cut.2[kord.ord.rnaseq,3],legend = F,
                        breaks = bks.ratio,
                        color = cols.ratio, 
                        gaps_row = kord.gap.rnaseq,
                        scale = 'none',cellwidth = 12,
                        cluster_cols = F,cluster_rows = F,
                        show_colnames = F,show_rownames=F)

hm.4.sp <-  pheatmap(pd.sp.log2.cut.2[kord.ord.rnaseq,1:2],legend = F,
                        breaks = bks.sp,
                        color = cols.sp, 
                        gaps_row = kord.gap.rnaseq,
                        scale = 'none',cellwidth = 12,
                        cluster_cols = F,cluster_rows = F,
                        show_colnames = F,show_rownames=F)


### plot examples 
colnames(pd3)
pd3.print<- signif(pd.sp.log2[,1:6],digits = 2)
plot.gene <- function(gene='Irg1'){
  #pd.rnaseq <- pd.rnaseq.cut
  pd.rnaseq <- pd.rnaseq[,c(1:14,22:35)]
  tm <- c(0,.5,1,2,3,5,8)
  plot(tm,pd.rnaseq[gene,1:7],pch=16,col=2,type='b',xlab='time (hr)',ylab='mRNA level',
       main=paste0("sp.ctrl=",pd3.print[gene,]$sp.LT.avg.ctrl,',sp.mt=',pd3.print[gene,]$sp.LT.avg.ko,
                   ',ratio=',pd3.print[gene,]$sp.avg.ratio), 
       ylim=c(0,1)) # LPS.ctrl
  text(.5,.9,gene)
  lines(tm,pd.rnaseq[gene,1:7 + 14],pch=1,col=2 ,type='b',lty=2) # lps.ko
  lines(tm,pd.rnaseq[gene,1:7 + 21],pch=1,col=1  ,type='b',lty=2) #tnf.ko
  lines(tm,pd.rnaseq[gene,1:7 + 7],pch=16,col=1  ,type='b',lty=1) #tnf.ctrl. 
  grid()
  sum(pd.rnaseq[gene,1:7])/sum(pd.rnaseq[gene,1:7+7])
  sum(pd.rnaseq[gene,1:7+14])/sum(pd.rnaseq[gene,1:7+21])
}

plot.gene('Rsad2')
plot.gene('Ifih1')

## group 2
## group 3
plot.gene.merge(c('Ccl2','Ccl7'))

## group 4
plot.gene('Ccl5')
plot.gene('Acpp')

## group 5
plot.gene('Cxcl1')
## group 6 
plot.gene("Relb")
plot.gene("Vcam1")
plot.gene("Nfkb2")
plot.gene("Nfkbie")
plot.gene("Tnfaip3")
plot.gene("Rel")

## group 7 
plot.gene("Ifih1")
plot.gene("Oasl1")

par(mfrow=c(2,1))
plot.gene.merge <- function(x) sapply(x,plot.gene)
plot.gene.merge(c("Ccl5","Acpp"))
plot.gene.merge(c("Irg1","Rsad2"))

sapply(c("Tnfaip3", "Stx11"), plot.gene)
sapply(c("Vcam1", "Relb"), plot.gene)
sapply(c("Ifih1", "Oasl1"), plot.gene)
sapply(c("Ccl1", "Cx3cl1"), plot.gene)
sapply(c("Gbp2", "Il34"), plot.gene) #1


#


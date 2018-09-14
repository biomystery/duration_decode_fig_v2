setwd('~/Dropbox/Projects/DurationDecoding-code/Fig_code/')

# mRNA Sp calculation ----------------------------------------------------
#at 2016-10-24
mRNA.rpkm <- read.csv(file = '../data/mRNA.nfkbgene.rpkm.csv',stringsAsFactors = F,row.names = 1,header = T)
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)


### average metrics 
rsums <- sapply(0:11, function(i)
  rowSums(mRNA.rpkm[,(i*7+1):((i+1)*7)])); 
mRNA.sp.avg <- log2(getSp(metrics = rsums))
mRNA.sp.avg.2 <- read.csv(file='./data/sp.avg.csv',stringsAsFactors = F,row.names = 1)
all.equal(mRNA.sp.avg,mRNA.sp.avg.2) #TRUE 
write.csv(file='./data/mRNA.sp.avg.csv',mRNA.sp.avg)

### peak metrics 
rsums <- sapply(0:11, function(i)
  apply(mRNA.rpkm[,(i*7+1):((i+1)*7)], 1,max))
mRNA.sp.avg.peak <- getSp(metrics = rsums)
mRNA.sp.avg.peak <- log2(mRNA.sp.avg.peak)
write.csv(file='./data/mRNA.sp.peak.csv',mRNA.sp.avg.peak)


mRNA.sp.avg.peak["Gda",]


# caRNA sp calculation ----------------------------------------------------

caRNA.max <- read.csv(file='./data/caRNA.cpm.max.geno.scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
all.equal(rownames(caRNA.max),k_ord$gene.2)
caRNA.sum<- data.frame(ctrl.L = rowSums(caRNA.max[,1:6]),
                       ctrl.T = rowSums(caRNA.max[,7:12]),
                       mt.L = rowSums(caRNA.max[,13:18]),
                       mt.T = rowSums(caRNA.max[,19:24]))
sp.caRNA <- data.frame(LT.ctrl=caRNA.sum$ctrl.L/caRNA.sum$ctrl.T,
                         LT.mt = caRNA.sum$mt.L/caRNA.sum$mt.T)
sp.caRNA[apply(sp.caRNA,2,is.infinite)] <- NA
sp.caRNA <- log2(sp.caRNA)
rownames(sp.caRNA) <- rownames(caRNA.sum);all.equal(rownames(caRNA.sum),rownames(sp.mRNA))
write.csv(file = './data/caRNA.sp.avg.csv',sp.caRNA)


caRNA.sum<- data.frame(ctrl.L = apply(caRNA.max[,1:6],1,max),
                       ctrl.T = apply(caRNA.max[,7:12],1,max),
                       mt.L = apply(caRNA.max[,13:18],1,max),
                       mt.T = apply(caRNA.max[,19:24],1,max))
sp.caRNA <- data.frame(LT.ctrl=caRNA.sum$ctrl.L/caRNA.sum$ctrl.T,
                       LT.mt = caRNA.sum$mt.L/caRNA.sum$mt.T)
sp.caRNA[apply(sp.caRNA,2,is.infinite)] <- NA
sp.caRNA <- log2(sp.caRNA)
rownames(sp.caRNA) <- rownames(caRNA.sum);all.equal(rownames(caRNA.sum),rownames(sp.mRNA))
write.csv(file = './data/caRNA.sp.peak.csv',sp.caRNA)


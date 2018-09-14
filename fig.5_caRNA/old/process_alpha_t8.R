
fun.calAssignRatio <- function(file='../data//counts-gene.txt.summary.txt'){
  counts.summary<- read.delim(file=file,header = T,row.names = 1)
  counts.summary[1,]/colSums(counts.summary)
}

assignRatio.mine <- fun.calAssignRatio()
assignRatio.r2000 <- fun.calAssignRatio(file = "../data//counts.txt.summary-roberto-2000.txt")
assignRatio.r4000 <- fun.calAssignRatio(file = "../data//counts.txt.summary-roberto-4000.txt")
assignRatio.m4000 <- fun.calAssignRatio(file = "../data//counts.gene.summary-me-4000.txt")
assignRatio.m2000 <- fun.calAssignRatio(file = "../data//counts-gene.txt.summary-me-2000.txt")
lapply(list(assignRatio.mine,assignRatio.r2000,assignRatio.r4000), length)


par(mfrow=c(1,1))

plot(as.numeric(assignRatio.mine),main="me(k) vs. Roberto(r)",type='b',ylim=c(0,1),
     xlim=c(0,22))
lines(as.numeric(assignRatio.r2000),type='b',col=2)
lines(as.numeric(assignRatio.r4000),type='b',col=2)
lines(as.numeric(assignRatio.m4000),type='b',col=1)
lines(as.numeric(assignRatio.m2000),type='b',col=1,pch=16)
par(mfrow=c(1,1))



# process -----------------------------------------------------------------
rm(list = ls())
counts <- read.delim(file="../data//counts-gene.txt",header = T,row.names = 1,skip = 1)
dim(counts)
colnames(counts[,c(1:5,9)])
counts <- counts[,-c(1:5,9)]
fun.cpm <- function(counts) counts/sum(counts)*10^6 
hist(log2(counts[,6]+.5))
counts$ifnarikba.dko.T8.filtered.bam  <- fun.cpm(counts$ifnarikba.dko.T8.filtered.bam)
sum(counts$ifnarikba.dko.T8.filtered.bam)

data_folder <- "~/Dropbox/Projects/DurationDecoding/data/"
ev <- new.env()
load(file=paste0(data_folder,"caRNA_RNAseq_kdeg.RData"),envir = ev)
load (file=paste0(data_folder,'avg_cpm_caRNA.RData'))
nfkbgenes <- ev$nfkbgenes; coldata_RNAseq_b1 <- ev$coldata_RNAseq_b1
cpm.avg.gene <- cpm.avg.gene[as.character(nfkbgenes),]
rownames(cpm.avg.gene) <- names(nfkbgenes)
all.equal(rownames(cpm.avg.gene), names(nfkbgenes))

plot(colSums(cpm.avg.gene))

idx <- which(substr(rownames(counts),1,18) %in% as.vector(nfkbgenes))
length(idx)
counts <- counts[idx,]
rownames(counts) <- substr(rownames(counts) ,1,18)
all.equal(rownames(counts),as.character(nfkbgenes))
counts<- counts[nfkbgenes,]
sum(counts$ifnarikba.dko.T8.filtered.bam)
head(cpm.avg.gene)
sum(cpm.avg.gene[,"D_T_8"])

pd <- data.frame(p= cpm.avg.gene[,"D_T_8"],
                 b3= counts$ifnarikba.dko.T8.filtered.bam)
pd <- log2(pd+.5)

plot(pd,main=paste0('cor=',signif( cor(pd,use = "na.or.complete")[1,2],2)))
abline(a=0,b=1,lwd=3,col=2)
cpm.avg.gene[,"D_T_8"] <- (cpm.avg.gene[,"D_T_8"]+ counts$ifnarikba.dko.T8.filtered.bam)/2
points(pd$p,log2(cpm.avg.gene[,"D_T_8"]+.5),pch=16)

write.csv(file='caRNA-cpm-avg.csv',cpm.avg.gene)


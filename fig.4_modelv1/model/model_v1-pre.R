## load data 
rm(list=ls())
load(file='mRNA-Fit-avg-sp-v1c-pre.Rdata')


# cmp two version of half-life --------------------------------------------
v1.hf <- read.csv(file='v4-hf-final.csv',row.names = 1)
genes <- read.csv(file='mRNA.cluster.csv',row.names = 2,stringsAsFactors = F)
all(rownames(v1.hf) %in% genes$gene)
v2.hf <- read.csv(file='Table.S.genome.hf.csv',row.names = 1,stringsAsFactors = F)
v2.hf <- v2.hf[rownames(v2.hf)%in% rownames(genes),]
v2.hf <- v2.hf[rownames(genes),]
rownames(v2.hf) <- rownames(genes)
v2.hf$Symbol <- genes$gene
v1.hf <- v1.hf[genes$gene,]
data.frame(round(10^v1.hf$LPS0),v2.hf$mean,v2.hf$Symbol)
sum(is.na(v2.hf$mean)) # 74 
sum(is.na(v1.hf$LPS0)) # 22

# construct kdeg ranges  --------------------------------------------------
ev <- new.env()
load("whole.genome.error.model.hf.RData",ev)

predict.rlm <- function (object, newdata = NULL, scale = NULL, ...){
  ## problems with using predict.lm are the scale and
  ## the QR decomp which has been done on down-weighted values.
  object$qr <- qr(sqrt(object$weights) * object$x)
  predict.lm(object, newdata = newdata, scale = object$s, ...)
}

final.hf <- data.frame(hf=round(10^v1.hf$LPS0),
                       #hf.sd = 10^predict.rlm(ev$rlm.fit,newdata = data.frame(x=v1.hf$LPS0)),
                       Symbol= rownames(v1.hf),stringsAsFactors = F)
for(i in 1:nrow(final.hf)) final.hf$hf[i] <- ifelse(is.na(v2.hf$mean[i]),final.hf$hf[i],v2.hf$mean[i])
final.hf$sd <- 10^predict.rlm(ev$rlm.fit,newdata = data.frame(x=log10(final.hf$hf)))
genes.kdeg <- log(2)/final.hf$hf;genes.kdeg.lb <- log(2)/(final.hf$hf+final.hf$sd);genes.kdeg.ub<-log(2)/(final.hf$hf-final.hf$sd)
names(genes.kdeg) <- names(genes.kdeg.lb) <- names(genes.kdeg.ub) <- final.hf$Symbol

# assign NA to max range 
idx.na <- which(is.na(genes.kdeg))
genes.kdeg[idx.na] <- median(genes.kdeg,na.rm = T)
genes.kdeg.lb[idx.na] <- min(genes.kdeg.lb,na.rm = T)
genes.kdeg.ub[idx.na] <- max(genes.kdeg.ub,na.rm = T)


# construct normalized counts  --------------------------------------------
genes.rpkm <- read.csv(file="mRNA.nfkbgene.rpkm.csv",stringsAsFactors = F,row.names = 1)
all.equal(rownames(genes.rpkm),names(genes.kdeg)) #TRUE 
genes.rpkm <- genes.rpkm[,c(1:14,22:35)]
data.frame(colnames(genes.rpkm),colnames(genes.normCount))
tmp <- cbind(t(apply(genes.rpkm[,1:14],1,function(x) x/max(x))),
             t(apply(genes.rpkm[,1:14 +14],1,function(x) x/max(x))))
all.equal(rownames(tmp),names(genes.kdeg)) #TRUE 

genes.normCount <- tmp 

# peak genes-------------------------------------------------------------------
genes.cat <- read.csv(file = "mRNA.cat.csv",stringsAsFactors = F)
table(genes.cat$cate)


# save the data -----------------------------------------------------------


save(file='mRNA-Fit-avg-sp-v1c-pre.Rdata',list=c("genes.kdeg","genes.kdeg.lb","genes.kdeg.ub",
                                                 "genes.normCount","genes.cat",
                                            "input.emsa.ifnarikba","idx.na",
                                            "input.emsa.ifnar", "genes.timepoints"
))


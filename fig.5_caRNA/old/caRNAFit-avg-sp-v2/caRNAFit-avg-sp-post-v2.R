rm(list=ls())
require(RColorBrewer)
require(pheatmap)

load(file = './caRNAFit-avg-sp-pre-v2.Rdata')
load(file = './caRNAFit-avg-sp-run-v2.Rdata')
source(file = 'caRNAFit-avg-sp-funs-v2.R')

nm <- 1:length(caFit)
for(i in 1:max(nm)) nm[i] <- caFit[[i]]$gene
names(caFit) <- nm

caFit.R2 <- sapply(1:length(caFit), calR2)
names(caFit.R2) <- nm

ord <- order(caFit.R2,decreasing = T) 

pdf(file='caRNAFit-avg-sp-v2.pdf')
hist(caFit.R2,breaks = c(min(caFit.R2),seq(0,1,length.out = 11)),freq = T,xlim = c(0,1))
lapply(caFit[ord], plotFitInd)
dev.off()




### correlation - ordering by caRNA time & mRNA time 
genes.timepoints
caRNA.time
tms <- seq(0,480,by = 30)
cor.mat <- rep(0,79)
for(i in 1:79){
  pd.caRNA <- approxfun(caRNA.time,caRNA.frac.avg[i,1:length(caRNA.time)])(tms)
  pd.mRNA <- approxfun(genes.timepoints,genes.normCount[i,1:length(genes.timepoints)])(tms)
  cor.mat[i] <- cor.test(pd.caRNA,pd.mRNA,method = "pearson")$estimate
} #LPS-contrl first 
names(cor.mat) <- rownames(caRNA.frac.avg)

which.max
all.equal(rownames(caRNA.frac.avg),rownames(genes.normCount))
peak.time <- data.frame(caRNA = apply(caRNA.frac.avg[,1:length(caRNA.time)], 1, which.max),
                        mRNA = apply(genes.normCount[,1:length(genes.timepoints)], 1, which.max),
                        gene = rownames(caRNA.frac.avg))
library(dplyr)
new.ord <- as.character(arrange(peak.time,caRNA,mRNA)$gene)

require(pheatmap)
dim()
pd <- cbind(caRNA.frac.avg[,1:length(caRNA.time)]/rowMax(caRNA.frac.avg[,1:length(caRNA.time)]),
            genes.normCount[,1:length(genes.timepoints)] )
timepoints <- union(caRNA.time,genes.timepoints)
timepoints <- timepoints[order(timepoints)]

pd.new<- sapply(1:nrow(pd),function(i)
  c(approxfun(caRNA.time,pd[i,1:length(caRNA.time)])(timepoints),
  approxfun(genes.timepoints,pd[i,(length(caRNA.time)+1):ncol(pd)])(timepoints)))
pd.new <- t(pd.new)
rownames(pd.new) <- rownames(pd)
pd <- pd.new

ann.row <- data.frame(cor =cut(cor.mat[new.ord],breaks = c(-1,seq(0,1,by = .2)))) 
ann.cols <- list(cor= colorRampPalette(colors = c("black","yellow"))(length(unique(ann.row$cor))))
rownames(ann.row) <- new.ord
names(ann.cols$cor) <- as.character(unique(ann.row$cor)[order(unique(ann.row$cor))])
pheatmap(pd[new.ord,],color = colorRampPalette(c("blue", "white", "red"))(50),
         scale = 'none',cluster_rows = F,cluster_cols = F,gaps_col = length(timepoints),#cellwidth = 12,cellheight = 4,
         annotation_row = ann.row,annotation_colors = ann.cols)

mean(cor.mat)
pheatmap(cor.mat[new.ord],scale = 'none',cluster_rows = F,cluster_cols = F,cellwidth = 12,cellheight = 4)

plot(cor.mat[new.ord],type = 'b',main='spearman rho, ctrl.lps')

plot(cor.mat[new.ord],type = 'b',main='pearson c, ctrl.lps')

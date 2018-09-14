rm(list = ls())
# half life data 
load(file='./caRNAFit-avg-sp-v6/nfkb-kdeg-final.Rdata')

# v1 data : emsa -> mRNA
tmp <- read.csv(file='./mRNA-Fit-avg-sp-v1c/result.csv',stringsAsFactors = F,row.names = 5)
v1.r2 <- tmp$r2; names(v1.r2) <- rownames(tmp)

# v2 data : emsa -> mRNA (kdeg(t))
tmp <- read.csv(file='./mRNA-Fit-avg-sp-v2/results.csv',stringsAsFactors = F,row.names = 4)
v2.r2 <- tmp$r2; names(v2.r2) <- rownames(tmp)


# v3 data: emasa -> caRNA 
tmp <- read.csv(file = './emsa-caRNA-v1/v1-cmp-N.csv')
v3.r2 <- tmp$r2; names(v3.r2) <- tmp$gene

# v4 data (caFit ) caRNA ->mRNA
v4.r2 <- read.csv(file='./caRNAFit-avg-sp-v6/caRNAFit-avg-sp-funs-v6-r2.csv',stringsAsFactors = F,row.names = 1)


genes.kdeg.isNA[,1]



# pdf 

require(pheatmap)
anno.row <- data.frame()
pd <- data.frame(basalKd = !genes.kdeg.isNA[names(v1.r2),1],
                 v1.r2 = v1.r2,
                 v2.r2 = v2.r2[names(v1.r2)],
                 v3.r2 = v3.r2[names(v1.r2)],
                 v4.r2 = v4.r2[names(v1.r2),1])
write.csv(file='result.csv',pd)


pd <- read.csv(file = 'result.csv',row.names = 1)
sp.avg <- read.csv(file = "../fig.1_wt/b1.sp.csv",row.names = 1)
pd <- (cbind(pd,sp.avg[rownames(pd),1:6]))



# read rnaseq data  -------------------------------------------------------
anno.row <- read.csv(file='result.csv',row.names = 1,stringsAsFactors = F)

pd.rnaseq.max <- read.csv(file='../fig.1_wt/data/maxscale.csv',row.names = 1,stringsAsFactors = F)
pd.rnaseq.clust <- read.csv(file = '../fig.1_wt/data/cluster.csv',stringsAsFactors = F)
all.equal(pd.rnaseq.clust$gene , rownames(pd.rnaseq.max))

pd.rnaseq.max <- pd.rnaseq.max[rownames(pd.rnaseq.max)%in% rownames(anno.row),]
anno.row <- anno.row[rownames(anno.row)%in% rownames(pd.rnaseq.max),]
anno.row <- anno.row[rownames(pd.rnaseq.max), ]
pd.rnaseq.clust <- subset(pd.rnaseq.clust,gene %in% rownames(pd.rnaseq.max))
all.equal((pd.rnaseq.clust$gene),rownames(pd.rnaseq.max))

anno.row[is.na(anno.row)] <- 0 
anno.row$basalKd <- as.factor(anno.row$basalKd)

require(pheatmap)
pheatmap(pd.rnaseq.max,
         cols <- colorRampPalette(c( "white", "firebrick1"))(5),
         cluster_rows = F,cluster_cols = F,
        gaps_row = sapply(2:8, function(k) max(which(pd.rnaseq.clust$cluster==k))),
        gaps_col = seq(0,ncol(pd.rnaseq.max),by=7),
        annotation_row =  anno.row
       )

rm(list=ls())
require(RColorBrewer)
require(pheatmap)
sp.caRNA <- read.csv(file='./data/caRNA-max-avg-sp.csv',header = T,
                     row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../fig.1_wt/data/sp.fraction.maxscale.csv",header = T,
                    row.names = 1,stringsAsFactors = F)




library(biomaRt)
ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
tmp.2 <- rownames(sp.mRNA)
tmp <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
             values = grep('ENSMUSG',tmp.2,value = T), mart = ensembl)
rownames(tmp) <- tmp$ensembl_gene_id
tmp.2[grep('ENSMUSG',tmp.2)] <- tmp[tmp.2[grep('ENSMUSG',tmp.2)],"mgi_symbol"]
rownames(sp.mRNA) <-  tmp.2 

sp.caRNA <- sp.caRNA[rownames(sp.mRNA),]
all.equal(rownames(sp.caRNA),rownames(sp.mRNA))
sum(is.na(sp.caRNA))


pdf(file='./fig/cmp-caRNA-mRNA-avgsp.pdf')
par(mfrow=c(2,2))
pd <- data.frame(tlsp.mRNA.ctrl = sp.mRNA$,
                 tlsp.caRNA.ctrl=sp.caRNA$wt,
                tlsp.mRNA.ko=sp.mRNA$sp.LT.avg.ko,
                tlsp.caRNA.ko= sp.caRNA$mt)

pd <- log2(pd)                
pd[125,2] <- NA; pd[125,2] <- min(pd,na.rm = T)
xylims <- c(min(pd),max(pd))                 
plot(pd[,1:2],pch=21,bg="grey",xlim=xylims,ylim=xylims,
     main=paste("cor =",signif(cor(pd[,1],y=pd[,2]),2)))
abline(a=0,b=1)

abline(v=c(-.5,.5))
abline(h=c(-.5,.5))

plot(pd[,3:4],pch=21,bg="grey",xlim=xylims,ylim=xylims,
     main=paste("cor =",signif(cor(pd[,3],y=pd[,4]),2)))
abline(a=0,b=1)
cor(pd[,3],y=pd[,4])
abline(v=c(-.5,.5))
abline(h=c(-.5,.5))




plot(pd[,c(2,4)],pch=21,bg=ifelse(pd[,1]<=-.5 & pd[,3]> -.5 & pd[,3]<.5,'coral1',"grey"),
     xlim=xylims,ylim=xylims,
     main=paste0('cor=',signif(cor(pd[,2],pd[,4]),2)))
abline(a=0,b=1)
abline(v=c(-.5,.5))
abline(h=c(-.5,.5))

plot(pd[,c(1,3)],pch=21,bg=ifelse(pd[,1]<=-.5 & pd[,3]> -.5 & pd[,3]<.5,'coral1',"grey"),
     xlim=xylims,ylim=xylims,
     main=paste0('cor=',signif(cor(pd[,1],pd[,3]),2)))
abline(a=0,b=1)
abline(v=c(-.5,.5))
abline(h=c(-.5,.5))

par(mfrow=c(1,1))
range(pd)
bks <- c(-1.5,-1,-0.5,0,.5,1,1.5)
cols <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                    "RdYlBu")))(length(bks)-1)
cols[4]<- cols[3]
pheatmap(pd, breaks = bks,color = cols,cutree_rows = 8, #kmeans_k = 8,
         scale = 'none',cluster_cols = F,show_rownames = F)

dev.off()

# plot caRNA seq (max scale)----------------------------------------------------------
rm(list=ls());

cpm.avg.gene <- read.csv(file='../data/caRNA-cpm-avg.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
load("../fig.1_wt/data/kcluster_b1.Rdata")
all.equal(rownames(cpm.avg.gene),k_ord$gene)

seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.5_caRNA/subfigs/'


# plot data 

pd <-cpm.avg.gene;
fun.scale <- function(pd){
  idx.I <- (substr(colnames(pd),1,1)== "I")
  rmax.I <- apply(pd[,idx.I],1,max)
  idx.D <- (substr(colnames(pd),1,1)== "D")
  rmax.D <- apply(pd[,idx.D],1,max)
  return(cbind(pd[,idx.I]/rmax.I,pd[,idx.D]/rmax.D))
}

fun.scale.2 <- function(pd){
  idx.I <- (substr(colnames(pd),1,1)== "I")
  rmax.I <- apply(pd[,idx.I],1,max)
  idx.D <- (substr(colnames(pd),1,1)== "D")
  rmax.D <- apply(pd[,idx.D],1,max)
  return(cbind(pd[,idx.I]/rmax.I,pd[,idx.D]/rmax.I))
}

pd <- fun.scale(cpm.avg.gene) # scale to each genotpye
pd <- fun.scale.2(cpm.avg.gene) # scale to wt only 

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
tmp.2 <- rownames(pd)
tmp <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
             values = grep('ENSMUSG',tmp.2,value = T), mart = ensembl)
rownames(tmp) <- tmp$ensembl_gene_id
tmp.2[grep('ENSMUSG',tmp.2)] <- tmp[tmp.2[grep('ENSMUSG',tmp.2)],"mgi_symbol"]
rownames(pd) <-  tmp.2 


head(pd)
write.csv(file='../data/caRNA-cpm-max-scale.csv',pd)

#pd[pd>1] <- 1
dim(pd)
dim(k_ord)
pd <- pd[k_ord$gene,]
bks <- seq(0,1,by = 0.2)
cols <-colorRampPalette(c("azure","blue"))(6)

require(pheatmap)
pheatmap(pd,cluster_cols = F,cluster_rows = F,scale = 'none',show_rownames = F,color = cols,
         breaks = bks)

#cols <- colorRampPalette(c("black","firebrick1"))(6)
source('../fig.1_wt/auxilary_functions.R')
fbks <- seq(0,1,length.out = 6)
sapply(1:2, 
       function(i){
         #i <-2
         setEPS()
         postscript(paste0(fig_folder,"caRNA_max",i,".eps"),onefile = F,width = 6,height = 6)
         pheatmap.2(pd=pd[,((i-1)*12+1):(i*12)], bks = bks,gaps_col = c(6),cwid = 6,
                    legend = F,cols = cols)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"caRNA_max_scale.eps"))


# plot caRNA zscale -------------------------------------------------------
pd <- cpm.avg.gene
?scale

fun.zscale <- function(pd){
  idx.I <- (substr(colnames(pd),1,1)== "I")
  idx.D <- (substr(colnames(pd),1,1)== "D")

  return(cbind(t(scale(t(pd[,idx.I]))),
                 t(scale(t(pd[,idx.D])))))
}

pd <- fun.zscale(cpm.avg.gene)
pd[pd>3] <- 3; pd[pd< -3] <- -3;pd <- pd[k_ord$gene,]

cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(pd),max(pd),length.out = 21)

pheatmap(pd,color = cols,breaks = bks,
         cluster_rows = F,cluster_cols = F)

sapply(1:2, 
       function(i){
         #i <-2
         setEPS()
         postscript(paste0(fig_folder,"caRNA_zscale",i,".eps"),onefile = F,width = 6,height = 6)
         pheatmap.2(pd=pd[,((i-1)*12+1):(i*12)], bks = bks,gaps_col = c(6),cwid = 6,
                    legend = F,cols = cols)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"caRNA_z_scale.eps"))


# caRNA zscale sp ---------------------------------------------------------

rsums <- sapply(0:3, function(i)
  rowSums(pd[,(i*6+1):((i+1)*6)]))

deltaAuc<-sapply(c(1,3), function(i) (rsums[,i+1] - rsums[,i]) )# TNF -LPS 
deltaAuc <- deltaAuc/6
range(deltaAuc);dim(deltaAuc)
colnames(deltaAuc) <- c('ctrl.tlsp','ko.tlsp')
write.csv(file='caRNA-zscale-sp.csv',deltaAuc)

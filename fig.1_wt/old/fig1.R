rm(list=ls())
dataFolder<-"~/Dropbox/Projects/DurationDecoding/data/RNAseq/processed_data/"
load(file=paste0(dataFolder,'lrpkmD_clean_data_all.Rdata'))
dim(lrpkm_all)
require(RColorBrewer)
require(pheatmap)
source('./auxilary_functions.R')

fig_folder <-'/Users/zhangcheng/Dropbox/Projects/DurationDecoding/figure/tmp_figs/'
k_ord <- read.csv(file=paste0(fig_folder,'clustering_b1.csv'),stringsAsFactors = F)

all.equal(levels(as.factor(rownames(lrpkm_all))),
          levels(as.factor(k_ord$gene)))
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

#-----plot z rpkm (individual genoScaling)-----
## batch 1 & batch2 

plot_data1 <- 2^lrpkm_all[,1:42]
rord <-c(1:7,15:21,29:35,8:14,22:28,36:42)
plot_data1<-plot_data1[k_ord$gene,rord]

plot_data2 <- 2^lrpkm_all[,43:ncol(lrpkm_all)]
rord <-c(1:8,17:24,33:40,9:16,25:32,41:48)
plot_data2<- plot_data2[k_ord$gene,rord]
plot_data2<- plot_data2[,-grep(6,colnames(plot_data2))]
 
plot_data <- cbind(plot_data1,plot_data2); 
rm(plot_data1);rm(plot_data2)
colnames(plot_data)

## individual geno scaling 

plot_data2 <- cbind(t(scale(t(plot_data[,1:21]))),
                    t(scale(t(plot_data[,22:42]))),
                    t(scale(t(plot_data[,43:63]))),
                    t(scale(t(plot_data[,64:84]))))

plot_data2[plot_data2>3] <- 3;plot_data2[plot_data2< -3] <- -3
colnames(plot_data2)

fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(plot_data2),max(plot_data2),length.out = 21)
sapply(1:4, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"z_rpkm_genoScale_",i,".eps"),onefile = F,width = 6,height = 6)
         pd <- plot_data2[,((i-1)*21+1):(i*21)]
         pheatmap.2(pd=pd,bks=bks,cols = cols,cwid = 4,
                    legend=F,gaps_col = grep('8',colnames(pd)))
         dev.off()
       })
plotLegend(cols = cols,
           bks = bks,
           fnames = paste0(fig_folder,"z_rpkm_batchScale_lg.eps"))

# AUC of z score -------
rsums <- sapply(0:11, function(i)
  rowSums(plot_data2[,(i*7+1):((i+1)*7)]))

cols <- colorRampPalette(c( "white", "blueviolet"))(5)
pd <- rsums + 7*abs(min(rsums)) #shift the auc to make them positive
setEPS()
postscript(paste0(fig_folder,"auc_zrpkm.eps"),onefile = F,width = 6,height = 6)
par(mar=rep(0,4)+0.1)
bks <- seq(round(min(pd)),round(max(pd))+1,length.out = 6)
pheatmap.2(pd=pd,bks=bks,cols=cols,gaps_col = c(3,6,9))
dev.off()

deltaAuc<-sapply(c(1,4,7,10), 
       function(i){
         deltaAuc <- pd[,i+1] - pd[,i] # TNF -LPS 
         #deltaAuc[deltaAuc>0] <- 0 
         #deltaAuc <- deltaAuc/pd[,i]
         #deltaAuc <- pd[,i+1] /pd[,i] # TNF/LPS 
         #deltaAuc[deltaAuc>1] <- 1
         deltaAuc
       })

bks <- seq(min(deltaAuc),max(deltaAuc),length.out = 6)
cols <- colorRampPalette(c( brewer.pal(n=11,name = 'PiYG')[6], "blueviolet"))(5)
cols <- rev(cols)
#pheatmap.2(pd=deltaAuc,bks=bks,cols=cols,legend=T)
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
sapply(1:4,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"dAUC_zrpkm_genoScale",i,".eps"),
                    onefile = F,width = 6,height = 6)
         pheatmap.2(pd=deltaAuc[,i],bks=bks,cols=cols,legend=F)  
         dev.off()
       })


plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"dAUC_zrpkm_genoScale_lg.eps"))
signif(bks,2)

#-----plot z rpkm (batch Scaling)-----
## batch 1 & batch2 

plot_data1 <- 2^lrpkm_all[,1:42]
rord <-c(1:7,15:21,29:35,8:14,22:28,36:42)
plot_data1<-plot_data1[k_ord$gene,rord]

plot_data2 <- 2^lrpkm_all[,43:ncol(lrpkm_all)]
rord <-c(1:8,17:24,33:40,9:16,25:32,41:48)
plot_data2<- plot_data2[k_ord$gene,rord]
plot_data2<- plot_data2[,-grep(6,colnames(plot_data2))]

plot_data <- cbind(plot_data1,plot_data2); 
rm(plot_data1);rm(plot_data2)
colnames(plot_data)

## individual batch scaling 

plot_data2 <- cbind(t(scale(t(plot_data[,1:42]))),
                    t(scale(t(plot_data[,43:84]))))

plot_data2[plot_data2>3] <- 3;plot_data2[plot_data2< -3] <- -3
colnames(plot_data2)

fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(plot_data2),max(plot_data2),length.out = 21)
sapply(1:4, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"z_rpkm_batchScale_",i,".eps"),onefile = F,width = 6,height = 6)
         pd <- plot_data2[,((i-1)*21+1):(i*21)]
         pheatmap.2(pd=pd,bks=bks,cols = cols,cwid = 4,
                    legend=F,gaps_col = grep('8',colnames(pd)))
         dev.off()
       })
plotLegend(cols = cols,
           bks = bks,
           fnames = paste0(fig_folder,"z_rpkm_batchScale_lg.eps"))


# AUC of z score (batch Scaling)-------
rsums <- sapply(0:11, function(i)
  rowSums(plot_data2[,(i*7+1):((i+1)*7)]))

cols <- colorRampPalette(c( "white", "blueviolet"))(5)
pd <- rsums + 7*abs(min(rsums))
bks <- seq(round(min(pd)),round(max(pd))+1,length.out = 6)

setEPS()
postscript(paste0(fig_folder,"auc_zrpkm_batchScale.eps"),onefile = F,width = 6,height = 6)
par(mar=rep(0,4)+0.1)
pheatmap.2(pd=pd,bks=bks,cols=cols,gaps_col = c(3,6,9))
dev.off()

deltaAuc<-sapply(c(1,4,7,10), 
                 function(i){
                   deltaAuc <- pd[,i+1] - pd[,i] # TNF -LPS 
                   deltaAuc[deltaAuc>0] <- 0 
                   deltaAuc <- deltaAuc/pd[,i]
                   deltaAuc
                 })

bks <- seq(min(deltaAuc),max(deltaAuc),length.out = 6)
cols <- colorRampPalette(c( brewer.pal(n=11,name = 'PiYG')[6], "blueviolet"))(5)
cols <- rev(cols)
#pheatmap.2(pd=deltaAuc,bks=bks,cols=cols,legend=T)
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
sapply(1:4,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"dAUC_zrpkm_batchScale",i,".eps"),
                    onefile = F,width = 6,height = 6)
         pheatmap.2(pd=deltaAuc[,i],bks=bks,cols=cols,legend=F)  
         dev.off()
       })


plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"dAUC_zrpkm_batchScale_lg.eps"))
signif(bks,2)


#-----plot fraction of max (genoScale)-----
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
plot_data1 <- 2^lrpkm_all[,1:42]
rord <-c(1:7,15:21,29:35,8:14,22:28,36:42)
plot_data1<- plot_data1[k_ord$gene,rord]

plot_data2 <- 2^lrpkm_all[,43:ncol(lrpkm_all)]
rord <-c(1:8,17:24,33:40,9:16,25:32,41:48)
plot_data2<- plot_data2[k_ord$gene,rord]
plot_data2<- plot_data2[,-grep(6,colnames(plot_data2))]

plot_data <- cbind(plot_data1,plot_data2)
rm(plot_data1);rm(plot_data2)

for(i in 1:4){
  cls <- ((i-1)*21+1):(i*21)
  rmx <- apply(plot_data[,cls], 1, max)  
  plot_data[,cls]<- plot_data[,cls]/rmx
}
plot_data[plot_data>1] <- 1

require(pheatmap)
cols <- colorRampPalette(c("black","firebrick1"))(6)
bks <- seq(0,1,length.out = 6)
sapply(1:4, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"z_maxFrac_genoScale_",i,".eps"),onefile = F,width = 6,height = 6)
         pd <- plot_data[,((i-1)*21+1):(i*21)]
         pheatmap.2(pd=pd, bks = bks,gaps_col = c(7,14),cwid = 4,
                  legend = F,cols = cols)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"z_maxFrac_genoScale_lg.eps"))

# AUC of max frac-----
rsums <- sapply(0:11, function(i)
  rowSums(plot_data[,(i*7+1):((i+1)*7)]))
range(rsums)
bks <- seq(0,5,length.out = 6)
cols <-colorRampPalette(c(  "deepskyblue","azure"))(5)

postscript(paste0(fig_folder,"auc_maxFrac_batchScale.eps"),
           onefile = F,width = 6,height = 6)
par(mar=rep(0,4)+0.1)
pheatmap.2(pd=rsums,bks = bks,legend = F,
         gaps_col = c(3,6,9),cols = cols)
dev.off()
plotLegend(cols = cols,bks = bks,
           fnames = paste0(fig_folder,"auc_maxFrac_batchScale_lg.eps"))


 
deltaAuc<-sapply(c(1,4,7,10), 
                 function(i){
                   deltaAuc <- rsums[,i+1]/rsums[,i] # TNF -LPS 
                   deltaAuc[deltaAuc>1] <- 1 
                   deltaAuc
                 })
range(deltaAuc)
bks <- seq(0,1,length.out = 6)
cols <-colorRampPalette(c(  "darkblue","azure"))(5)
pheatmap.2(pd=deltaAuc,bks=bks,cols=cols,legend=T)  
sapply(1:4,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"dAUC_maxFrac_genoScale",i,".eps"),
                    onefile = F,width = 6,height = 6)
         pheatmap.2(pd=deltaAuc[,i],bks=bks,cols=cols,legend=F)  
         dev.off()
       })

plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"dAUC_maxFrac_genoScale_lg.eps"))
signif(bks,2)

#-----plot fraction of max (batchScale)-----
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
plot_data1 <- 2^lrpkm_all[,1:42]
rord <-c(1:7,15:21,29:35,8:14,22:28,36:42)
plot_data1<- plot_data1[k_ord$gene,rord]

plot_data2 <- 2^lrpkm_all[,43:ncol(lrpkm_all)]
rord <-c(1:8,17:24,33:40,9:16,25:32,41:48)
plot_data2<- plot_data2[k_ord$gene,rord]
plot_data2<- plot_data2[,-grep(6,colnames(plot_data2))]

plot_data <- cbind(plot_data1,plot_data2)
rm(plot_data1);rm(plot_data2)

for(i in 1:2){
  cls <- ((i-1)*42+1):(i*42)
  rmx <- apply(plot_data[,cls], 1, max)  
  plot_data[,cls]<- plot_data[,cls]/rmx
}
plot_data[plot_data>1] <- 1

require(pheatmap)
cols <- colorRampPalette(c("black","firebrick1"))(6)
bks <- seq(0,1,length.out = 6)
sapply(1:4, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"z_maxFrac_bScale_",i,".eps"),onefile = F,width = 6,height = 6)
         pd <- plot_data[,((i-1)*21+1):(i*21)]
         pheatmap.2(pd=pd, bks = bks,gaps_col = c(7,14),cwid = 4,
                    legend = F,cols = cols)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"z_maxFrac_bScale_lg.eps"))

# AUC of max frac (batch scaleing)-----
rsums <- sapply(0:11, function(i)
  rowSums(plot_data[,(i*7+1):((i+1)*7)]))
range(rsums)
bks <- seq(0,5,length.out = 6)
cols <-colorRampPalette(c(  "deepskyblue","azure"))(5)

postscript(paste0(fig_folder,"auc_maxFrac_batchScale.eps"),
           onefile = F,width = 6,height = 6)
par(mar=rep(0,4)+0.1)
pheatmap.2(pd=rsums,bks = bks,legend = F,
           gaps_col = c(3,6,9),cols = cols)
dev.off()
plotLegend(cols = cols,bks = bks,
           fnames = paste0(fig_folder,"auc_maxFrac_batchScale_lg.eps"))



deltaAuc<-sapply(c(1,4,7,10), 
                 function(i){
                   deltaAuc <- rsums[,i+1]/rsums[,i] # TNF -LPS 
                   deltaAuc[deltaAuc>1] <- 1 
                   deltaAuc
                 })
range(deltaAuc)
bks <- seq(0,1,length.out = 6)
cols <-colorRampPalette(c(  "darkblue","azure"))(5)
pheatmap.2(pd=deltaAuc,bks=bks,cols=cols,legend=T)  
sapply(1:4,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"dAUC_maxFrac_genoScale",i,".eps"),
                    onefile = F,width = 6,height = 6)
         pheatmap.2(pd=deltaAuc[,i],bks=bks,cols=cols,legend=F)  
         dev.off()
       })

plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"dAUC_maxFrac_genoScale_lg.eps"))
signif(bks,2)





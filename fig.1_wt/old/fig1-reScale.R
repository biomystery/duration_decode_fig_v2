rm(list=ls())

require(RColorBrewer)
require(pheatmap)
source('./auxilary_functions.R')

load(file=paste0('./data/lrpkmD_clean_data_all.Rdata'))
load(file="./data/kcluster_b1.Rdata") # k_ord
all.equal(levels(as.factor(rownames(lrpkm_all))),
          levels(as.factor(k_ord$gene)))
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])


#-----plot fraction of max (genoScale)-----
#fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'
fig_folder <- './'
plot_data1 <- 2^lrpkm_all[,1:42]
rord <-c(1:7,15:21,29:35,8:14,22:28,36:42)
plot_data1<- plot_data1[k_ord$gene,rord]

plot_data2 <- 2^lrpkm_all[,43:ncol(lrpkm_all)]
rord <-c(1:8,17:24,33:40,9:16,25:32,41:48)
plot_data2<- plot_data2[k_ord$gene,rord]
plot_data2<- plot_data2[,-grep(6,colnames(plot_data2))]

plot_data <- cbind(plot_data1,plot_data2)
rm(plot_data1);
plot_data2 <- plot_data
eps <- 1e-20 #avoid 0 
fun.rescale <- function(x) (x-min(x)+eps)/(max(x)-min(x))
for(i in 1:4)
  plot_data2[,((i-1)*21+1):(i*21)]<- t(apply(plot_data[,((i-1)*21+1):(i*21)], 1, fun.rescale))

range(plot_data2);colnames(plot_data2)
dim(plot_data2)
write.csv(file='./rescaled.counts.csv',plot_data2)
fig_folder <- "./"


require(pheatmap)
plot_data <- plot_data2
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


# AUC of z score -------
rsums <- sapply(0:11, function(i)
  rowSums(plot_data2[,(i*7+1):((i+1)*7)]))


deltaAuc<-sapply(c(1,4,7,10), function(i) (rsums[,i+1]/rsums[,i]) )# TNF -LPS 
deltaAuc.2 <-sapply(c(1,4,7,10), function(i) (rsums[,i+2]/rsums[,i]) )# IL1 -LPS 
deltaAuc.3 <-sapply(c(1,4,7,10), function(i) (rsums[,i+1]/rsums[,i+2]) )# TNF - IL1


deltaAuc.all <- data.frame(
  b1.TL.ctrl = deltaAuc[,1],
  b1.IL.ctrl = deltaAuc.2[,1],
  b1.TI.ctrl = deltaAuc.3[,1],
  b1.TL.ko = deltaAuc[,2],
  b1.IL.ko = deltaAuc.2[,2],
  b1.TI.ko = deltaAuc.3[,2],
  b2.TL.ctrl = deltaAuc[,3],
  b2.IL.ctrl = deltaAuc.2[,3],
  b2.TI.ctrl = deltaAuc.3[,3],
  b2.TL.ko = deltaAuc[,4],
  b2.IL.ko = deltaAuc.2[,4],
  b2.TI.ko = deltaAuc.3[,4]
) 

write.csv(file="sp.log2fraction.rescaled.csv",log2(deltaAuc.all))


# plot deltaAuc -----------------------------------------------------------


fig_folder <- './'#~/Dropbox/Projects/DurationDecoding/figure/Fig.1/subfigs/'

pd <- log2(deltaAuc)
range(pd)
#bks <- c(0,0.25,0.5,.75,1,1.25,1.5,2,3.5)
bks <- c(-8,-4,-1,0,1,1.5,2)
pd[pd< -8] <- -8
#bks <- c(0,0.25,0.5,1,1.5,2,3.5)
cols <-colorRampPalette(c(  "darkblue","azure","darkorange"))(length(bks)-1)
#cols[4]<- cols[5]
cols[3]<- cols[4] <- "azure"
pheatmap.2(pd=pd,bks=bks,cols=cols,legend=T)  
sapply(1:4,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"dAUC_maxFrac_genoScale",i,".eps"),
                    onefile = F,width = 6,height = 6)
         pheatmap.2(pd=pd[,i],bks=bks,cols=cols,legend=F)  
         dev.off()
       })

plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"dAUC_maxFrac_genoScale_lg.eps"))
signif(bks,2)


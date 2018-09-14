source('../auxilary_functions.R')
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.2/subfigs/'
require(reshape)
require(ggplot2)
require(plyr)
require(pheatmap)

# cluster vs. half-life  --------------------------------------------------
genes.hf <- read.csv(file="../data/v4-hf-final.csv",row.names = 1,stringsAsFactors = F)
all.equal(rownames(sp.mat[kord.sp.idx,]),names(kord.sp.idx))
all.equal(rownames(pd.scale[kord.sp.idx,]),names(kord.sp.idx))
sum(names(kord.sp.idx) %in% rownames(genes.hf)) ==length(kord.sp.idx) # check if all have hf
pd.hf <- data.frame(hf=genes.hf[names(kord.sp.idx),"LPS0"],
                    clust=c(rep(1,kord.sp.gaps[1]),rep(2,kord.sp.gaps[2]-kord.sp.gaps[1]),
                            rep(3,kord.sp.gaps[3]-kord.sp.gaps[2]),
                            rep(4,kord.sp.gaps[4]-kord.sp.gaps[3])))
ggplot(pd.hf,aes(factor(clust),hf)) + geom_boxplot(aes(fill=clust))

# fig2b RNAseq heatmap (mutant vs ctl.) -------------------------------------------------------------------

### sp plot 
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
sp.glist <- (getSplist(sp.mat = sp.mat,sp.th = 0.5))$LT.L
sp.mat <- sp.mat[sp.glist,1:6]
sp.mat[sp.mat>3] <- 3; sp.mat[sp.mat < -3] <- -3 

set.seed(10)
kord.sp <- kmeans(sp.mat[,c(1,4)],nstart = 25,iter.max = 4000,centers = 4)
r_kord = c(2,4,1,3)
kord.sp.idx <- unlist(sapply(r_kord,function(x) which(kord.sp$cluster==x)))
kord.sp.gaps <- cumsum(sapply(r_kord,function(x) sum(kord.sp$cluster==x)))
seps <- kord.sp.gaps + 1 # for pheatmap.2 
cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)

sapply(c(1,4),
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"subfig2c_peak_sp_",i,".eps"),
                    onefile = F,width = 1.5,height = 6)
         pheatmap.2(pd=sp.mat[kord.sp.idx,i],bks=bks,cols=cols,legend=F,
                    cwid = 10,border_color=NA)#,gaps_col=c(1,2))  
         dev.off()
       })


### expression plot 
pd <- read.csv(file = '../data/mRNA.nfkbgene.rpkm.csv',stringsAsFactors = F,row.names = 1,header = T)

# include IL1 
pd.scale <- cbind(t(scale(t(pd[,1:21]))),
                  t(scale(t(pd[,22:42]))),
                  t(scale(t(pd[,43:63]))),
                  t(scale(t(pd[,64:84]))))
# not include IL1
pd.scale <- cbind(t(scale(t(pd[,1:14]))),
                  t(scale(t(pd[,22:35]))))

pd.scale[pd.scale>3] <- 3; pd.scale[pd.scale< -3] <- -3
pd.scale <- pd.scale[sp.glist,]

cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(pd.scale),max(pd.scale),length.out = 21)
sapply(c(1,15),  #c(1,22) for scaling including IL1
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"subfig2c_zscore_hm_",i,".eps"),onefile = F,width = 2,height = 6)
         dat <- pd.scale[kord.sp.idx,i:(i+13)]
         pheatmap.2(pd=dat,bks=bks,cols = cols,cwid = 5,
                    legend=F,gaps_col = grep('8',colnames(dat)),border_color=NA)
         dev.off()
       })


# old fig2F venn of sp-------------------------------------------------------------------
## venn plot 
require(limma)

glist <- data.frame(LPS.specific.ctrl. = pd$b1.LT.ctrl>= sp.th,
                    TNF.specific.ctrl. = pd$b1.LT.ctrl <= -sp.th,
                    dynamic.dependent = ((pd$b1.LT.ko -pd$b1.LT.ctrl) <= -sp.th) )
setEPS()
postscript(file=paste0(fig_folder,'subfig2f_venn.eps'),width = 3,height = 3)
par(mar=rep(0,4))
vennDiagram(vennCounts(glist),cex = .75,mar = rep(0,4),
            circle.col = c('red','green','blue'))
dev.off()

# old fig2F scatter of sp-------------------------------------------------------------------
sp.th <- 0.5
# 
sp.mat <- read.csv(file='../data/sp.avg.csv',row.names = 1,stringsAsFactors = F)
sp.mat.log2 <- log2(sp.mat)
colnames(sp.mat.log2)
pd.1<- pd <- sp.mat.log2[,c(1,4)]
pd[pd< -3] <- -3 ; pd[pd > 3] <- 3
head(pd)

setEPS()
postscript(file=paste0(fig_folder,'subfig2f.eps'),width = 2.5,height = 2.5)
par(mar=rep(0,4))
plot(pd,xlim=c(-3,3),ylim=c(-3,3),pch=21,
     xlab='LvsT.sp.ctrl',ylab='LvsT.sp.mutant',
     bg=ifelse(pd[,1]>= sp.th  & pd[,1] -pd[,2]>=sp.th, 'blue','grey'))
axis(2,tck=0.02);axis(1,tck=0.02)
abline(a=0,b=1,col='grey80')
abline(a = -sp.th,b=1,col='blue',lty=2)
abline(v=c(-.5,.5),col='grey80'); abline(h=c(-.5,.5),col='grey80')
dev.off()

setEPS()
postscript(file=paste0(fig_folder,'subfig2f_a.eps'),width = 2.5,height = 2.5)
par(mar=rep(0,4))
plot(pd,xlim=c(-3,3),ylim=c(-3,3),pch=21,
     xlab='LvsT.sp.ctrl',ylab='LvsT.sp.mutant',
     bg=ifelse(pd[,1]>= sp.th , col.map['LPS'],'grey'))
abline(v=c(.5),col=col.map['LPS']); 
grid()
dev.off()

identify(pd,labels = rownames(pd))
# fig2b RNAseq heatmap (mutant only) -------------------------------------------------------------------
ev <- new.env()
load(file = '../data/lrpkmD_clean_data_all.Rdata',ev)
source(file = '../auxilary_functions.R')
require(pheatmap)

pd <- ev$lrpkm_all[,!grepl('[.]ifnr[.]',colnames(ev$lrpkm_all))]
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene)
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

pd.scale <- cbind(t(scale(t(pd[,1:21]))),
                  t(scale(t(pd[,22:42]))))

pd.scale[pd.scale>3] <- 3; pd.scale[pd.scale< -3] <- -3
cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(pd.scale),max(pd.scale),length.out = 21)
sapply(3:4, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"z_rpkm_genoScale2_",i,".eps"),onefile = F,width = 2,height = 6)
         dat <- pd.scale[,((i-3)*14+1):((i-2)*14)]
         pheatmap.2(pd=dat,bks=bks,cols = cols,cwid = 5,
                    legend=F,gaps_col = grep('8',colnames(dat)))
         dev.off()
       })


# fig2Sxx  ChipSeq locations  -------------------------------------------------------------------
pd <- read.csv(file = '../data/mRNA.nfkbgene.rpkm.csv')
par.windowsize <- 10000

# gene list 
pd <- pd[,!grepl('[.]ifnr[.]',colnames(ev$lrpkm_all))]
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene)
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

# annotations 
genes.anno <- read.delim(file='~/Dropbox/Projects/DurationDecoding/data/RNAseq/MEF_count_set1/counts-roberto.txt',
                         stringsAsFactors = F,skip = 1)
genes.anno <- genes.anno[,1:6]; genes.anno$Geneid <- substr(genes.anno$Geneid,1,18)
rownames(genes.anno) <- genes.anno$Geneid; genes.anno$Geneid <- NULL 

# read peak bed  
rela.genes <- k_ord$gene_id
peak.bed.raw <- read.delim(file = '../data/chip-seq_b1b2_overlap_narrowpeak.bed',header = F,
                           stringsAsFactors = F)
colnames(peak.bed.raw) <- c("Chr",'Start',"End")
rela.genes.anno <- genes.anno[rela.genes,]
rela.genes.anno$TSS <- unlist(lapply(rela.genes.anno$Start,
                                     function(x) as.numeric(unlist(strsplit(x,split = ";"))[1])))
rela.genes.anno$Chr <- unlist(lapply(rela.genes.anno$Chr,
                                     function(x) unlist(strsplit(x,split = ";"))[1]))

# reduce the bed 
chr.common <- intersect(unique(rela.genes.anno$Chr),unique(peak.bed.raw$Chr))
peak.bed.raw <- subset(peak.bed.raw, Chr %in% chr.common)
peak.bed.raw$length <- peak.bed.raw$End -  peak.bed.raw$Start + 1

# run find peaks in the windows
system.time(rela.peak.mat <- getPeakMat(par.windowsize = par.windowsize))

# plot the results  
png(file=paste0("rela.peak.mat.",par.windowsize,'.png'))
require(pheatmap)
pheatmap(rela.peak.mat,scale = 'none',
         show_rownames = F,color = c("grey90","darkmagenta"),show_colnames = F,
         cluster_rows = F,cluster_cols = F,legend = F,
         gaps_row = sapply(1:8, function(i) max(which(k_ord$cluster==i))),
         gaps_col = (ncol(rela.peak.mat-1)/2))
dev.off()
write.csv(file=paste0('rela.peak.mat.',par.windowsize,'.csv')
          ,rela.peak.mat)

# check peak number
any.peak <- apply(rela.peak.mat==1,1,any)
sum(any.peak)

write.csv(file='../data/cluster.csv',data.frame(k_ord,chip.seq.10k.peak=any.peak))

# fig2S-B -------------------------------------------------------------------
# add a layer of peaks into sp scatter
sp.th <- 0.5

# scatter of sp
sp.mat <- read.csv(file='../data/sp.avg.csv',row.names = 1,stringsAsFactors = F)
sp.mat.log2 <- log2(sp.mat)
colnames(sp.mat.log2)
all.equal(rownames(sp.mat),k_ord$gene)
pd.1<- pd <- sp.mat.log2[,c(1,4)]
pd[pd< -3] <- -3 ; pd[pd > 3] <- 3
head(pd)
#par(mar=rep(0,4))
plot(pd,xlim=c(-3,3),ylim=c(-3,3),pch=ifelse(pd.1[,2]>3,24,
                                             ifelse(pd.1[,2]< -3,25,
                                                    ifelse(pd.1[1,]>3,'>',
                                                           ifelse(pd.1[1,]< -3, '<',21)))),
     #bg=ifelse(pd[,1]>= sp.th & pd[,2] <sp.th & pd[,2]>-sp.th & pd[,1] -pd[,2]>=sp.th, 'blue','grey'))
     bg=ifelse(any.peak, 'blue','grey'))
axis(2,tck=0.02);axis(1,tck=0.02)
abline(a=0,b=1,col='grey80')
abline(a = -sp.th,b=1,col='blue',lty=2)
abline(v=c(-.5,.5),col='grey80'); abline(h=c(-.5,.5),col='grey80')



# fig2S - plot Time courses -----------------------
ev <- new.env()
load(file = '../data/lrpkmD_clean_data_all.Rdata',ev)
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F,row.names = 2)

pd <- ev$lrpkm_all[,grepl('[.]1$',colnames(ev$lrpkm_all))]
pd.scale <- cbind(t(scale(t(pd[,1:21]))), # z scale
                  t(scale(t(pd[,22:42]))))

pd.scale.1 <- cbind(apply(pd[,1:21],2,function(x) x/rowMax(pd[,1:21])), # max frac scale 
                    apply(pd[,22:42],2,function(x) x/rowMax(pd[,22:42]))) # max frac scale 


plotKmeansk(df=pd.scale.1,seed = 40)
set.seed(40)
kmeans.ord <- kmeans(pd.scale.1,9,iter.max = 2000,nstart = 25)$cluster


all.equal(rownames(pd.scale),rownames(k_ord))
pd.scale <- pd.scale[,grepl('[.]l|t',colnames(pd))]
head(cbind(pd.scale),rownames(pd.scale))
pd.scale <- melt(data.frame(pd.scale,gene=rownames(pd.scale)),varnames = gene)
pd.scale<- pd.scale[order(pd.scale$gene),]
n.gene <- length(unique(pd.scale$gene))
pd.scale$time <- rep(c(0,.5,1,2,3,5,8),4*n.gene)
n.time <- length(unique(pd.scale$time))
pd.scale$stimuli <- rep(rep(c("LPS","TNF"),each=n.time),2*n.gene)
pd.scale$genotype <- rep(rep(c('ctrl.','mt'),each=n.time*2),n.gene)
pd.scale$clust <- k_ord[pd.scale$gene,'cluster']
pd.scale$clust <- kmeans.ord[pd.scale$gene]

p <- ggplot(pd.scale,aes(time,value,group=gene,colour=factor(clust))) + 
  geom_line(position = position_jitter(width = 0.25, height = 2.5))

p + facet_grid(genotype ~ stimuli)
pd.scale$genotype <- factor(pd.scale$genotype)
p <- ggplot(pd.scale,aes(time,value,colour=stimuli,group=gene)) + 
  geom_line(position = position_jitter(width = 0.5, height = .15),alpha=.6)
p + facet_wrap(~ clust)

require(lattice)

cols <- colorRampPalette(c('lightgreen','orange','darkblue'))(9)
tmp.pd <- subset(pd.scale,clust %in% c(1,3))
tmp.pd <- pd.scale
col <-by(tmp.pd,tmp.pd$gene,
         FUN=function(x){
           v <- cols[x$clust]
           v[1]  ## I return one parameter , since I need one color
         }
)

xyplot(value ~ time| genotype + stimuli,tmp.pd,groups = gene,type='l',alpha=.5,
       par.settings=
         list(superpose.line   = list(col = col),                 ## color of lines
              superpose.symbol = list(col=col),                   ## colors of points
              add.text = list(col=c('darkgreen','orange'))))

# Fig.2s- plot-tc-by-cluster  -----------------------------------------------------
require(dplyr)
require(tidyr)
pd <- read.csv(file='./data/maxscale.csv',stringsAsFactors = F,row.names = 1)

#ggplot(pd.here,aes(time,rpkm,colour=stimuli))+ geom_line(aes(linetype=genotype)) + facet_wrap(batch~stimuli)
pd.tc_stat <- getAvgData()

for(b in c('b1',"b2")){
  ggplot(pd.tc_stat%>% filter(batch==b), aes(time,mean_val,colour=stimuli)) + geom_line(aes(linetype=genotype)) + geom_point(aes(shape=genotype))+ 
    geom_errorbar(aes(ymax=mean_val+std_value,ymin=mean_val-std_value,linetype=genotype)) + facet_grid(cluster~stimuli)
  ggsave(file=paste0('fig.S1_',b,'.b2.avg.by.cluster.pdf'),width = 6,height = 9)  
}

pd.tc_stat <- getAvgData(do.mergebatch = T)
ggplot(pd.tc_stat, aes(time,mean_val,colour=stimuli)) + geom_line(aes(linetype=genotype)) + geom_point(aes(shape=genotype))+ 
  geom_errorbar(aes(ymax=mean_val+std_value,ymin=mean_val-std_value,linetype=genotype)) + facet_grid(cluster~stimuli)
ggsave(file=paste0('fig.S1_b1.b2.avg.by.cluster.pdf'),width = 6,height = 9)  

# fig2S.test_induction  ---------------------------------------------------
range(pd)

pd.tbl <- pd %>% as.data.frame %>% mutate(gene=rownames(pd)) %>%
  gather(key=sample,value=lfc,-gene)
coldata <-pd.g[,3:6]

tmp.idx <- seq(1,85,by=7) 
for(i in 1:(length(tmp.idx)-1))
  pd.log2.fold[,tmp.idx[i]:(tmp.idx[i+1]-1)] <- pd.log2.fold[,tmp.idx[i]:(tmp.idx[i+1]-1)] -pd.log2.fold[,tmp.idx[i]]
range(pd.log2.fold)


# mRNA:plot heatmap by group  --------------------------------------------------
pd <- read.csv(file='../data/mRNA.nfkbgene.rpkm.csv',header = T,row.names = 1)
pd <- pd[,1:42]
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene)
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

# max scale 
pd.scale <- cbind(t(apply(pd[,1:21],1,function(x) x/max(x))),
                  t(apply(pd[,22:42],1,function(x) x/max(x))))


# category 
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../data/mRNA.sp.avg.csv",header = T,
                    row.names = 1,stringsAsFactors = F)

all.equal(rownames(sp.mRNA),rownames(pd.scale))
cateogry = ifelse(sp.mRNA$b1.LT.ctrl<0.5,'non-sp.',ifelse(sp.mRNA$b1.LT.ko-sp.mRNA$b1.LT.ctrl<=-0.5,
                                                                "dy","LPS-sp-no-dy")) #no-specific, LPS-sp but not dy, LPS-sp &dy 
rev(table(cateogry))
names(cateogry )<- rownames(sp.mRNA)

pd.scale.new<- rbind(pd.scale[cateogry=='non-sp.',],
                 pd.scale[cateogry=='LPS-sp-no-dy',],
                 pd.scale[cateogry=='dy',])
pd.scale.new <- pd.scale.new[,-c(15:21,36:42)]
pheatmap(pd.scale.new,cluster_cols = F,cluster_rows = F,scale = 'none',show_rownames = F,
         gaps_col =  seq(7,28,by=7),gaps_row = cumsum(table(cateogry)),
         show_colnames = F)
cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)
bks <- seq(0,1.0001,by = .1)
seps <- cumsum(rev(table(cateogry)))+1 # order: non-sp, lps-sp-no-dy,dy
sapply(1:2, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"fig2s_mRNA_by_cat",i,".eps"),onefile = F,width = 2,height = 6)
         dat <- pd.scale.new[,((i-1)*14+1):(i*14)]
         pheatmap.2(pd=dat,bks=bks,cols = cols,cwid = 5,
                    legend=F,gaps_col = grep('8',colnames(dat)))
         dev.off()
       })

### sp plot 
all.equal(rownames(pd.scale),rownames(sp.mRNA))
all.equal(rownames(pd.scale),names(cateogry))
sp.mRNA.new<- rbind(sp.mRNA[cateogry=='non-sp.',],
                    sp.mRNA[cateogry=='LPS-sp-no-dy',],
                    sp.mRNA[cateogry=='dy',])


cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)
sp.mRNA.new[sp.mRNA.new>3]<- 3 ;sp.mRNA.new[sp.mRNA.new< -3]<- -3 
sapply(1:2,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"fig2s_mRNA_sp_by_cat",i,".eps"),
                    onefile = F,width = 1.5,height = 6)
         pheatmap.2(pd=sp.mRNA.new[,i*3-2],bks=bks,cols=cols,legend=F,
                    cwid = 10)#,gaps_col=c(1,2))  
         dev.off()
       })

signif(bks,2)
# fig2E donut plot ------------------------------------------------------------------
sp.th <- 2^0.5

pd.venn.ko <- getVennData(gtype = 'ko',simplify = T)
pd.venn.ko.th <- getVennData(sp.th.val = 2,gtype = 'ko',simplify = T)

p.1 <- donutPlot(pd.venn.ko,inside.v = "LT.nosp")
p.2 <- donutPlot(pd.venn.ko.th,inside.v = "LT.nosp")

pd <- rbind(p.1$data,p.2$data)
pd$th <- c(rep(0.5,nrow(p.1$data)),rep(1,nrow(p.2$data)))
pd <- donutPlot(pd,trans = F,group.v = 'th') 
pd <- pd + theme(strip.text.x=element_blank())

ggsave(filename = paste0(fig_folder,'subfig2C.eps'),pd,width = 4,height = 2)


# caRNA:heatmap by group  -------------------------------------------------
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../data/mRNA.sp.avg.csv",header = T,
                    row.names = 1,stringsAsFactors = F)

cateogry = ifelse(sp.mRNA$b1.LT.ctrl<0.5,'non-sp.',ifelse(sp.mRNA$b1.LT.ko-sp.mRNA$b1.LT.ctrl<=-0.5,
                                                          "dy","LPS-sp-no-dy")) #no-specific, LPS-sp but not dy, LPS-sp &dy 
rev(table(cateogry))
names(cateogry )<- rownames(sp.mRNA)
seps <- cumsum(rev(table(cateogry)))+1 # order: non-sp, lps-sp-no-dy,dy

cpm.avg.gene <- read.csv(file='../data/caRNA.cpm.max.geno.scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)
sum(rownames(cpm.avg.gene) %in% k_ord$gene.2)
rownames(cpm.avg.gene) <- k_ord$gene
all.equal(rownames(cpm.avg.gene),rownames(sp.mRNA))

pd <- cpm.avg.gene
pd.new<- rbind(pd[cateogry=='non-sp.',],
                     pd[cateogry=='LPS-sp-no-dy',],
                     pd[cateogry=='dy',])


bks <- seq(0,1,by = 0.1)
cols <- colorRampPalette(brewer.pal(9,"Blues"))(10)

require(pheatmap)
pheatmap(pd.new,cluster_cols = F,cluster_rows = F,scale = 'none',show_rownames = F,color = cols,
         breaks = bks)

sapply(1:2, 
       function(i){
         #i <-2
         setEPS()
         postscript(paste0(fig_folder,"caRNA_max_bycat",i,".eps"),onefile = F,width = 2,height = 6)
         pheatmap.2(pd=pd.new[,((i-1)*12+1):(i*12)], bks = bks,gaps_col = c(6),cwid = 7,
                    legend = F,cols = cols)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"caRNA_max_scale.eps"))

# sp. caRNA 

sp.caRNA <- read.csv('../data/caRNA.sp.peak.csv',stringsAsFactors = F,row.names = 1)
all.equal(rownames(sp.caRNA),rownames(sp.mRNA))

sp.caRNA.log2.new <- rbind(sp.caRNA[cateogry=='non-sp.',],
                           sp.caRNA[cateogry=='LPS-sp-no-dy',],
                           sp.caRNA[cateogry=='dy',])


cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3.001)

sp.caRNA.log2.new[sp.caRNA.log2.new>3] <- 3
sapply(1:2,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"sp.avg.caRNA.cat",i,".eps"),
                    onefile = F,width = 1.5,height = 6)
         pheatmap.2(pd=sp.caRNA.log2.new[,i],bks=bks,cols=cols,legend=F,
                    cwid = 10)  
         dev.off()
       })

write.csv(file='./data/caRNA-max-avg-sp.csv',sp.caRNA.log2)



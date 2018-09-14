source(file='../auxilary_functions.R')
require("dplyr")
require("plyr")
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.3-4/subfigs/'
# residual plot ------------------------------------------------------
head(pd)

ph.pd<- (pd %>% spread(type,normCnt.frac) 
         %>% group_by(geno,sti,time,gene)
         %>% summarise(residues = Sim. - Exp. )
         %>% unite(key,geno,sti,time)
         %>% spread(key,residues))  
tmp.cols <- colorRampPalette(c('slateblue','white','sienna1'))(8); tmp.cols[c(4,5)] <- 'white'
pd.2 <- data.frame(ph.pd[,-1]); rownames(pd.2) <- ph.pd$gene
setEPS()
postscript(paste0(fig_folder,"v1_residuals.eps"),onefile = F,width = 6,height = 8)
pheatmap(pd.2[ord,c(15:28,1:14)],breaks = c(-0.7,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.9),
         scale = 'none',cluster_cols = F,cluster_rows = F,
         color = tmp.cols,gaps_col = seq(0,ncol(pd.2),by=7),fontsize_row = 8,
         annotation_row = anno.row,annotation_colors = anno.row.col) 
dev.off()

# all sim vs. exp ---------------------------------------------------------
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F)
sp.mRNA <- read.csv(file="../data/sp.avg.csv",header = T,
                    row.names = 1,stringsAsFactors = F)
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)

v1.simData <- read.csv(file='../fig.4_modelfit/mRNA-Fit-avg-sp-v1c/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
fit.r2 <- read.csv(file = '../fig.4_modelfit/mRNA-Fit-avg-sp-v1c/result.csv',row.names = 1,
                   stringsAsFactors = F)
v1.simData$r2 <- fit.r2[v1.simData$gene,'r2']
maxScale.mRNA$r2 <- fit.r2[maxScale.mRNA$gene,'r2']
require(dplyr)
v1.simData <- (v1.simData %>% arrange(desc(r2)))
maxScale.mRNA <- arrange(maxScale.mRNA,desc(r2))
genes <- unique(v1.simData$gene)

n.gene <- 69
tmp <- c(seq(1,n.gene,by=24),n.gene)

pdf(file='v1.cmp.pdf',width = 12)
for(i in 1:(length(tmp)-1)){
  sub.genes <- genes[tmp[i]:(tmp[i+1]-1)]
  pd.exp <- subset(maxScale.mRNA,gene %in% sub.genes) 
  pd.sim <- subset(v1.simData,gene %in% sub.genes)
  pd.exp$gene <- sapply(pd.exp$gene,function(g) paste(g,signif(fit.r2$r2[fit.r2$gene==g],2)))
  pd.sim$gene <- sapply(pd.sim$gene,function(g) paste(g,signif(fit.r2$r2[fit.r2$gene==g],2))) 
  pd.exp$gene <- factor(pd.exp$gene,levels = pd.exp$gene)
  pd.sim$gene <- factor(pd.sim$gene,levels = pd.sim$gene)
  
  p <- ggplot(pd.exp,aes(x=time,y=normCnt.frac,colour=sti)) + 
    geom_point(aes(shape=geno)) 
  p <- p+ geom_line(data = pd.sim,aes(linetype=geno)) + facet_wrap(~ gene, ncol = 8)
  print(p)
}
dev.off()
system("open v1.cmp.pdf")


# second version with specificity 
require(tidyr)
require(grid) 
rownames(sp.mRNA) <- sapply(rownames(sp.mRNA), function(x) k_ord$gene.2[k_ord$gene==x])

tmp.idx <- c(seq(1,69,by=6),70)
pdf(file='v1.cmp2.pdf',width = 8,height = 8)
sapply(1:(length(tmp.idx)-1),plotAllEgFit)
dev.off()
system('open v1.cmp2.pdf')

plotAllEgFit(1)

# peak vs. half-life  -----------------------------------------------------
pd.hf <- read.csv(file='../data/v4-hf-final.csv',row.names = 1,stringsAsFactors = F)
ev <- new.env()
load(file = '../data/lrpkmD_clean_data_all.Rdata',ev)
pd <- ev$lrpkm_all
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene,rownames(pd.hf))


pd.peak <- NULL 
for(i in 1:nrow(pd)){
  p <- plotExpSingleGene(i,maxFrac = F)  
    pd.peak <- rbind(pd.peak,p$data %>% group_by(gene,stimuli,genotype,batch) %>% summarise(peak_rpkm=max(value)))  
}

sum(rownames(pd.hf) %in% rownames(pd))
pd.peak <- pd.peak %>% left_join(data.frame(hf=round(10^pd.hf[,1]),
                                 gene=rownames(pd.hf),
                                 stringsAsFactors = F),by="gene") %>% arrange(hf)

ggplot(pd.peak,aes(hf,peak_rpkm,colour=stimuli)) + geom_point(aes(shape=genotype))# + geom_line(aes(linetype=genotype)) +
  coord_trans(x='log10',y='log2')
  
pd.peak <- pd.peak[complete.cases(pd.peak$hf),]
pd.peak$hf_bin <- cut(x=pd.peak$hf,c(0,60,480,max(pd.peak$hf)+1))

pd.peak.sum <- pd.peak %>% group_by(stimuli,genotype,batch,hf_bin) %>% summarise(mean_rpkm=mean(peak_rpkm),
                                                           std_rpkm = sd(peak_rpkm))

# violin plot 
p.v <- ggplot(filter(pd.peak,batch=='b1'),aes(hf_bin,peak_rpkm,colour=stimuli))
p.v <- p.v + geom_violin() + facet_wrap(~genotype)  + scale_y_log10() + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2,position = position_dodge(width = .9))
ggsave(file='fig3s_peak_vs_hf_violin.pdf',p.v,height = 4,width = 8)
  #geom_boxplot(width=0.1,aes(),position = position_dodge(width = .9))

# bar plot 
p <- ggplot(filter(pd.peak.sum,batch=='b1'),aes(hf_bin,mean_rpkm,colour=stimuli))
p <- p+ geom_bar(aes(fill=stimuli),stat = 'identity',position = "dodge")+ 
  facet_wrap(~genotype) 
p <- p + geom_errorbar(aes(ymin=mean_rpkm+std_rpkm,ymax=mean_rpkm+std_rpkm),
                  position = position_dodge(width = .9),width=.2) + 
  geom_linerange(aes(ymin=mean_rpkm,ymax=mean_rpkm+std_rpkm),position = position_dodge(width = .9))
ggsave(file='fig3s_peak_vs_hf.pdf',p)

# dot point / it will become fig 3 
m <- ggplot(pd.peak%>%filter(batch=="b1"),aes(hf,peak_rpkm,colour=stimuli))
m<-m + geom_point(aes(type=genotype),alpha=.5) + facet_wrap(~genotype)  + 
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')
m<- m + stat_smooth(se = T,aes(fill=stimuli))
ggsave(file='fig3s_peak_vs_hf_dot.pdf',m,height = 8,width = 8)

# sp vs. hf only control --------------------------------------------------

genes.hf <- read.csv(file="../data/v4-hf-final.csv",row.names = 1,stringsAsFactors = F)
genes.sp <- read.csv(file="../data/mRNA.sp.peak.csv",row.names = 1,stringsAsFactors = F)
genes.sp <- genes.sp[rownames(genes.hf),]
all.equal(rownames(genes.hf),rownames(genes.sp))

pd <- data.frame(hf_log10=genes.hf[,1],
                 lt_sp = genes.sp[,1],
                 lt_sp_mt = genes.sp$b1.LT.ko)
pd <- pd[complete.cases(pd),]
pd$hf <- 10^(pd$hf_log10)/60 
pd.cor <- cor.test(x = pd$hf,y = pd$lt_sp)

bks.hf <- c(.15,.5,1,3,5,14)
cols.hf <- colorRampPalette(c("azure","chartreuse4"))(length(bks.hf)-1)

setEPS()
postscript(file=paste0(fig_folder,'subfig3c.eps'),width = 2,height = 2)
par(mar=rep(0,4))
boxplot(pd$lt_sp~ cut(pd$hf,
                      bks.hf),
        #seq(min(pd$hf_log10),max(pd$hf_log10),length.out = 5)),
        outline = F,
        col= cols.hf,
        xlab = 'half-life(hrs)',
        ylab="TL specificity")
abline(h=.5,lty=2,col=grey(.5))
axis(2,tck=0.02);axis(1,tck=0.02)
grid()

# cluster based plot  -----------------------------------------------------

# version 2 
pd.hf <- read.csv(file='../data/v4-hf-final.csv',row.names = 1,stringsAsFactors = F)
ev <- new.env()
load(file = '../data/lrpkmD_clean_data_all.Rdata',ev)
pd <- ev$lrpkm_all
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene,rownames(pd.hf))

pd.peak <- NULL 
for(i in 1:nrow(pd)){
  p <- plotExpSingleGene(i,maxFrac = F)  
  pd.peak <- rbind(pd.peak,p$data %>% group_by(gene,stimuli,genotype,batch) %>% summarise(peak_rpkm=max(value),
                                                                                          avg_rpkm = mean(value)))  
}

pd.peak <- pd.peak %>% left_join(data.frame(hf=round(10^pd.hf[,1]),
                                            gene=rownames(pd.hf),
                                            stringsAsFactors = F),by="gene") %>% arrange(hf)
pd.peak <- pd.peak[complete.cases(pd.peak$hf),]
pd.peak <- pd.peak %>% gather(key=feature,value=rpkm,peak_rpkm,avg_rpkm)

# dot point / it will become fig 3 
m <- plotFeaturevsHf(b='b1',a = .5) + xlab('Half-life (mins)') +theme_bw()+ 
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'top',
        axis.title=element_blank())
ggsave(file=paste0(fig_folder,'fig3s_peak_vs_hf_b1.eps'),m,height = 4,width = 5,device = cairo_ps)


#install.packages("plot3D")
require(plot3D)
d.3d <- m$data %>% filter(feature=="peak_rpkm") %>% 
  mutate(lrpkm = log2(rpkm)) %>% select(-rpkm) %>%
  spread(key=stimuli,value=lrpkm) 
with(d.3d,scatter3D(tnf,il1,lps,colvar = NULL, col='blue',pch=16))

d.3d <- m$data %>% filter(genotype=='mt',feature=="peak_rpkm") %>% 
  mutate(lrpkm = log2(rpkm)) %>% select(-rpkm) %>%
  spread(key=stimuli,value=lrpkm) 
with(d.3d,scatter3D(tnf,il1,lps,colvar = NULL, col='blue',add = T))

require(plotly)
p <- plot_ly(d.3d, x= ~tnf,y= ~il1,z= ~lps,color= ~genotype,colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers()
p

# error model and final half-life -----------------------------------------
require(dplyr)
require(ggplot2)
pd.hf.all <- read.csv(file='../../half-life/Supriya/ActDBrowser/final.genome.kdeg.csv',stringsAsFactors = F)
#pd.hf.all <- subset(pd.hf.all,ensembleID %in% kb.genes$X)
pd.hf.all <- subset(pd.hf.all,sample %in% c('b1_0.0_U','b2_0.0_U'))
pd.hf.all <- pd.hf.all %>% mutate(hf=1/kdeg*60) %>% filter(hf>=0)
ggplot(pd.hf.all,aes(hf,adjR2)) + geom_point(alpha=.2) + scale_x_log10() +
  geom_vline(xintercept = as.numeric(quantile(pd.hf.all$hf, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)),
             linetype=2,colour="#3366ff")
require(LSD)
setEPS()
postscript(file= paste0(fig_folder,'whole_geome_adjR2_vs_hf.eps'))
heatscatter(x=log10(pd.hf.all$hf),y=pd.hf.all$adjR2,xlab = "hf(log10, mins)",ylab='adjR2',
            main = 'Whole genome scale half-life')
#lines(loess.smooth(x=log10(pd.hf.all$hf),y=pd.hf.all$adjR2),col=2,lwd=2)
abline(v = log10(as.numeric(quantile(pd.hf.all$hf, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))),
       lty=2,col="chartreuse4")
sapply(c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
       function(x) text(log10(as.numeric(quantile(pd.hf.all$hf,x))),
                        -0.4,paste0(x*100,'%, ',signif(as.numeric(quantile(pd.hf.all$hf,x)),2)),
                        pos = 4,srt=90,col = "chartreuse4"))
abline(h = (as.numeric(quantile(pd.hf.all$adjR2, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
sapply(c(0., 0.25, 0.5, 0.75, 0.99),
       function(x) text(4,(as.numeric(quantile(pd.hf.all$adjR2,x))),
                        paste0(x*100,'%, ',signif(as.numeric(quantile(pd.hf.all$adjR2,x)),2)),
                        pos=4,col = "#3366ff"))

dev.off()


# filter 0.86 adjR2 -------------------------------------------------------
require(tidyr)
pd.hf.all <- pd.hf.all %>% filter(adjR2>=0.86)
pd.hf.pair <- pd.hf.all%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
na.fil <- complete.cases(pd.hf.pair)
pd.hf.pair <- pd.hf.pair[na.fil,]
setEPS()
postscript(file= paste0(fig_folder,'whole_geome_adjR2_0.86.hf.eps'))
x <- log10(pd.hf.pair$b1_0.0_U[na.fil]); y <- log10(pd.hf.pair$b2_0.0_U[na.fil])
heatscatter(x,y,
            xlab = 'half-life (log10, mins,rep1)',
            ylab = 'half-life (log10, mins,rep2)',
            main = 'Half-life with adjR2>=.86')
lines(loess.smooth(x,y),col=2,lwd=2)
abline(b=1,a=log10(2),lty=2)
abline(b=1,a=-log10(2),lty=2)
#abline(b=1,a=0,lty=2)
#lm.fit <-lm(y~x)
#abline(lm.fit)
#abline(rlm(log10(pd.hf.pair$b2_0.0_U)~log10(pd.hf.pair$b1_0.0_U)),col=2)
#rlm(log10(pd.hf.pair$b2_0.0_U)~log10(pd.hf.pair$b1_0.0_U))
abline(h = (as.numeric(quantile(y, c(0.,.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
abline(v = (as.numeric(quantile(x, c(0.,.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
sapply(c(0.,.05, 0.25, 0.5, 0.75, 0.95,1.0),
       function(x1) text(1, as.numeric(quantile(y,x1,na.rm = T)),
                         paste0(x1*100,'%, ',signif(10^as.numeric(quantile(y,x1,na.rm = T)),2)),
                         pos=4,col = "#3366ff"))
sapply(c(0.,.05, 0.25, 0.5, 0.75, 0.95,1.0),
       function(x1) text(as.numeric(quantile(x,x1,na.rm = T)),1.2,
                         paste0(x1*100,'%, ',signif(10^as.numeric(quantile(x,x1,na.rm = T)),2)),
                         pos=4,srt=90,col = "#3366ff"))
#legend(x=2.2,y=3.2,legend = paste0("y=",
#                                   signif(as.numeric(lm.fit$coefficients[2]),3),'x+',
#                                   signif(as.numeric(lm.fit$coefficients[1]),2)),
#       lty = 1,col=1)

dev.off()


idx <- which(with(pd.hf.pair,b1_0.0_U > quantile(b1_0.0_U,.05,na.rm = T) & b1_0.0_U < quantile(b1_0.0_U,.95,na.rm = T)
                  &b2_0.0_U > quantile(b2_0.0_U,.05,na.rm = T) & b2_0.0_U < quantile(b2_0.0_U,.95,na.rm = T)))

pd.hf.pair$var <- apply(pd.hf.pair[,2:3],1,var)
pd.hf.pair$sd <- apply(pd.hf.pair[,2:3],1,sd)
pd.hf.pair$mean <- apply(pd.hf.pair[,2:3],1,mean)

ggplot(log10(pd.hf.pair[idx,c("mean","var")]),aes(mean,var)) + geom_point(alpha=.2) + 
  geom_smooth() #+ scale_x_log10() + scale_y_log10()

#abline(lm(log10(pd.hf.pair$sd[idx])~log10(pd.hf.pair$mean[idx])))
y <- log10(pd.hf.pair$sd[idx]); x<- log10(pd.hf.pair$mean[idx])
setEPS()
postscript(file= paste0(fig_folder,'whole_geome_adjR2_0.86_90pecent.errmodel.eps'))
heatscatter(x,y,xlab = 'mean half-life (log10,mins)',ylab='SD (log10,mins)',main='Error model (%5~95% genes)')
#lines(loess.smooth(x,y),col=2,lwd=2)
require(MASS)
rlm.fit <- rlm(y~x)
save(file='whole.genome.error.model.hf.RData',rlm.fit)
abline(rlm.fit,col=2)

legend(x=2.2,y=-2,legend = paste0("y=",
                                  signif(as.numeric(rlm.fit$coefficients[2]),2),'x',
                                  signif(as.numeric(rlm.fit$coefficients[1]),2)),
       lty = 1,col=2)
dev.off()


#matlines(new$x,pred.c.clim[,-1],lty = 1,col = 1)

#abline(lsfit(y=log10(pd.hf.pair$var),x=log10(pd.hf.pair$mean)),col=3)
#abline(lm(log10(pd.hf.pair$var[idx])~poly(log10(pd.hf.pair$mean[idx]),2)),col=2,lwd=2)

#lines(loess.smooth(y=log10(pd.hf.pair$var[idx]),x=log10(pd.hf.pair$mean[idx])),col=3,cex=3)
#identify(log10(pd.hf.pair$var)~log10(pd.hf.pair$mean))
#abline(lm(log10(pd.hf.pair$var[c(5865,7089)])~log10(pd.hf.pair$mean[c(5865,7089)])))
#abline(lm(log10(pd.hf.pair$var[c(5865,7723)])~log10(pd.hf.pair$mean[c(5865,7723)])))
#abline(a = 0,b = 1,col=2)

# Predicted final half-life  ----------------------------------------------
pd.hf.pair <- pd.hf.all%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
pd.hf.pair$mean <- apply(pd.hf.pair[,2:3],1,mean,na.rm=T)
new <- data.frame(x=log10(pd.hf.pair$mean))
sd.pred <- 10^as.numeric(predict.rlm(rlm.fit,newdata = new))
pd.hf.pair$lb <- pd.hf.pair$mean-sd.pred;pd.hf.pair$ub <- pd.hf.pair$mean+sd.pred
pd.hf.pair[,2:6] <- round(pd.hf.pair[,2:6])
pd.hf.pair$Symbol <- getSymbol(pd.hf.pair$ensembleID)
write.csv(file=paste0(fig_folder,'Table.S.genome.hf.csv'),
          pd.hf.pair,quote = F,row.names = F)

## Half-life length correlation 

anno.gene <- read.delim(file='../../../DurationDecoding/data/RNAseq/ActD/actd_b1-roberto.txt',
                      stringsAsFactors = F,sep = "\t")[,c(1,6)]
anno.gene$Geneid <- substr(anno.gene$Geneid,1,18)
rownames(anno.gene) <- anno.gene$Geneid
pd.hf.pair$Length <- anno.gene[pd.hf.pair$ensembleID,"Length"]

x<- log10(pd.hf.pair$Length);y <- log10(pd.hf.pair$mean)
cor <- cor.test(x,y)
setEPS()
postscript(file=paste0(fig_folder,'whole_geome_hf_vs_l.eps'))
heatscatter(x,y,xlab='Transcript Length (bp,log10)',
            ylab='Half-life (mins,log10)',main = "")
abline(h = (as.numeric(quantile(y, c(0.,.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
abline(v = (as.numeric(quantile(x, c(0.,.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
sapply(c(0.,.05, 0.25, 0.5, 0.75, 0.95,1.0),
       function(x1) text(2.5, as.numeric(quantile(y,x1,na.rm = T)),
                         paste0(x1*100,'%, ',signif(10^as.numeric(quantile(y,x1,na.rm = T)),2)),
                         pos=4,col = "#3366ff"))
sapply(c(0.,.05, 0.25, 0.5, 0.75, 0.95,1.0),
       function(x1) text(as.numeric(quantile(x,x1,na.rm = T)),1.0,
                         paste0(x1*100,'%, ',signif(10^as.numeric(quantile(x,x1,na.rm = T)),2)),
                         pos=4,srt=90,col = "#3366ff"))
text(4.5,3.0,paste0("Pearson's cor = ",signif(cor(x,y),2)),col=2)
dev.off()

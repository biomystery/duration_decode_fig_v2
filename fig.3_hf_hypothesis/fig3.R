source('../auxilary_functions.R')
subfig_dir <- '../figures/Fig.3/subfigs/'

# Fig3A: wt  ---------------------------------------------------------------

# do simulation for heatmap 
if(F){
  hfs<- 10^seq(0,3,length.out = 50)
  sim.data <- fun.runSim(hf = hfs,n = 6)
  save(file = 'half-life-sim-ctrl.Rdata',sim.data)
}
 
# heatmap
figs.ctrl <-plotHMsim(geno='ctrl',forFig = T)
setEPS()
postscript(file=paste0(subfig_dir,'fig3A_hm_ctrl.eps'),width = 4,height = 4)
with(figs.ctrl,multiplot(tf.p,p,layout = matrix(c(1,rep(2,9)),ncol = 1)))
dev.off()

# fig3A_lines.

hf<- c(60*8,60,5)
p.lines <- plotSimEg(hfs=hf)
ggsave(filename = paste0(subfig_dir,'fig3A_lines.eps'),plot = p.lines,width = 2,height = 4)

if(F){
  pd.2 <- get_fig3_pdata(hfs = hfs,is_ggplot = T,simdata = fun.runSim(hf=hfs,n=3))
  pd.3 <- get_fig3_pdata(hfs = hfs,is_ggplot = T,simdata = fun.runSim(hf=hfs,n=6))
  
  pd <- rbind(cbind(pd.1,n=1),
              cbind(pd.2,n=3),
              cbind(pd.3,n=6))
  pd$hf <- 10^pd$hf
  ggplot(data=pd,aes(x=time,y=value,colour=stimuli,group=n))+ geom_line(aes(linetype=factor(n))) + facet_grid(hf ~ stimuli)
}

# fig3A_dynamics-summary 
# peak fraction
# avg/auc fraction 

sum.sim1 <- figs.ctrl$p$data %>% 
  group_by(hf,stimuli) %>%
  summarise(peak=max(value),mean=mean(value)) %>% 
  group_by(hf) %>%
  summarise(ratio.peak=first(peak)/last(peak),
            ratio.mean=first(mean)/last(mean)) %>%
    gather(ratio_type,ratio,-hf)


bks<-c(0,.5,1,2,3)
cols <- colorRampPalette(col.sp[5:8])(5)
bks.2 <- c('(0,0.5]','(0.5,1]','(1,1.5]','(1.5,2]','(2,3]')
fig3a_dySp <-ggplot(sum.sim1,aes(x=ratio_type,y=hf,group=ratio_type)) +
  geom_tile(aes(fill=cut(log2(ratio),bks))) +
  scale_fill_manual(values = cols,breaks =bks.2) +
  theme_bw() + theme(strip.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.title.y=element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     legend.position="none",
                     plot.margin=margin()) +
  geom_hline(yintercept = c(log10(30),2),linetype=2,color='grey60')
ggsave(filename = paste0(subfig_dir,'fig3a_ctrl_peak.eps'),
       plot = fig3a_dySp %+% subset(sum.sim1,ratio_type=='ratio.peak'),width = 1,height = 4,scale = 2.5/4)

# Fig3B:mt hm ---------------------------------------------------------

figs <-plotHMsim(geno='mt',forFig = T)
setEPS()
postscript(file=paste0(subfig_dir,'fig3A_hm_mt.eps'),width = 4,height = 4)
with(figs,multiplot(tf.p,p,layout = matrix(c(1,rep(2,9)),ncol = 1)))
dev.off()


hf<- c(60*8,60,5)
p.lines <- plotSimEg(hfs=hf,genotype = 'mt')
ggsave(filename = paste0(subfig_dir,'fig3A_lines_mt.eps'),plot = p.lines,width = 2,height = 4)

## mt.
sum.sim2 <- figs$p$data %>% 
  group_by(hf,stimuli) %>%
  summarise(peak=max(value),mean=mean(value)) %>% 
  group_by(hf) %>%
  summarise(ratio.peak=first(peak)/last(peak),
            ratio.mean=first(mean)/last(mean)) %>%
  gather(ratio_type,ratio,-hf)


#last_plot()+ scale_fill_gradientn(colours=colorRampPalette(col.sp[5:8])(11))
fig3a_dySp <-ggplot(sum.sim2,aes(x=ratio_type,y=hf,group=ratio_type))  + geom_tile(aes(fill=cut(log2(ratio),bks))) +
 scale_fill_manual(values = cols,breaks=bks.2)+ 
  theme_bw() + theme(strip.text.x=element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.y=element_blank(),
                                axis.title.y=element_blank(),
                                axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                legend.position="none",
                                plot.margin=margin()) + 
  geom_hline(yintercept = c(log10(30),2),linetype=2,color='grey60')
ggsave(filename = paste0(subfig_dir,'fig3a_mt_peak.eps'),
       plot = fig3a_dySp %+% subset(sum.sim2,ratio_type=='ratio.peak'),width = 1,height = 4,scale = 2.5/4)

# Fig3C: sim cmp ctrl and mt ---------------------------------------------------------
pd <- data.frame(rbind(sum.sim1,sum.sim2),
                 geno=rep(c('ctrl.','mt.'),each=nrow(sum.sim1)))
pd$Sp <- log2(pd$ratio)

p<- ggplot(data.frame(hf=c(-0.5,-0.5,log10(30),log10(30),
                                log10(30),log10(30),2,2,
                                2,2,3.5,3.5),
                      Sp=rep(c(-0.5,3,3,-0.5),3),
                      gp=factor(rep(c(1,2,3),each=4))), # alpha(c("green"), c(.05,.1,.2))
           aes(hf,Sp))  +  geom_polygon(aes(fill=gp)) + scale_fill_manual(values =c('#F2FFF2','#E5FEE5','#CCFECC'))
  
p<- p+geom_line(data=pd%>% filter(ratio_type=="ratio.peak"),
                aes(linetype=geno))+ xlab("log10(mins)")+
  geom_vline(xintercept = c(log10(30),2),linetype=2,color='grey60')+
  geom_hline(yintercept = 0.5,linetype=2,color='grey60')
p<- p  + theme_bw()+ theme(axis.title = element_blank(),
                       axis.text = element_blank(),
                       legend.position = 'none') + 
  coord_cartesian(xlim=range(pd$hf),ylim=range(pd$Sp),expand = T) + 
  geom_vline(xintercept = seq(0,3,by=.5),colour=grey(.9),size=.2)+
  geom_hline(yintercept = seq(0,2,by=.5),colour=grey(.9),size=.2)

ggsave(filename = paste0(subfig_dir,'fig3_simCmp.eps'),height = 2,width = 2,plot = p)

# fig3C-half-life ---------------------------------------------------------

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

# include mts _ ggplot 
pd.gg <- pd %>% mutate(halflife=cut(hf,bks.hf)) %>% gather(genotype,sp,c(2,3))
ggplot(pd.gg,aes(halflife,sp)) + geom_boxplot(aes(fill=genotype),outlier.shape=NA) + scale_fill_manual(values=c('gray','white'))+
   geom_hline(yintercept = 0.5,linetype=2,colour=grey(0.5))+theme_bw() + 
  scale_y_continuous(limits = quantile(pd.gg$sp, c(0.1, 0.9)))+
  theme(axis.title = element_blank(),legend.position = 'none',axis.text = element_blank())

ggsave(filename = paste0(subfig_dir,'subfig3c_2.eps'),width = 3,height = 3)


# fig3D- ActD-seq heatmaps  -----------------------------------------------
# read the half-life data 
pd.hf.raw <- read.csv(file='../../half-life/Supriya/ActDBrowser/allNormCount.csv',
                      stringsAsFactors = F)
kb.genes <- read.csv(file='../data/mRNA.cluster.csv',stringsAsFactors = F)
pd.hf.raw <- subset(pd.hf.raw, X %in% kb.genes$X)
rownames(pd.hf.raw) <- pd.hf.raw$X; pd.hf.raw$X <- NULL 
pd.hf.raw <- pd.hf.raw[,grep('_U_',colnames(pd.hf.raw))]
pd.hf.proc <- log2(pd.hf.raw+1)
pd.hf.proc[,1:5] <- pd.hf.proc[,1:5]-pd.hf.proc[,1]
pd.hf.proc[,6:11] <- pd.hf.proc[,6:11]-pd.hf.proc[,6]
pd.hf.proc[pd.hf.proc>0] <- 0
ord <- order(apply(pd.hf.proc[,1:5],1,mean))
names(ord) <- rownames(pd.hf.proc)[ord]
pheatmap(pd.hf.proc[ord,],scale = 'none',
         cluster_cols = F,show_rownames = F,
         cluster_rows = F,gaps_col = 5)
pd.hf.pair.anno <- pd.hf.all%>% mutate(hf=1/kdeg)%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
mis <- setdiff(rownames(pd.hf.proc),pd.hf.pair.anno$ensembleID)
pd.hf.pair.anno <- rbind(pd.hf.pair.anno,
                         data.frame(ensembleID=mis,
                                    b1_0.0_U = NA,
                                    b2_0.0_U = NA))
rownames(pd.hf.pair.anno) <- pd.hf.pair.anno$ensembleID;
pd.hf.pair.anno$ensembleID <-NULL
pd.hf.pair.anno <- pd.hf.pair.anno[rownames(pd.hf.proc),]

pd.hf.bks <- c(7.5,15,30,60,120,240,480, 960) # categories 
pd.hf.cut <- apply(pd.hf.pair.anno*60,2,function(x) as.numeric(cut(x,breaks = pd.hf.bks)))
rownames(pd.hf.cut) <- rownames(pd.hf.pair.anno)
cols <- colorRampPalette(brewer.pal(n = 7, name =
                                          "RdYlGn"))(length(pd.hf.bks)*2-1)
cols<- cols[8:15]

for(i in 1:2){
  setEPS(); postscript(file=paste0(subfig_dir,'hf-val-rp',i,'.eps'),height = 14,width = 3)  
  pheatmap(pd.hf.cut[ord,i],scale = 'none',breaks= 1:(length(pd.hf.bks)-1),
           cluster_cols = F,show_rownames = F,color = cols,border_color = "grey60",
           cluster_rows = F,cellwidth = 12,show_colnames = F,legend = F)
  dev.off()  
}
plotLegend(cols,bks = pd.hf.bks,fnames = paste0(subfig_dir,'hf-heatmap-lg.eps'))

bks <- seq(0,-7,by=-0.5)
cols <- colorRampPalette(brewer.pal(n = 11, name ="BrBG")[6:11])(length(bks)-1)

setEPS(); postscript(file=paste0(subfig_dir,'hf-hm-rp',1,'.eps'))  
pheatmap(pd.hf.proc[ord,1:5],scale = 'none',
         cluster_cols = F,show_rownames = F,color = cols,
         cluster_rows = F,cellwidth = 24,show_colnames = F,legend = F)
dev.off()  

setEPS(); postscript(file=paste0(subfig_dir,'hf-hm-rp',2,'.eps'))  
pheatmap(pd.hf.proc[ord,6:11],scale = 'none',
         cluster_cols = F,show_rownames = F,color = cols,
         cluster_rows = F,cellwidth = 24,show_colnames = F,legend = F)
dev.off()  

plotLegend(cols,bks = bks,fnames = paste0(subfig_dir,'hf-heatmap-cnt-lg.eps'))
# fig3D-ActD example ------------------------------------------------------
source("../../half-life/Supriya/ActDBrowser/auxfunctions.R")
pd.hf.raw <- read.csv(file='../../half-life/Supriya/ActDBrowser/allNormCount.csv',
                      stringsAsFactors = F,row.names = 1)
coldata <- read.csv(file='../../half-life/Supriya/ActDBrowser/coldata.csv',stringsAsFactors = F,row.names = 1)
pd.hf.raw <- log2(pd.hf.raw+1)
res.hf.fit <- read.csv(file = "../../half-life/Supriya/ActDBrowser/final.genome.kdeg.csv",
                       stringsAsFactors = F)
gene.dic <- read.csv(file='../../half-life/Supriya/ActDBrowser/gene.dic.csv',
                     stringsAsFactors = F, row.names = 1)
# heatmap order ord
idx <- 1; 
g.target <- "Fos"; g<- rownames(gene.dic)[which(gene.dic$x==g.target)]
idx <- which(names(ord)==g); 
plotActD_eg <- function(idx){
  g <- names(ord)[idx]; g.target <- gene.dic[g,1]
  cat(paste("No.",idx,"gene is:",g.target))
  pd <-getGdata(g.target,fit = res.hf.fit,cnt = pd.hf.raw)
  print(pd$g)
  print(pd$kdeg.fit%>% filter(sample %in% c("b1_0.0_U","b2_0.0_U")))
  pd.2 <- plotActD.gene.log2.gg.key(pd,"U")
  ggplot(pd.2$pd.data,aes(AcDTime,log2normCnt,group=key,colour=key))+
    geom_point() + 
    geom_point(data=pd.2$pd.fit,shape=1,size=3) + 
    stat_smooth(data=pd.2$pd.fit,method = 'lm',se = F) + ggtitle(paste(idx,g.target)) 
}
p.all <- lapply(round(seq(3,176,length.out = 4)),plotActD_eg)
p.all <- lapply(p.all,function(x) x + theme_bw()+
                    scale_color_brewer(palette = "Paired") +
                    ylim(0,8.2))
require(gridExtra)
pdf(file=paste0(subfig_dir,"./subfig3d_actD_eg.pdf"),width = 4)
grid.arrange(p.all[[1]],p.all[[2]],p.all[[3]],p.all[[4]],ncol=1)
dev.off()

p.all <- lapply(p.all,function(x) x + theme(text = element_blank(),
                                            legend.position = "none"))

setEPS()
postscript(paste0(subfig_dir,"./subfig3d_actD_eg.eps"),width = 2.5,height = 5)
grid.arrange(p.all[[1]],p.all[[2]],p.all[[3]],p.all[[4]],ncol=1)
dev.off()
# half-life scatter  ------------------------------------------------------
# load half-life results 
pd.hf.all <- read.csv(file='../../half-life/Supriya/ActDBrowser/final.genome.kdeg.csv',stringsAsFactors = F)
pd.hf.all <- subset(pd.hf.all,ensembleID %in% kb.genes$X)
pd.hf.all <- subset(pd.hf.all,sample %in% c('b1_0.0_U','b2_0.0_U'))
pd.hf.pair <- pd.hf.all%>% mutate(hf=1/kdeg)%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
pd.hf.pair <-pd.hf.pair[complete.cases(pd.hf.pair[,2:3]),2:3]
#cor.test(pd.hf.pair$b1_0.0_U,pd.hf.pair$b2_0.0_U)

cor.test(log10(pd.hf.pair$b1_0.0_U*60),log10(pd.hf.pair$b2_0.0_U*60))
#require(LSD)

#heatpairs(pd.hf.pair)
#require(GGally)
#ggpairs(as.data.frame(pd.hf.pair))
require(dplyr)

pd <- pd.hf.all%>% mutate(hf=1/kdeg*60)%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
pd <- pd[complete.cases(pd),]
p1 <- ggplot(pd,aes(b1_0.0_U,b2_0.0_U)) + geom_point(alpha=0.2,colour='red') + 
  geom_smooth(method = 'lm') + xlab('Half-life (mins, rep.1)')+
  ylab('Half-life (mins, rep.2)') + scale_x_log10() + scale_y_log10()  +
  annotate('text',x=20,y=10000,label='cor = 0.85',colour='#3366ff')

p2 <- ggplot(pd%>% select(b1_0.0_U),aes(b1_0.0_U))+ 
  geom_density() + geom_vline(data = pd%>% select(b1_0.0_U) %>% summarise(median=median(b1_0.0_U)),
                              aes(xintercept=median),linetype=2)+
  theme(axis.title.x = element_blank()) + scale_x_log10()#+ coord_trans(x="log10")

p3 <- ggplot(pd%>%select(b2_0.0_U),aes(b2_0.0_U))+ 
  geom_density()+  geom_vline(data = pd%>% select(b2_0.0_U) %>% summarise(median=median(b2_0.0_U)),
                              aes(xintercept=median),linetype=2)+
  coord_flip() + theme(axis.title.y = element_blank())+ scale_x_log10()

p4 <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

require(gridExtra)
g<-grid.arrange(p2, p4, p1, p3, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
ggsave(filename = paste0(subfig_dir,'hf-scatter-whole.eps'),width = 6,height = 6,device =cairo_ps,plot = g)

p2 <- p2 + theme(axis.title = element_blank(),axis.text = element_blank())
p3 <- p3 + theme(axis.title = element_blank(),axis.text = element_blank())
p1 <- p1 + theme(axis.title = element_blank(),axis.text = element_blank())

# plot hf-fit examples 
#gid <- kb.genes$X[kb.genes$gene=='Ccl5']
gid <- kb.genes$X[kb.genes$gene=='Fos']
#gid <- kb.genes$X[kb.genes$gene=='Il1a']
setEPS();postscript(file=paste0(subfig_dir,'fos-hf.eps'),width = 4,
                    height = 4)
x=c(0,.5,1,3,6);y=as.numeric(pd.hf.proc[gid,1:5])
plot(x, y,pch=16,col=brewer.pal(12,'Paired')[1],xlab='Time(hr)',
     ylab='Log2(normalised counts)')
abline(lm(y[1:3]~x[1:3]),
       col=brewer.pal(12,'Paired')[1])
x<-c(0,.5,1,2,4,6); y<- as.numeric(pd.hf.proc[gid,6:11])
points(x, y,pch=16,col=brewer.pal(12,'Paired')[2])
abline(lm(y[1:3]~x[1:3]),
       col=brewer.pal(12,'Paired')[2])
legend(3,0,pch = 16,col=brewer.pal(12,'Paired')[1:2],
       legend = c('rep1','rep2'))
dev.off()
pd.hf.all %>% filter(ensembleID==gid)





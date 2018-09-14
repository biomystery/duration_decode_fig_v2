# focus on v1 unfitted genes+ caRNA 
rm(list=ls())
source(file='../auxilary_functions.R')
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.5_caRNA/subfigs/'
pd.fig5 <- list()
# fig.5A - heatmap  -------------------------------------------------
ev <- new.env();load(file="./old/caRNAFit-avg-sp-v6c/caRNAFit-avg-sp-pre-v6.Rdata",ev)
cpm.peak.gene.2 <- read.csv(file="./old/caRNAFit-avg-sp-v6c/caRNA.cpm.max.geno.scale.csv",header = T,
                          row.names = 1,stringsAsFactors = F)
cpm.peak.gene <- read.csv(file='../data/caRNA.cpm.max.geno.scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
all(rownames(cpm.peak.gene) %in% rownames(cpm.peak.gene.2))
cpm.peak.gene <- cpm.peak.gene.2 # exon reads 
ev<- new.env();load(file="../fig.4_modelv1/data/clustering_cateIII_clean.Rdata",ev)
g.dic <- read.csv(file = '../data/mRNA.cluster_old.csv',stringsAsFactors = F,row.names = 1)
dic.clust<-unique.data.frame(ev$pd %>% dplyr::select(gene,clust))
dic.clust <- dic.clust %>% arrange(clust,gene)
rownames(dic.clust) <-dic.clust$gene
cpm.peak.gene <- cpm.peak.gene[dic.clust$gene,]

seps <- cumsum(sapply(c(1,2,4,5,7,8,9),
                      function(x) sum(dic.clust$clust==x)))
ord <- unlist(sapply(c(1,2,4,5,7,8,9),
                     function(x) which(dic.clust$clust==x)))
pd.ca <- cpm.peak.gene
bks <- seq(0,1,by = 0.2)
cols <-colorRampPalette(c("azure","blue"))(6)
seps.2 <- seps


# combine with cytoRNA data 
ev <- new.env();load('../fig.4_modelv1/data//clustering_cateIII_clean.Rdata',ev)
pd.cyto<- ev$pd %>%filter(type=='Exp.')%>%
  #ph.pd<- pd %>%filter(type=='Sim.')%>%
  unite(key,geno,type,sti,time) %>%
  spread(key,normCnt.frac)
rownames(pd.cyto)  <- pd.cyto$gene; pd.cyto$gene<-NULL 
cl_ord <- c(1,2,4,7,5,9,8)
ord <- unlist(sapply(cl_ord,function(x) which(pd.cyto$clust==x)))
ord.2 <- unlist(sapply(cl_ord,function(x) which(dic.clust$clust==x)))
cord <- c(15:28,1:14)+2
all.equal(rownames(pd.cyto)[ord],rownames(pd.ca)[ord.2])

pd.merge <- cbind(pd.ca[ord.2,],pd.cyto[ord,cord])



# Calc sp  ----------------------------------------------------------------
pd.merge.sp <- data.frame(cyto.ctrl.sp=log2(rowMax(pd.merge[,25:31])/rowMax(pd.merge[,32:38])),
                          cyto.mt.sp=log2(rowMax(pd.merge[,39:45])/rowMax(pd.merge[,46:52])),
                          ca.ctrl.sp = log2(rowMax(pd.merge[,1:6])/rowMax(pd.merge[,7:12])),
                          ca.mt.sp = log2(rowMax(pd.merge[,13:18])/rowMax(pd.merge[,19:24]))) 
pd.merge.sp <- pd.merge.sp %>%  mutate(gene=rownames(pd.merge.sp))

# Sp pieplot  -------------------------------------------------------------
if(T){
  setEPS()
  postscript(file =  paste0(fig_folder,"subfig5A_caRNA_ctrl.anno.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.ctrl.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels =table(pd.merge.sp$ca.ctrl.sp>=0.5))
  dev.off()
  
  setEPS()
  postscript(file =  paste0(fig_folder,"subfig5A_caRNA_ctrl.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.ctrl.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels ="")
  dev.off()
  
  
  setEPS()
  postscript(file =  paste0(fig_folder,"subfig5A_caRNA_mt.anno.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.mt.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels =table(pd.merge.sp$ca.mt.sp>=0.5))
  dev.off()
  
  setEPS()
  postscript(file =  paste0(fig_folder,"subfig5A_caRNA_mt.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.mt.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels ="")
  dev.off()
}


# heatmap - order ---------------------------------------------------------
pd.merge.sp <- pd.merge.sp%>%
  mutate(cat=ifelse(ca.ctrl.sp<0.5,"I",ifelse(ca.mt.sp>=0.5,'II','III'))) %>%
  arrange(cat)

pd.merge.sp <- pd.merge.sp%>% mutate(ord =apply(pd.merge[pd.merge.sp$gene,1:6],1,which.max))%>%
  arrange(cat,ord)
#pd.merge.sp <- pd.merge.sp%>% arrange(cat,ca.ctrl.sp,ca.mt.sp,cyto.ctrl.sp,cyto.mt.sp)

# examine
pheatmap(pd.merge[pd.merge.sp$gene,],scale = "none",
         breaks =  seq(0,1.0001,by = 0.1),
         cluster_rows = F,cluster_cols = F,gaps_row = cumsum(table(pd.merge.sp$cat)),
         gaps_col = c(6,12,18,24,31,38,45),show_rownames = T,legend = F,
         color = colorRampPalette(brewer.pal(9,"Blues"))(10))

# real plot 
bks <- seq(0,1.0001,by = 0.1)
seps <- cumsum(table(pd.merge.sp$cat))+1;

   
if(T){
  cols=colorRampPalette(brewer.pal(9,"Blues"))(10);
  setEPS()
  postscript(paste0(fig_folder,"subfig5c_caRNA_ctrl.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,1:12], bks = bks,gaps_col = c(6),cwid = 7,
             legend = F,cols = cols)#ca_contorl 
  dev.off()
  
  setEPS()
  postscript(paste0(fig_folder,"subfig5c_caRNA_mt.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,13:24], bks = bks,gaps_col = c(6),cwid = 7,
             legend =  F,cols = cols)#ca_contorl 
  dev.off()
  plotLegend(cols=cols,bks=bks,
             fnames = paste0(fig_folder,"subfig5c_caRNA_scale.eps"))
}

# cytoRNA 
if(T){
  cols=colorRampPalette(brewer.pal(9,"Reds"))(10);
  setEPS()
  postscript(paste0(fig_folder,"subfig5c_cytoRNA_ctrl.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,25:38], bks = bks,gaps_col = c(7),cwid = 7,
             legend = F,cols = cols) #ca_contorl 
  dev.off()
  
  setEPS()
  postscript(paste0(fig_folder,"subfig5c_cytoRNA_mt.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,39:52], bks = bks,gaps_col = c(7),cwid = 7,
             legend = F,cols = cols,rn = T,fontsize=6) #ca_contorl 
  dev.off()
  
  plotLegend(cols=cols,bks=bks,
             fnames = paste0(fig_folder,"subfig5c_cytoRNA_scale.eps"))
  
}

# heatmap 1b- sp bar -----------------------------------------------------------------
cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)
pd.merge.sp <- pd.merge.sp %>% 
  mutate(cyto.ctrl.sp = ifelse(cyto.ctrl.sp>3,3,ifelse(cyto.ctrl.sp< -3, -3 ,cyto.ctrl.sp)))%>%
  mutate(cyto.mt.sp = ifelse(cyto.mt.sp>3,3,ifelse(cyto.mt.sp< -3, -3 ,cyto.mt.sp)))%>%
  mutate(ca.ctrl.sp = ifelse(ca.ctrl.sp>3,3,ifelse(ca.ctrl.sp< -3, -3 ,ca.ctrl.sp)))%>%
  mutate(ca.mt.sp = ifelse(ca.mt.sp>3,3,ifelse(ca.mt.sp< -3, -3 ,ca.mt.sp)))

if(T){
  nms <- c("subfig5c_cytoRNA_ctrl_sp.eps","subfig5c_cytoRNA_mt_sp.eps","subfig5c_caRNA_ctrl_sp.eps","subfig5c_caRNA_mt_sp.eps")
  for(i in 1:4){
    setEPS()
    postscript(paste0(fig_folder,nms[i]),
               onefile = F,width = 1.5,height = 6)
    pheatmap.2(pd=pd.merge.sp[,i],bks=bks,cols=cols,legend=F,
               cwid = 10)  
    dev.off()
  }
  plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"subfig5c_sp_scale.eps"))
}

signif(bks,2)

# save order 
pd.fig5$row_names <- (pd.merge.sp$gene)
pd.fig5$row_seps <- seps

saveRDS(file = "pd.fig5.Rds",object = pd.fig5)
pd.fig5 <-readRDS(file = "pd.fig5.Rds")
#pd.merge <- pd.merge[,c(1:12,25:38)] #control only 
seps.2 <- seps
cols.2 <- colorRampPalette(brewer.pal(9,"Reds"))(10)





# old ---------------------------------------------------------------------







require(plotly); ggplotly(p)
require(scales)
p %+% (p$data%>% filter(!clust %in% 1:3)) + scale_color_manual(values = c('FALSE'=hue_pal()(2)[2])) 
p.2 <- p %+% (p$data%>% filter(clust %in% 1:3))
p.2 + theme_bw()+coord_equal(xlim=c(min(p.2$data$ca.sp),max(p.2$data$cyto.sp)),
                                    ylim=c(min(p.2$data$ca.sp),max(p.2$data$cyto.sp)))+
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank())+
  scale_color_manual(values = c('TRUE'=grey(.2))) 
ggsave(filename = paste0(fig_folder,'ctrl.sp.scatter.eps'),width = 3,height = 3)
pd.merge.sp <- pd.merge.sp %>% mutate(is.cyto.sp = cyto.sp>=0.5,is.ca.sp=ca.sp>=0.5)
pd.merge.sp <- pd.merge.sp %>% mutate(is.v1= !clust %in% 1:3)

ggplot(pd.merge.sp,aes(is.ca.sp,group=is.v1)) + geom_bar(aes(fill=is.v1),
                                                         position = "dodge")




# average by clustering  --------------------------------------------------
pd.sum <- pd.ca %>% mutate(gene= rownames(pd.ca))%>%
  mutate(clust=dic.clust[gene,"clust"])%>% 
  gather(key = samples,value = normCnt.frac,1:24) %>%
  separate(samples,c("geno","sti","time"),sep = "_")%>%
  mutate(key=paste(time,sti,geno,clust,sep="_")) %>%
  group_by(key)%>%
  summarise(m = mean(normCnt.frac),
            e = std_err(normCnt.frac))%>% 
  separate(key,into=c("time","sti","geno","clust"),sep = "_")
pd.sum<- pd.sum%>% mutate(time=as.numeric(time))

pd.sum$sti<- toupper(sub("T","tnf",sub("L","lps",pd.sum$sti)))
pd.sum$geno<- sub("D","mt",sub("I","wt",pd.sum$geno))

names(col.map) <- toupper(names(col.map))
pd.sum$geno <- factor(pd.sum$geno,levels = c("wt","mt"))
pd.sum$clust <- factor(pd.sum$clust,1:10)

ggplot(pd.sum,aes(time,m,colour=sti))  + 
  geom_line(aes(linetype=geno))+
  geom_point(aes(shape=geno))+
  facet_wrap(~clust,ncol = 1)+ 
  geom_errorbar(aes(ymin=m-e,ymax=m+e),width=.2)+ 
  scale_color_manual(values = col.map)+ ylab("caRNA")+
  theme_bw()+theme(strip.text.x = element_blank())


# cmp sp ------------------------------------------------------------------

pd.ca %>%mutate(gene=rownames(pd.ca))%>% 
  gather(key = sample,value = caRNA,1:24)%>%
  separate(col=sample,c('geno',"sti",'time'),sep = "_")%>%
  group_by(interaction(gene,geno,sti))%>% 
  summarise(m=max(caRNA))



sapply(1:2, 
       function(i){
         #i <-2
         setEPS()
         postscript(paste0(fig_folder,"caRNA_max",i,".eps"),onefile = F,width = 2,height = 6)
         pheatmap.2(pd=pd[,((i-1)*12+1):(i*12)], bks = bks,gaps_col = c(6),cwid = 7,
                    legend = F,cols = cols)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"caRNA_max_scale.eps"))

# fig.5A - sp comparason (peak) -------------------------------------------------

sp.caRNA <- read.csv(file='./data/caRNA-max-avg-sp.csv',header = T,
                     row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)
sum(rownames(sp.caRNA) %in% k_ord$gene)
sp.caRNA <- sp.caRNA[k_ord$gene,]
all.equal(rownames(sp.caRNA) , k_ord$gene)
all.equal(rownames(cpm.peak.gene) , k_ord$gene.2)

caRNA.sum<- data.frame(ctrl.L = rowSums(cpm.peak.gene[,1:6]),
           ctrl.T = rowSums(cpm.peak.gene[,7:12]),
           mt.L = rowSums(cpm.peak.gene[,13:18]),
           mt.T = rowSums(cpm.peak.gene[,19:24]))

sp.caRNA.2 <- data.frame(LT.ctrl=caRNA.sum$ctrl.L/caRNA.sum$ctrl.T,
                       LT.mt = caRNA.sum$mt.L/caRNA.sum$mt.T)
sp.caRNA.2[apply(sp.caRNA.2,2,is.infinite)] <- NA
sp.caRNA.2[sp.caRNA.2==0] <- NA
sp.caRNA.log2 <- log2(sp.caRNA.2)
rownames(sp.caRNA.log2) <- rownames(caRNA.sum)
all.equal(rownames(caRNA.sum),rownames(sp.mRNA))
data.frame(rownames(caRNA.sum),rownames(sp.mRNA))
#sp.caRNA.log2[sp.caRNA.log2>3] <- 3 ; sp.caRNA.log2[sp.caRNA.log2 < -3] <- -3 

cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)


sapply(1:2,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"dAUC_maxFrac_genoScale",i,".eps"),
                    onefile = F,width = 1.5,height = 6)
         pheatmap.2(pd=sp.caRNA.log2[,i],bks=bks,cols=cols,legend=F,
                    cwid = 10)  
         dev.off()
       })


plotLegend(cols = cols,bks = bks,fnames = paste0(fig_folder,"dAUC_maxFrac_genoScale_lg.eps"))
signif(bks,2)

## scatter 
all.equal(rownames(sp.mRNA),k_ord$gene)
sp.mRNA.log2 <- sp.mRNA
range(sp.caRNA,na.rm = T)
sp.caRNA.log2 <- sp.caRNA

plot(sp.mRNA.log2$b1.LT.ctrl,sp.caRNA.log2$LT.ctrl)
abline(h = c(-.5,.5),v=c(-.5,.5))
cor(sp.mRNA.log2$b1.LT.ctrl,sp.caRNA.log2$LT.ctrl,use = "na.or.complete")

g.dic <- k_ord$gene.2;names(g.dic)<- k_ord$gene
pd <- data.frame(mRNA.sp = c(sp.mRNA.log2$b1.LT.ctrl, sp.mRNA.log2$b1.LT.ko),
                 caRNA.sp =c( sp.caRNA.log2$LT.ctrl,sp.caRNA.log2$LT.mt))
pd[pd[,1:2]>3] <- 3 ; pd[pd< -3] <- -3 
pd$genotype= factor(rep(c('ctrl','mt'),each=nrow(sp.mRNA.log2)))
pd$gene = as.character(g.dic[rownames(sp.mRNA)])
pd$isDy = pd$gene %in% unique(fit.r2.caRNA$gene)
pd$isFit = pd$gene %in% unique(fit.r2$gene[fit.r2$r2>0.7])
cateogry = ifelse(sp.mRNA.log2$b1.LT.ctrl<=0.5,'non-sp.',ifelse(sp.mRNA.log2$b1.LT.ko-sp.mRNA.log2$b1.LT.ctrl<=-0.5,
                                                                "dy","LPS-sp-no-dy")) #no-specific, LPS-sp but not dy, LPS-sp &dy 
table(cateogry)
pd$clust= rep(k_ord$cluster,2)
pd$cat <- rep(cateogry,2)

pd.cors <- pd %>% group_by(genotype,clust) %>% summarise(cor=signif(cor(mRNA.sp,caRNA.sp,use = "na.or.complete"), 2))
pd.cors <- pd %>% group_by(genotype,cat) %>% summarise(cor=signif(epi.ccc(mRNA.sp,caRNA.sp)$rho.c$est, 2))


require(epiR)
epi.ccc()
pd$cat <- factor(pd$cat,levels = c("non-sp." ,    "LPS-sp-no-dy",  "dy"          ))


#pd.2%>%filter(cat=='dy'& type=='caRNA.sp') %>% summarise(n=sum(ctrl>=.5 &mt-ctrl<=-0.5,na.rm = T)/length(ctrl)*100)
#pd.2%>%filter(cat=='non-sp.'& type=='caRNA.sp') %>% summarise(n=sum(ctrl<.5,na.rm = T))

p <- ggplot(pd,aes(mRNA.sp,caRNA.sp))+geom_point() + facet_grid(cat~genotype)
pd.3 <- pd.2 %>% group_by(gene)%>%select(-c(ctrl,mt))%>%spread(type,delta)


p <- ggplot(pd.3, aes(mRNA.sp,caRNA.sp))+geom_point() + facet_grid(cat~.)
p + geom_abline(intercept = 0,slope=1,linetype=2,colour='cornflowerblue')+ theme_bw() + geom_vline(xintercept = -.5,linetype=2,colour='cornflowerblue')+
  geom_hline(yintercept = -.5,linetype=2,colour='cornflowerblue')+xlab('delta_sp_mRNA')+ylab('delta_sp_caRNA')
#p<- p + 
p   + geom_smooth(method = 'lm',se = F) + theme_bw() + geom_hline(yintercept = c(.5),linetype=2,colour='cornflowerblue')+
  geom_vline(xintercept = c(.5),linetype=2,colour='cornflowerblue') + geom_abline(slope = 1,intercept = 0,linetype=2,colour='cornflowerblue')+ 
  geom_text(data=pd.cors, aes(label=paste("C=", cor, sep=""), x=-1.5, y=2.5),colour='red',                  hjust=0)

pd.2 <- pd %>% gather(type,sp,1:2) %>% spread(genotype,sp) %>% mutate(delta=mt-ctrl)
pd.2%>%filter(type=='caRNA.sp') %>% summarise(n=sum(ctrl>=.5,na.rm = T)) #LPS-specific 
pd.2%>%filter(type=='caRNA.sp'&cat!='non-sp.') %>% summarise(n=sum(ctrl>=.5,na.rm = T)) #LPS-specific 
pd.2%>%filter(type=='mRNA.sp') %>% summarise(n=sum(ctrl>=.5,na.rm = T)) #
pd.2%>%filter(type=='caRNA.sp') %>% summarise(n=sum(ctrl>=.5 & mt-ctrl <= -.5,na.rm = T)) #dy-specific 
pd.2%>%filter(type=='caRNA.sp'&cat=='dy') %>% summarise(n=sum(ctrl>=.5 & mt-ctrl <= -.5,na.rm = T)) #LPS-specific 

p <- ggplot(pd.2, aes(ctrl,mt))+geom_point() + facet_grid(cat~type)
pd.cors <- pd.2 %>% group_by(genotype,cat) %>% summarise(cor=signif(epi.ccc(mRNA.sp,caRNA.sp)$rho.c$est, 2))
p   + theme_bw() + geom_hline(yintercept = c(.5),linetype=2,colour='cornflowerblue')+
  geom_vline(xintercept = c(.5),linetype=2,colour='cornflowerblue') + geom_abline(slope = 1,intercept = c(0,-.5),linetype=2,colour='cornflowerblue')
  
ggsave(file='fig3s_fitted_sp_sc.pdf')
ggsave(filename = paste0(fig_folder,'subfig5a_spscatter.eps'),width = 3,height = 5,p,scale = .75)


# raw-peak/avg vs. hf -----------------------------------------------------
pd.hf <- read.csv(file='../data/v4-hf-final.csv',row.names = 1,stringsAsFactors = F)
ev <- new.env()
pd <- read.csv(file='../data/caRNA-cpm-max-geno-scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F)
pd <- pd[k_ord$gene.2,]
all.equal(rownames(pd),k_ord$gene.2,rownames(pd.hf))

pd.peak <- NULL 
for(i in 1:nrow(pd)){
  p <- plotExpSingleGene.2(i,maxFrac = F)  
  pd.peak <- rbind(pd.peak,p$data %>% group_by(gene,stimuli,genotype) %>% summarise(peak_rpkm=max(value),
                                                                                          avg_rpkm = mean(value)))  
}

pd.peak <- pd.peak %>% left_join(data.frame(hf=round(10^pd.hf[,1]),
                                            gene=rownames(pd.hf),
                                            stringsAsFactors = F),by="gene") %>% arrange(hf)
pd.peak <- pd.peak[complete.cases(pd.peak$hf),]
pd.peak <- pd.peak %>% gather(key=feature,value=rpkm,peak_rpkm,avg_rpkm)

# dot point / it will become fig 3 
m <- plotFeaturevsHf.2(a = .5) + xlab('Half-life (mins)') +theme_bw()
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'top',
        axis.title=element_blank())
ggsave(file=paste0(fig_folder,'fig5s_peak_vs_hf.eps'),m,height = 4,width = 5,device = cairo_ps)


# ccl1-sim hf -------------------------------------------------------------

fit.r2.ca <- read.csv(file = '../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-funs-v6-r2.csv',row.names = 1)
round(log(2)/fit.r2.ca['Ccl1',]$k_deg)
pturbs <- c(1/6,1/5,1/4,1/3,1/2,1,2,3,4,5,6)
signif(94*c(1/6,1/5,1/4,1/3,1/2,1,2,3,4,5,6)/60,2)
pd.sp <- readRDS(file="../../notebooks/kdeg_v2.Rdata")
pd.sp.2 <- readRDS(file="../../notebooks/kdeg_v1.Rdata")
pd.sp$gene<-"Ccl1"
pd.sp.2$gene<- "Ccl2"
pd.sp.2$type <- pd.sp$type

pd.2 <- pd <- rbind(pd.sp,pd.sp.2)%>% filter(!type %in% c("Sim.norm","Exp."))
ggplot(pd, aes(type,LvT.sp,colour=gene))+
  geom_line(arrow=arrow(type = "closed",angle = 15,ends = "first"),
            position=position_dodge(width = 0.4))+
  geom_line(data=pd.2,aes(type,LvT.sp,group=interaction(geno,gene),
                          linetype=geno))+xlab("half-life fold-change")+
  theme_bw()+scale_x_discrete(labels=rev(c("1/6","1/5","1/4","1/3","1/2","1",
                                       "2","3","4","5","6")))+
  scale_color_brewer(palette = "Set1")

ggsave(file=paste0(fig_folder,'cc1-ccl2-sim.eps'),width = 5,height = 4)



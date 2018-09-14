source('../auxilary_functions.R')
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure//Fig.4/subfigs/'


# Fig4B  - model vs. data hm---------------------------------------------------------
# see preprocessing
load('../fig.4_modelv1/data/clustering_cateIII_clean.Rdata')
rord <- read.csv(file="../fig.2_knockout/gene_cat.csv",stringsAsFactors = F)
rord <- rord%>% filter(cate=="III")


ph.pd<- pd %>%filter(type=='Exp.')%>%
#ph.pd<- pd %>%filter(type=='Sim.')%>%
  unite(key,geno,type,sti,time) %>%
  spread(key,normCnt.frac)
rownames(ph.pd)  <- ph.pd$gene; ph.pd$gene<-NULL 
all(rownames(ph.pd) %in% rord$gene.2)
rord <- (rord %>% filter(gene.2 %in% rownames(ph.pd)))$gene.2


bks <- seq(0,1.0001,by = .1)
anno_row <- data.frame(v1=as.factor(ph.pd[,'v1']))
anno_row_color <- list(v1=c('TRUE'='yellow','FALSE'='black'))
rownames(anno_row) <-  rownames(ph.pd)
cord <- c(15:28,1:14)+2
ph.pd.sub <-  ph.pd[rord,cord]
cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)

#postscript(file = paste0(fig_folder,'sim_hm.eps'),
postscript(file = paste0(fig_folder,'exp_hm.eps'),
           onefile = F,width = 2.5,height = 5)
pheatmap(ph.pd.sub,scale = "none",cluster_rows = F,
         cluster_cols = F,show_rownames = F,
         color = cols,breaks = bks,
         fontsize_row = 5,cellwidth = 4,
         #annotation_row = genes.hf,
         annotation_colors = anno_row_color,annotation_row = anno_row,
         legend = F,annotation_names_row = F,
         show_colnames = F,annotation_legend = F,
         gaps_col = c(7,14,21))
dev.off()

pheatmap(ph.pd.sub,scale = "none",cluster_rows = F,
         cluster_cols = F,show_rownames = T,
         color = cols,breaks = bks,
         fontsize_row = 5,
         #annotation_row = genes.hf,
         annotation_colors = anno_row_color,annotation_row = anno_row,
         legend = F,annotation_names_row = F,
         show_colnames = F,annotation_legend = F,
         gaps_col = c(7,14,21))

if(T){
  fig4.setting <- list(ord= ord,
                       gord = rownames(ph.pd)[ord],
                       seps= seps,
                       bks = bks,
                       cols =cols)
  fig4.setting$hm_data <- ph.pd
  saveRDS(file='fig4.setting.rds',fig4.setting)  
}


# Fig4F -------------------------------------------------------------------
sp.th <- 2^0.5
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
pd <- data.frame(ctrl.sp= sp.mat$b1.LT.ctrl,
                 mt.sp = sp.mat$b1.LT.ko,
                 gene = rownames(sp.mat),
                 stringsAsFactors = T) 
pd$gene <- gene.dic[as.character(pd$gene),"gene.2"]
pd <- subset(pd,gene %in% rownames(anno_row))
selc <- rownames(anno_row)[anno_row$v1=="TRUE"]

p<- ggplot(data.frame(ctrl.sp=c(0.5,0.5,8,8,0.5,0.5,8,8,-3,-3,0.5,0.5),
                      mt.sp=c(0.5/sp.th,8,8,8/sp.th, 
                              -3,0.5/sp.th,8/sp.th,-3,
                              -3,8,8,-3),
                      gp=factor(rep(c(2,3,1),each=4))),
           aes(ctrl.sp,mt.sp))+   coord_cartesian(xlim=range(pd[,1:2]),ylim=range(pd[,1:2])) +
  geom_polygon(aes(fill=gp)) + scale_fill_manual(values = alpha(c("blue"), c(.05,.1,.2)))+
  theme_bw()+
  geom_abline(slope = 1,intercept = 0,colour=grey(.8)) +
  geom_abline(slope = 1/sp.th,intercept = 0,colour="blue",linetype=2) +
  geom_hline(yintercept = c(-.5,.5),colour=grey(.8),linetype=2)+
  geom_vline(xintercept = c(-.5,.5),colour=grey(.8),linetype=2) +
  geom_point(data = pd,aes(ctrl.sp,mt.sp),size=2,shape=21,fill="grey70",colour="black") + 
  annotate('point',x=pd[pd$gene%in%selc,1],y=pd[pd$gene%in%selc,2],fill='yellow',colour="black",size=2,shape=21)+ theme(legend.position = 'none')+
 coord_cartesian(xlim=c(0.5,5),ylim=c(-1,4))
ggsave(file=paste0(fig_folder,'subfig4f_w_ticks.eps'),plot = p,width = 2.5,height = 2.5,device=cairo_ps)
p +  theme(legend.position = 'none',axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())

ggsave(file=paste0(fig_folder,'subfig4f.eps'),width = 2.5,height = 2.5,device=cairo_ps)

# explore
pd <- pd %>% mutate(V1 = anno_row[gene,])
pd <- pd %>% mutate(dynamics_dependence = ctrl.sp-mt.sp)
p<- ggboxplot(pd,x="V1",y="dynamics_dependence",color="V1",palette = c("black","#E7B800"),
          ylab = "Dynamics dependence",add="jitter") 
p+stat_compare_means(method.args = list(alternative="less"),method="t.test",label.y = 4.5) +theme_bw()+ 
  theme(legend.position = "none",strip.background  = element_rect(fill="grey100" ))
ggsave(file=paste0(fig_folder,'subfig4f_explore.eps'),width = 2,height = 4,device=cairo_ps)

# explore 2
pd.2 <- pd%>% gather(key="Geno",value = "sp",1:2)

#  unite("Geno_Fit",c("Geno","isFit"),remove = F)
p.2<- ggboxplot(pd.2,x="V1",y="sp",color="V1",palette = c("black","#E7B800"),
                add="jitter",facet.by = "Geno")
p.2+theme_bw()+stat_compare_means(method.args = list(alternative="greater"),method = "t.test",label.y = 5)+
  theme(legend.position = "none",strip.background  = element_rect(fill="grey96" ))
ggsave(file=paste0(fig_folder,'subfig4f_explore2.eps'),width = 4.2,height = 4,device=cairo_ps)

# Fig4C examples- model vs. data -------------------------------------------------------
eg.genes <- c('Rel',"Nfkb1","Ccl1","Ccl2","Gsap",'Rab15','Mmp3')
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 

v1.simData <- read.csv(file='./models/mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
v1.simData$gene <- gene.dic[v1.simData$gene,'gene.2']

tmp.simData <- read.csv(file = "./models/mRNA-Fit-avg-sp-v1e/bestFit.tc.csv",stringsAsFactors = F)
colnames(tmp.simData) <- sub('mRNA','normCnt.frac',colnames(tmp.simData)) 
v1.simData <- (rbind(v1.simData,tmp.simData))

# 

# 
col.map <- c(lps="#f8766d",tnf="#00ba38",il1="#619cff")
p <- plotEgFit(eg.genes = eg.genes,wide = F)
p<- p + theme_bw() + theme(legend.position = "none",strip.text.x = element_blank())
p + theme(text = element_blank())
ggsave(filename = paste0(fig_folder,'v1_eg_anno.eps'),width = 2,height = 6)
ggsave(filename = paste0(fig_folder,'v1_eg.eps'),width = 2,height = 6)

# Fig4D plot perturb -------------------------------------------------------------
# subset the genes passed fitting 
pd.tmp <- ph.pd[ord,]; 
target.genes <- rownames(pd.tmp)[pd.tmp$v1]; rm(pd.tmp)

# loading simulatino data 
envs <- new.env()
pd.p <- readRDS('~/Dropbox/Projects/DurationDecoding-code/notebooks/kdeg_15m_new.Rdata')

pd.p <- (pd.p %>% filter(type!='Exp.')) #exclude exp. data 
pd.p$type <- sub("k_deg","k_deg_up",pd.p$type) # rename 
pd.p$gene <- as.character(pd.p$gene) 

# change name from ENSEM.. to gene symbol 
all(unique(pd.p$gene)%in% rownames(gene.dic))
pd.p$gene <- gene.dic[pd.p$gene,"gene.2"]
all( target.genes%in% pd.p$gene)

# plot
if(T){
  cols <- c(brewer.pal(3,"Set1")[1:2],'black')
  eps = 0.2
  pd.p<- pd.p %>% filter(gene %in% target.genes)
  pd.p$gene <- factor(pd.p$gene,levels = target.genes)
  
  pd.p <- pd.p %>%arrange(gene,type,geno)
  ggplot(pd.p,
         aes(ifelse(type=="Sim.norm",as.numeric(gene)+eps,
                    as.numeric(gene)-eps),
             LvT.sp,group=interaction(gene,type))) +
    #coord_flip(xlim = c(0,length(unique(pd$gene))+1),expand = F)+
    geom_point(aes(shape=geno,colour=type))+
    geom_line(aes(colour=type))+
    #arrow=arrow(type = "closed",angle = 15,ends = "last",
    #          length = unit(.1,"inches")))  +
    scale_color_manual(name="half-life",values = cols[-1],labels=c("15m","norm."))+
    scale_shape_manual(values = c(16,1))+
    scale_x_continuous(breaks = 1:length(unique(pd.p$gene)),
                       labels = unique(pd.p$gene))+ 
    theme_bw()+theme(panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_line(size = .5),
                     axis.text.x =element_text(angle = 90,hjust = 1,vjust = 0.5),
                     axis.title.x = element_blank())+
    geom_hline(yintercept = .5,linetype=2,colour="blue")+
    coord_cartesian(xlim = c(0,length(unique(pd.p$gene))+1),expand = F,
                    ylim = c(-0.2,2.6)) 
  
}
ggsave(filename = paste0(fig_folder,"feature_plot_down.pdf"),width = 4.43,height = 1.63,scale = 1.5)
ggsave(filename = paste0(fig_folder,"subfig4e_perturb.eps"),width = 7,height = 4)

# test Fig4D --------------------------------------------------------------
test <- readRDS('~/Dropbox/Projects/DurationDecoding-code/notebooks/kdeg_15m_tc_test.Rdata')

pd.tc <- test %>% 
  mutate(time=ifelse(type=="Sim.norm",time,time/60))
pd.tc$type <- factor(pd.tc$type,levels = c("Sim.norm","k_deg"))
pd.tc.2 <- pd.tc %>% 
  group_by(gene,geno,type) %>% 
  dplyr::mutate(normCnt.frac.scale_geno=normCnt.frac/max(normCnt.frac))%>%
  ungroup()%>%
  group_by(gene,type)%>%
  dplyr::mutate(normCnt.frac.scale=normCnt.frac/max(normCnt.frac))

ggplot(pd.tc.2,aes(time,normCnt.frac.scale,colour=sti,linetype=geno)) + geom_line() + 
  facet_grid(gene~type)+ scale_color_manual(values = col.map) + theme_bw()+ 
  theme(legend.position = "none")

pd.tc.2 %>%
  group_by(gene,geno,type,sti)%>%
  dplyr::summarise(m = max(normCnt.frac))%>%
  ungroup()%>%
  group_by(gene,geno,type)%>%
  dplyr::summarise(log2(m[1]/m[2]))%>%
  arrange(type)



# Fig4E perturb plot with lines  ------------------------------------------------
eg.genes <- c("Rab15","Sod2","Ccl5", "Gsap",  "Mmp3",  "Il1rl1") 
#eg.genes <- c("Slc6a12", "Gsap", "Cd274", "Exoc3l4") 
#eg.genes  <- c("Bcl3", "Ccl7", "Tmem132e", "Slco3a1","Mmp9") 
pd.tc <- readRDS(file = "~/Dropbox/Projects/DurationDecoding-code/notebooks/kdeg_15m_tc.Rdata")
#pd.tc <- readRDS(file = "../../notebooks/kdeg_15m_tc_test.Rdata")
pd.tc <- pd.tc %>% 
  mutate(time=ifelse(type=="Sim.norm",time,time/60))
pd.tc$type <- factor(pd.tc$type,levels = c("Sim.norm","k_deg"))
pd.tc.2 <- pd.tc %>% 
  group_by(gene,geno,type) %>% 
  dplyr::mutate(normCnt.frac=normCnt.frac/max(normCnt.frac))

pd.tc.2$gene <- factor(pd.tc.2$gene,levels = eg.genes)
ggplot(pd.tc.2,aes(time,normCnt.frac,colour=sti,linetype=geno)) + geom_line() + 
  facet_grid(gene~type)+ scale_color_manual(values = col.map) + theme_bw()+ 
  theme(legend.position = "none")
ggsave(filename = paste(fig_folder," subfig4f_tc.eps"),height = 7,width = 4)

# add v2 (caRNA) ----------------------------------------------------------
v2 <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6b/caRNAFit-avg-sp-funs-v6-r2.csv',
         stringsAsFactors = F,row.names = 1)
v2$gene[grep('ENSMU',v2$gene)] <- dic.gene[v2$gene[grep('ENSMU',v2$gene)]]
rownames(v2) <- v2$gene
dic.v2 <- v2[dic.v1$gene,]
dic.v12 <- data.frame(v1=dic.v1$nrmsd,
                      v2=dic.v2$nrmsd)
dic.v12.fit <- dic.v12<=0.13
venn(list(v1=which(dic.v12.fit[,1]),
          v2=which(dic.v12.fit[,2])))

pd <- pd%>% mutate(v2=v2[gene,'nrmsd']<=0.13)



source('../auxilary_functions.R')
subfig_folder <- '../figures/Fig.4/subfigs/'


# Fig4B: Exp hm---------------------------------------------------------
# see preprocessing
load('../fig.4_modelv1/data/clustering_cateIII_clean.Rdata')
rord <- read.csv(file="../fig.2_knockout/data/mRNA.cat.new.csv",stringsAsFactors = F)
rord <- rord%>% filter(cate=="III")


## exp_hm.eps
ph.pd<- pd %>%filter(type=='Exp.')%>%
#ph.pd<- pd %>%filter(type=='Sim.')%>%
  unite(key,geno,type,sti,time) %>%
  spread(key,normCnt.frac)
rownames(ph.pd)  <- ph.pd$gene; ph.pd$gene<-NULL 
all(rownames(ph.pd) %in% rord$gene)
all(rord$gene%in% rownames(ph.pd) )
rord <- (rord %>% filter(gene %in% rownames(ph.pd)))$gene


bks <- seq(0,1.0001,by = .1)
anno_row <- data.frame(v1=as.factor(ph.pd[,'v1']))
anno_row_color <- list(v1=c('TRUE'='yellow','FALSE'='black'))
rownames(anno_row) <-  rownames(ph.pd)
cord <- c(15:28,1:14)+2
ph.pd.sub <-  ph.pd[rord,cord]
cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)

#postscript(file = paste0(subfig_folder,'sim_hm.eps'),
postscript(file = paste0(subfig_folder,'exp_hm.eps'),
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


# Fig4B: Sim hm---------------------------------------------------------

ph.pd<- pd %>%filter(type=='Sim.')%>%
  unite(key,geno,type,sti,time) %>%
  spread(key,normCnt.frac)
rownames(ph.pd)  <- ph.pd$gene; ph.pd$gene<-NULL 

rord <- read.csv(file="../fig.2_knockout/data/mRNA.cat.new.csv",stringsAsFactors = F)
rord <- rord%>% filter(cate=="III")

all(rownames(ph.pd) %in% rord$gene)
all(rord$gene%in% rownames(ph.pd) )
rord <- (rord %>% filter(gene %in% rownames(ph.pd)))$gene


bks <- seq(0,1.0001,by = .1)
anno_row <- data.frame(v1=as.factor(ph.pd[,'v1']))
anno_row_color <- list(v1=c('TRUE'='yellow','FALSE'='black'))
rownames(anno_row) <-  rownames(ph.pd)
cord <- c(15:28,1:14)+2
ph.pd.sub <-  ph.pd[rord,cord]
cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)

postscript(file = paste0(subfig_folder,'sim_hm.eps'),
           onefile = F,width = 2.5,height = 5)
pheatmap(ph.pd.sub,scale = "none",cluster_rows = F,
         cluster_cols = F,show_rownames = T,
         color = cols,breaks = bks,
         fontsize_row = 5,cellwidth = 4,
         #annotation_row = genes.hf,
         annotation_colors = anno_row_color,
         annotation_row = anno_row,
         legend = F,annotation_names_row = F,
         show_colnames = F,annotation_legend = F,
         gaps_col = c(7,14,21))
dev.off()


if(T){
  fig4.setting <- list(ord= ord,
                       gord = rownames(ph.pd)[ord],
                       seps= seps,
                       bks = bks,
                       cols =cols)
  fig4.setting$hm_data <- ph.pd
  saveRDS(file='fig4.setting.rds',fig4.setting)  
}


# Fig4F: specificity vs. fit -------------------------------------------------------------------
sp.th <- 2^0.5
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
pd <- data.frame(ctrl.sp= sp.mat$b1.LT.ctrl,
                 mt.sp = sp.mat$b1.LT.ko,
                 gene = rownames(sp.mat),
                 stringsAsFactors = F) 
pd <- subset(pd,gene %in% rownames(anno_row))


# explore
pd <- pd %>% mutate(V1 = anno_row[gene,])
pd <- pd %>% mutate(dynamics_dependence = ctrl.sp-mt.sp)
pd.2 <- pd%>% gather(key="Geno",value = "sp",1:2)

p.2<- ggboxplot(pd.2,x="V1",y="sp",color="V1",palette = c("black","#E7B800"),
                add="jitter",facet.by = "Geno")
p.2+theme_bw()+stat_compare_means(method.args = list(alternative="greater"),method = "t.test",label.y = 5)+
  theme(legend.position = "none",strip.background  = element_rect(fill="grey96" ))
ggsave(file=paste0(subfig_folder,'subfig4f_explore2.eps'),width = 4.2,height = 4,device=cairo_ps)

# Fig4C: examples- model vs. data -------------------------------------------------------
eg.genes <- c('Rel',"Nfkb1","Ccl1","Ccl2","Gsap",'Rab15','Mmp3')
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 
v1.simData <- read.csv(file='./models/mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 

tmp.simData <- read.csv(file = "./models/mRNA-Fit-avg-sp-v1e/bestFit.tc.csv",stringsAsFactors = F)
colnames(tmp.simData) <- sub('mRNA','normCnt.frac',colnames(tmp.simData)) 
v1.simData <- (rbind(v1.simData,tmp.simData))


p <- plotEgFit(eg.genes = eg.genes,wide = F)
p<- p + theme_bw() + theme(legend.position = "none",strip.text.x = element_blank())
p + theme(text = element_blank())
ggsave(filename = paste0(subfig_folder,'v1_eg_anno.eps'),width = 2,height = 6)
ggsave(filename = paste0(subfig_folder,'v1_eg.eps'),width = 2,height = 6)

# Fig4D: plot perturb -------------------------------------------------------------
# subset the genes passed fitting 
pd.tmp <- ph.pd[rord,]; 
target.genes <- rownames(pd.tmp)[pd.tmp$v1]; rm(pd.tmp)

# loading simulatino data 
pd.p <- readRDS('./data/kdeg_15m_new.Rdata')

pd.p <- (pd.p %>% filter(type!='Exp.')) #exclude exp. data 
pd.p$type <- sub("k_deg","k_deg_up",pd.p$type) # rename 
pd.p$gene <- as.character(pd.p$gene) 

# change name from ENSEM.. to gene symbol 
all( target.genes%in% pd.p$gene)

# plot
if(T){
  cols <- c(brewer.pal(3,"Set1")[1:2],'black')
  eps = 0.2
  pd.p<- pd.p %>% filter(gene %in% target.genes)
  pd.p$gene <- factor(pd.p$gene,levels = target.genes)
  
  pd.p <- pd.p %>%arrange(gene,type,geno)
  p<- ggplot(pd.p,
         aes(LvT.sp,ifelse(type=="Sim.norm",as.numeric(gene)-eps,
                    as.numeric(gene)+eps),
             group=interaction(gene,type))) +
    #coord_flip(xlim = c(0,length(unique(pd$gene))+1),expand = F)+
    geom_point(aes(shape=geno,colour=type))+
    geom_line(aes(colour=type))+
    #arrow=arrow(type = "closed",angle = 15,ends = "last",
    #          length = unit(.1,"inches")))  +
    scale_color_manual(name="half-life",values = cols[-1],labels=c("15m","norm."))+
    scale_shape_manual(values = c(16,1))+
    scale_y_continuous(breaks = 1:length(unique(pd.p$gene)),
                       labels = unique(pd.p$gene))+ 
    theme_bw()+theme(panel.grid.major.y = element_blank(),
                     panel.grid.minor.y = element_line(size = .5),
                     #axis.text.y =element_text(angle = 90,hjust = 1,vjust = 0.5),
                     axis.title.y = element_blank())+
    geom_vline(xintercept = .5,linetype=2,colour="blue")+
    coord_cartesian(ylim = c(0,length(unique(pd.p$gene))+1),expand = F,
                    xlim = c(-0.2,2.6)) 
  
}
ggsave(filename = paste0(subfig_folder,"subfig4e_perturb.eps"),
       width = 4.5,height = 7,scale = .75,
       p + theme(legend.position = "none")+xlab(""))

#Fig4E: examples  --------------------------------------------------------------
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
if(T){
  eg.genes  <- c("Bcl3", "Ccl7", "Rab15", "Mmp3") 
  pd.tc <- readRDS(file = "./data/kdeg_15m_tc.Rdata")
  
  pd.tc <- pd.tc %>% 
    mutate(time=ifelse(type=="Sim.norm",time,time/60))
  pd.tc$type <- factor(pd.tc$type,levels = c("Sim.norm","k_deg"))
  pd.tc.2 <- pd.tc %>% 
    filter(gene %in% eg.genes)%>%
    group_by(gene,geno,type) %>% 
    dplyr::mutate(normCnt.frac=normCnt.frac/max(normCnt.frac))
  
  pd.tc.2$gene <- factor(pd.tc.2$gene,levels = eg.genes)
  p<-ggplot(pd.tc.2,aes(time,normCnt.frac,colour=sti,linetype=geno)) + geom_line() + 
    facet_grid(gene~type)+ scale_color_manual(values = col.map) + theme_bw()+ 
    theme(legend.position = "none")
}

ggsave(filename = paste(subfig_folder," subfig4e_tc.eps"),height = 4,width = 2.6,
       p + theme(strip.text = element_blank(),text = element_blank()))





# load data ---------------------------------------------------------------
# focus on v1 unfitted genes+ caRNA 
rm(list=ls())
source(file='../auxilary_functions.R')
subfig_folder <- '../figures/Fig.5/subfigs/'
pd.fig5 <- list()


ev <- new.env();load(file="./data/caRNAFit-avg-sp-pre-v6.Rdata",ev)
cpm.peak.gene <- read.csv(file='../data/caRNA.cpm.max.geno.scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
ev<- new.env();
load(file="../fig.4_modelv1/data/clustering_cateIII_clean.Rdata",ev)

dic.clust<-unique.data.frame(ev$pd %>% dplyr::select(gene,clust))
dic.clust <- dic.clust %>% arrange(clust,gene)
rownames(dic.clust) <-dic.clust$gene
cpm.peak.gene <- cpm.peak.gene[dic.clust$gene,]
pd.ca  <- cpm.peak.gene

# combine with cytoRNA data 
ev <- new.env();
load('../fig.4_modelv1/data//clustering_cateIII_clean.Rdata',ev)
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


## Calc sp  
pd.merge.sp <- data.frame(cyto.ctrl.sp=log2(rowMax(pd.merge[,25:31])/rowMax(pd.merge[,32:38])),
                          cyto.mt.sp=log2(rowMax(pd.merge[,39:45])/rowMax(pd.merge[,46:52])),
                          ca.ctrl.sp = log2(rowMax(pd.merge[,1:6])/rowMax(pd.merge[,7:12])),
                          ca.mt.sp = log2(rowMax(pd.merge[,13:18])/rowMax(pd.merge[,19:24]))) 
pd.merge.sp <- pd.merge.sp %>%  mutate(gene=rownames(pd.merge.sp))
saveRDS(pd.merge.sp,"./data/pd.merge.sp.Rds")
# Fig5B: Sp pieplot  -------------------------------------------------------------
if(T){
  setEPS()
  postscript(file =  paste0(subfig_folder,"subfig5A_caRNA_ctrl.anno.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.ctrl.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels =table(pd.merge.sp$ca.ctrl.sp>=0.5))
  dev.off()
  
  setEPS()
  postscript(file =  paste0(subfig_folder,"subfig5A_caRNA_ctrl.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.ctrl.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels ="")
  dev.off()
  
  
  setEPS()
  postscript(file =  paste0(subfig_folder,"subfig5A_caRNA_mt.anno.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.mt.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels =table(pd.merge.sp$ca.mt.sp>=0.5))
  dev.off()
  
  setEPS()
  postscript(file =  paste0(subfig_folder,"subfig5A_caRNA_mt.eps"),width = 3,height = 3)
  pie(table(pd.merge.sp$ca.mt.sp>=0.5),col = c(gray(.8),col.map['LPS']),labels ="")
  dev.off()
}

# Fig5C: Sp scatter  ------------------------------------------------------

sp.th <- 2^0.5
selc <- c("Cgn",  "Ccl5" ,   "Fpr1") 
a <- data.frame(pd.merge.sp[,1:2],type="cyto")
b <- data.frame(pd.merge.sp[,3:4],type="ca")
colnames(a)<- colnames(b) <- c("ctrl.sp","mt.sp","type")
pd <- rbind(a,b)
pd$gene <- pd.merge.sp$gene
p<- ggplot(data.frame(ctrl.sp=c(0.5,0.5,8,8,0.5,0.5,8,8,-3,-3,0.5,0.5),
                      mt.sp=c(0.5/sp.th,8,8,8/sp.th, 
                              -3,0.5/sp.th,8/sp.th,-3,
                              -3,8,8,-3),
                      gp=factor(rep(c(2,3,1),each=4))),
           aes(ctrl.sp,mt.sp))+ 
  coord_cartesian(xlim=range(pd[,1:2]),ylim=range(pd[,1:2])) +
  geom_polygon(aes(fill=gp)) + 
  scale_fill_manual(values = alpha(c("blue"), c(.05,.1,.2)))+
  theme_bw()+
  geom_abline(slope = 1,intercept = 0,colour=grey(.8)) +
  geom_hline(yintercept = c(-.5,.5),colour=grey(.8),linetype=2)+
  geom_vline(xintercept = c(-.5,.5),colour=grey(.8),linetype=2) +
  geom_line(data = pd,aes(ctrl.sp,mt.sp,group=gene),colour=grey(.8))+
  geom_point(data = pd,aes(ctrl.sp,mt.sp,color=type),alpha=.3,size=1) + 
  annotate('point',x=pd[pd$gene%in%selc,1],y=pd[pd$gene%in%selc,2],
           colour=rep(c(col.caRNA,col.cytoRNA),3),size=2)+
  geom_line(data = pd%>% filter(pd$gene%in% selc),
            aes(ctrl.sp,mt.sp,group=gene),colour="black")+
  theme(legend.position = 'none',axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())+
  scale_color_brewer(palette = "Set1")

ggsave(file=paste0(subfig_folder,'subfig5b_sp.scatter.eps'),plot = p,width = 2.5,height = 2.5,device=cairo_ps)

# Fig5D: heatmap ---------------------------------------------------------
pd.merge.sp <- pd.merge.sp%>%
  mutate(cat=ifelse(ca.ctrl.sp<0.5,"I",ifelse(ca.mt.sp>=0.5,'II','III'))) %>%
  arrange(cat)

pd.merge.sp <- pd.merge.sp%>% mutate(ord =apply(pd.merge[pd.merge.sp$gene,1:6],1,which.max))%>%
  arrange(cat,ord)
#pd.merge.sp <- pd.merge.sp%>% arrange(cat,ca.ctrl.sp,ca.mt.sp,cyto.ctrl.sp,cyto.mt.sp)


# real plot 
bks <- seq(0,1.0001,by = 0.1)
seps <- cumsum(table(pd.merge.sp$cat))+1;

   
if(T){
  cols=colorRampPalette(brewer.pal(9,"Blues"))(10);
  setEPS()
  postscript(paste0(subfig_folder,"subfig5c_caRNA_ctrl.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,1:12], bks = bks,gaps_col = c(6),cwid = 7,
             legend = F,cols = cols)#ca_contorl 
  dev.off()
  
  setEPS()
  postscript(paste0(subfig_folder,"subfig5c_caRNA_mt.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,13:24], bks = bks,gaps_col = c(6),cwid = 7,
             legend =  F,cols = cols)#ca_contorl 
  dev.off()
  plotLegend(cols=cols,bks=bks,
             fnames = paste0(subfig_folder,"subfig5c_caRNA_scale.eps"))
}

# cytoRNA 
if(T){
  cols=colorRampPalette(brewer.pal(9,"Reds"))(10);
  setEPS()
  postscript(paste0(subfig_folder,"subfig5c_cytoRNA_ctrl.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,25:38], bks = bks,gaps_col = c(7),cwid = 7,
             legend = F,cols = cols) #ca_contorl 
  dev.off()
  
  setEPS()
  postscript(paste0(subfig_folder,"subfig5c_cytoRNA_mt.eps"),onefile = F,width = 2,height = 6)
  pheatmap.2(pd=pd.merge[pd.merge.sp$gene,39:52], bks = bks,gaps_col = c(7),cwid = 7,
             legend = F,cols = cols,rn = T,fontsize=6) #ca_contorl 
  dev.off()
  
  plotLegend(cols=cols,bks=bks,
             fnames = paste0(subfig_folder,"subfig5c_cytoRNA_scale.eps"))
  
}

# sp bar
cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)
pd.merge.sp <- pd.merge.sp %>% 
  mutate(cyto.ctrl.sp = ifelse(cyto.ctrl.sp>3,3,ifelse(cyto.ctrl.sp< -3, -3 ,cyto.ctrl.sp)))%>%
  mutate(cyto.mt.sp = ifelse(cyto.mt.sp>3,3,ifelse(cyto.mt.sp< -3, -3 ,cyto.mt.sp)))%>%
  mutate(ca.ctrl.sp = ifelse(ca.ctrl.sp>3,3,ifelse(ca.ctrl.sp< -3, -3 ,ca.ctrl.sp)))%>%
  mutate(ca.mt.sp = ifelse(ca.mt.sp>3,3,ifelse(ca.mt.sp< -3, -3 ,ca.mt.sp)))

if(T){
  nms <- c("subfig5c_cytoRNA_ctrl_sp.eps",
           "subfig5c_cytoRNA_mt_sp.eps",
           "subfig5c_caRNA_ctrl_sp.eps",
           "subfig5c_caRNA_mt_sp.eps")
  for(i in 1:4){
    setEPS()
    postscript(paste0(subfig_folder,nms[i]),
               onefile = F,width = 1.5,height = 6)
    pheatmap.2(pd=pd.merge.sp[,i],bks=bks,cols=cols,legend=F,
               cwid = 10)  
    dev.off()
  }
  plotLegend(cols = cols,bks = bks,fnames = paste0(subfig_folder,"subfig5c_sp_scale.eps"))
}

signif(bks,2)

# save order 
pd.fig5$row_names <- (pd.merge.sp$gene)
pd.fig5$row_seps <- seps

saveRDS(file = "pd.fig5.Rds",object = pd.fig5)
pd.fig5 <-readRDS(file = "./data/pd.fig5.Rds")
#pd.merge <- pd.merge[,c(1:12,25:38)] #control only 
seps.2 <- seps
cols.2 <- colorRampPalette(brewer.pal(9,"Reds"))(10)

# Fig5C: time course ------------------------------------------------------
pd <- pd.merge[selc,]%>% 
  rownames_to_column("gene")%>%
  gather(key = "sample",value="expression",2:ncol(.))%>%
  mutate(sample=sub("Exp._","",sample))%>%
  group_by(gene)%>%
  separate(col = sample,into = c("geno","sti","time"),sep = "_")%>%
  mutate(type=ifelse(geno %in% c("I","D"),"ca","cyto"))%>%
  unite("geno_sti",c("geno","sti"))
pd$geno_sti<- recode(pd$geno_sti,I_L="wt_lps",I_T="wt_tnf",D_L="mt_lps",D_T="mt_tnf")
pd$geno_sti<- factor(pd$geno_sti,levels = c("wt_lps","wt_tnf","mt_lps","mt_tnf"))
pd$type <- factor(pd$type,levels = c("cyto","ca"))
p<- ggplot(pd,aes(time,expression,group=type,color=type))+ 
  geom_point()+
  geom_line()+
  facet_grid(gene ~ geno_sti)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()

ggsave(file=paste0(subfig_folder,'sufig5d_tc_eg.eps'),
       p+theme(legend.position = "none"),
       width = 5,height = 5,device=cairo_ps)
  


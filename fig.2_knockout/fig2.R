source('../auxilary_functions.R')
subfig_dir <- '../figures/Fig.2/subfigs/'
fig_dir <- '../figures/Fig.2/'

# Fig2A: EMSA ------------------------------------------------------------------
dat.emsa.ctrl <- read.csv(file = '../data/Final_EMSA_ifnar.csv')
dat.emsa.mt <- read.csv(file = '../data/Final_EMSA_ifnarikba.csv')
colnames(dat.emsa.mt) <- sub('time',"Time",colnames(dat.emsa.mt))
dat.emsa <- rbind(dat.emsa.ctrl,dat.emsa.mt)
dat.emsa$genotype=c(rep("ctrl.",nrow(dat.emsa.ctrl)),rep('mutant',nrow(dat.emsa.mt)))
pd <-melt(dat.emsa,id=c("Time","genotype"),variable_name = "Stimuli")

pd <-pd %>% group_by(genotype)%>%mutate(value=value/max(value))
pd$Time <- pd$Time/60

pd.2 <- pd[1,]
for(g in unique(pd$genotype))
  for(s in unique(pd$Stimuli)){
    tmp.pd <- pd %>% filter(genotype==g & Stimuli==s)
    pd.2<-rbind(pd.2,data.frame(Time=seq(0,8,by=.1),
                          genotype=g,
                          Stimuli =s,
                          value=pchipfun(tmp.pd$Time,tmp.pd$value)(seq(0,8,by=.1)),stringsAsFactors = F))
}
pd.2 <- pd.2[-1,]
p <- ggplot(data = pd.2 %>% filter(Stimuli%in%c('LPS','TNF')),
            aes(x=Time, y=value,colour=Stimuli)) 
p<- p+ geom_line(aes(linetype=genotype),size=1) + theme_bw() +
  geom_point(data=pd%>% filter(Stimuli%in%c('LPS','TNF')),aes(shape=genotype),size=2) 
p<- p + xlab('Time (hr)') + scale_x_continuous(breaks = unique(pd$Time)) +ylab("") 
p<-p + ylab(expression(paste("Nuclear NF",kappa,B)))
p<-p + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = .5))
p <- p + facet_wrap(~genotype,nrow = 2)
p <- p + scale_colour_manual(values = col.map.cap)
ggsave(filename = paste0(subfig_dir,'subfig2A.eps'),
       width = 5,height = 4,p,scale=.75)
ggsave(filename = 'fig2a.png',
       width = 5,height = 4)

saveRDS(file = './data/subfig2A.rds',p)



# Fig2B: pie-chart on sp. -------------------------------------------------------------------
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
pd.venn.ctrl <- getVennData(simplify = T,sp.th.val = .5,sp.mat = sp.mat,gtype = "ko")
setEPS()
postscript(paste0(subfig_dir,"subfig1c_pie_th1_peak_gtype_","ko.eps"),
           onefile = F,width = 3,height = 3)
pie(x = pd.venn.ctrl$Counts,labels = "",col = c(gray(.8),col.map['tnf'],col.map['lps']))
dev.off()

# Fig2C: scatter of sp-------------------------------------------------------------------
sp.th <- 2^0.5
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
pd <- data.frame(ctrl.sp= sp.mat$b1.LT.ctrl,
                 mt.sp = sp.mat$b1.LT.ko,
                 gene = rownames(sp.mat),
                 stringsAsFactors = T) 
selc <- c("Il1a",  "Ccl5" ,   "Fos") 
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
  geom_point(data = pd,aes(ctrl.sp,mt.sp),alpha=.3,size=1) + 
  annotate('point',x=pd[pd$gene%in%selc,1],y=pd[pd$gene%in%selc,2],colour='red',size=2)+
  theme(legend.position = 'none',axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())

ggsave(file=paste0(subfig_dir,'subfig2f.eps'),plot = p,width = 2.5,height = 2.5,device=cairo_ps)

## save cate 
pd <- pd %>% mutate(cate=ifelse(ctrl.sp<0.5,'I',
                                ifelse(mt.sp>ctrl.sp*2^(-.5),'II','III')),
                    cate2=ifelse(ctrl.sp<0.5,'I',
                                 ifelse(mt.sp-ctrl.sp>=-.5,'II','III')))
pd$gene <- as.character(pd$gene)
pd$gene[110] <- 'Ifi203'
write.csv(pd,file = "./data/mRNA.cat.csv",quote = F)

# Fig2C: examples  ---------------------------------------------------------
rpkm_all <- read.csv(file = '../data/mRNA.nfkbgene.rpkm.csv',stringsAsFactors = F,row.names = 1,header = T)
selc <- c("Fos","Il1a",  "Ccl5" ) 
#selc <- c("Nfkb1","Nfkbie",  "Slfn2",'Rel',"Nfkb2" ) 
p<- fun.plotExp(rpkm_all,selc,scale = T)
#saveRDS(p,file='subFig2D.rds')
p<- p + theme(strip.text.x = element_blank(),legend.position = 'none',
              axis.title = element_blank(),
              axis.text =  element_blank())
p$data$gene <- factor(p$data$gene,levels = selc)         

ggsave(filename = paste0(subfig_dir,'subfig2g.eps'),plot = p,width = 2,height = 3.5,scale = .75)

# Fig2S: FP genes & TP genes -------------------------------------------------------------------
ev<- new.env(); load('../fig.4_modelv1/data/clustering_cateIII_clean.Rdata',ev)
cat3.genes <- unique(ev$pd$gene)
old_cat <- read.csv("./data/gene_cat.csv",stringsAsFactors = F)
all(cat3.genes %in% old_cat$gene.2)

selc<- setdiff((old_cat%>%filter(cate=="III"))$gene.2,cat3.genes)

## A scatter 
sp.th <- 2^0.5
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
pd <- data.frame(ctrl.sp= sp.mat$b1.LT.ctrl,
                 mt.sp = sp.mat$b1.LT.ko,
                 gene = rownames(sp.mat),
                 stringsAsFactors = T) 

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
  geom_point(data = pd,aes(ctrl.sp,mt.sp),alpha=.3,size=1) + 
  annotate('point',x=pd[pd$gene%in%selc,1],y=pd[pd$gene%in%selc,2],colour='red',size=2)+
  annotate('text',x=pd[pd$gene%in%selc,1],y=pd[pd$gene%in%selc,2]-0.2,colour='red',size=2,label=selc)+
  ylab("L vs T in mutant") + xlab("L vs T in control")+
  theme(legend.position = 'none')

ggsave(file=paste0(fig_dir,'Fig2s_A.pdf'),plot = p,width = 4,height = 4,device=cairo_pdf)

### B line graph compare 

chg_legend <- guides(colour=guide_legend(title="Stimuli"),
                     shape=guide_legend(title="Genotype"),
                     linetype=guide_legend(title="Genotype"))
p<- fun.plotExp(rpkm_all,selc,scale = T)+ylab("Expression level (0-1)")  + chg_legend
p$facet$params$ncol=3
if(T){
  pdf(file=paste0(fig_dir,'Fig2s_B.pdf'))
  tmp <- p$data %>% group_by(gene,genotype,stimuli) %>%
    summarise(m = max(value))%>% 
    ungroup() %>% 
    mutate(genotype=sub("[.]","",genotype))%>%
    group_by(gene,genotype)%>%
    summarise(sp = log2(m[stimuli=="LPS"]/m[stimuli=="TNF"]))%>% 
    spread(key = genotype,value = sp)
  
  print(p + 
    geom_text(data=data.frame(tmp,Time=6,value=.75,stimuli="LPS"),
              aes(label=paste0('ctrl:',signif(ctrl,2),"\n",'mt:',signif(mt,2))),
              colour='blue'))
  dev.off()
}




if(T){
  pdf(file=paste0(fig_dir,'Fig2s_D.pdf'))
  selc.2<- intersect((old_cat%>%filter(cate=="III"))$gene.2,cat3.genes)
  for(i in 1:9){
    p.TP<- fun.plotExp(rpkm_all,selc.2[((i-1)*9+1):(i*9)],scale = T)+ylab("Expression level (0-1)")  + chg_legend
    p.TP$facet$params$ncol=3
    tmp <- p.TP$data %>% group_by(gene,genotype,stimuli) %>%
      summarise(m = max(value))%>% 
      ungroup() %>% 
      mutate(genotype=sub("[.]","",genotype))%>%
      group_by(gene,genotype)%>%
      summarise(sp = log2(m[stimuli=="LPS"]/m[stimuli=="TNF"]))%>% 
      spread(key = genotype,value = sp)
    
    print(p.TP + 
            geom_text(data=data.frame(tmp,Time=6,value=.75,stimuli="LPS"),
                      aes(label=paste0('ctrl:',signif(ctrl,2),"\n",'mt:',signif(mt,2))),
                      colour='blue'))
  }
  dev.off()
}


# Fig2E: motif & GO  --------------------------------------------

if(F){
  gdic <- read.csv('../data/mRNA.cluster_old.csv',stringsAsFactors = F)
  rownames(gdic) <- gdic$gene
  all.equal(rownames(rpkm_all),gdic$gene)
  rownames(rpkm_all) <- gdic[rownames(rpkm_all),'gene.2']
  write.csv(rpkm_all,file = '../data/mRNA.nfkbgene.rpkm.csv')
  
  all.equal(gdic$gene,pd$gene)
  pd <- cbind(pd,gdic[,c('X','gene.2')])
  write.csv(file='mRNA.cat.csv',pd)
  
  sapply(unique(pd$cate),function(x) write.table(file =paste0('cat',x,'.txt'),(pd %>% filter(cate==x))$X,
                                               quote = F,row.names = F,col.names = F))
  # for homer
  sapply(c('I',"II","III"),FUN = function(x) write.table(pd$ensembl[pd$cate==x],file=paste0(x,'.txt'),sep = '\n',quote = F))
  all.gene <- read.delim(file='../../../DurationDecoding/data/RNAseq/MEF_count_set1/counts-roberto.txt',stringsAsFactors = F,skip = 1)
  all.gene <- substr(all.gene$Geneid,1,18)
  set.seed(100)
  sapply(c('I',"II","III"),FUN = function(x) write.table(all.gene[sample(1:length(all.gene),sum(pd$cate==x))],
                                                         file=paste0('cat',x,'-bg.txt'),sep = '\n',quote = F,
                                                         row.names = F,col.names = F))
  #findMotifs.pl catI.txt mouse ./catI -len 8,10,12,14 -start -500 -end 200 -p 8 -bg catI-bg.txt 
}


# fig2E heatmap -----------------------------------------------------------

if(F){
  pd <- read.csv("./data/mRNA.cat.csv",stringsAsFactors = F)
  #fixed 7 genes 
  rownames(pd) <- pd$gene
  pd[selc,'cate'] <- 'II'
  pd['Gdf15','cate'] <- 'I'
  write.csv(pd,"./data/mRNA.cat.new.csv")
}

sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
sp.mat <- sp.mat[,1:6]
sp.mat[sp.mat>3] <- 3; sp.mat[sp.mat < -3] <- -3 
all.equal(rownames(pd),rownames(sp.mat))

kord.sp.idx <- unlist(sapply(unique(pd$cate),function(x) which(pd$cate==x)))
kord.sp.gaps <- cumsum(sapply(unique(pd$cate),function(x) sum(pd$cate==x)))
seps <- kord.sp.gaps + 1 # for pheatmap.2 
cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)

sapply(c(1,4),
       function(i){
         setEPS()
         postscript(paste0(subfig_dir,"subfig2c_peak_sp_",i,"2.eps"),
                    onefile = F,width = 1.5,height = 6)
         pheatmap.2(pd=sp.mat[kord.sp.idx,i],bks=bks,cols=cols,legend=F,
                    cwid = 10,border_color=NA)#,gaps_col=c(1,2))  
         dev.off()
       })


### expression plot 
pd.2 <- pd # backup 
pd <- read.csv(file = '../data/mRNA.nfkbgene.rpkm.csv',stringsAsFactors = F,row.names = 1,header = T)

# not include IL1
pd.scale <- cbind(t(scale(t(pd[,1:14]))),
                  t(scale(t(pd[,22:35]))))

pd.scale[pd.scale>3] <- 3; pd.scale[pd.scale< -3] <- -3
all.equal(rownames(pd.scale),rownames(pd.2))

cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(pd.scale),max(pd.scale),length.out = 21)
sapply(c(1,15),  #c(1,22) for scaling including IL1
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(subfig_dir,"subfig2c_zscore_hm_",i,".eps"),onefile = F,width = 2,height = 6)
         dat <- pd.scale[kord.sp.idx,i:(i+13)]
         pheatmap.2(pd=dat,bks=bks,cols = cols,cwid = 5,
                    legend=F,gaps_col = grep('8',colnames(dat)),border_color=NA)
         dev.off()
       })





source(file='../auxilary_functions.R')
pd <- read.csv(file = '../data/mRNA.nfkbgene.rpkm.csv',stringsAsFactors = F,row.names = 1,header = T)
k_ord<-read.csv(file = '../data/mRNA.cluster_old.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene)
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

fig_dir <- "../figures/Fig.1/"
subfig_dir <- "../figures/Fig.1/subfigs/"

# Fig1S: rep1 vs. rep2 peak scatter  ------------------------------------------------------------------
idx.final <- c(0,grep("8",colnames(pd)))

pd.max <- sapply(2:length(idx.final), function(i) apply(pd[(idx.final[i-1]+1):idx.final[i]], 1, max))
pd.max <- log2(pd.max)
colnames(pd.max) <- gsub("X8.","",colnames(pd)[idx.final[-1]])

colnames(pd.max) <-sub("a.2","a.i.2",
    sub("a.1","a.i.1",
        sub("ifnr.2","ifnr.i.2",
      sub("ifnr.1","ifnr.i.1",
          colnames(pd.max)))))
colnames(pd.max) <-sub("[.]t",".TNF",
    sub("[.]i",".IL1",
    sub(".l",".LPS",colnames(pd.max))))

head(pd.max)
pd.max <- as.data.frame(pd.max) %>% mutate(gene=rownames(pd.max))%>% 
  gather(key = "condi",value = 'rpkm',1:ncol(pd.max))%>%
  separate(condi,into = c("geno","sti","rep"))
pd.max$geno<- recode(pd.max$geno,"ifnr"="Control","ifnrikba"="Mutant")


plots <- lapply(unique(pd.max$geno), function(x){
  pd.concord= pd.max %>% mutate(rep=paste0("r",rep))%>%
     filter(geno==x)%>%
    spread(key = rep,value = rpkm)%>% 
    group_by(sti,geno) %>% 
    summarise(concordance=ccc_fun(r1,r2))%>%
    mutate(r1=1,r2=10)
  ggplot(pd.max %>% mutate(rep=paste0("r",rep))%>%
           filter(geno==x)%>%
           spread(key = rep,value = rpkm),
         aes(x = r1,y=r2)) + geom_point() + geom_abline(slope = 1,intercept =0,color="blue")+
    facet_grid(geno~sti)  + geom_text(data = pd.concord,color="blue",
                                      aes(label=paste0("concordance=",signif(concordance,2))))+
    xlab("rep1 (log2(rpkm))") + ylab("rep2 (log2(rpkm)") + theme_bw()
  
})

plots.2 <- lapply(unique(pd.max$geno), function(x){
  pd.concord= pd.max %>% mutate(rep=paste0("r",rep))%>%
    filter(geno==x)%>%
    spread(key = rep,value = rpkm)%>% 
    group_by(sti,geno) %>% 
    summarise(concordance=cor(r1,r2))%>%
    mutate(r1=1,r2=10)
  ggplot(pd.max %>% mutate(rep=paste0("r",rep))%>%
           filter(geno==x)%>%
           spread(key = rep,value = rpkm),
         aes(x = r1,y=r2)) + geom_point() + geom_smooth(method = "lm")+
    facet_grid(geno~sti)  + geom_text(data = pd.concord,color="blue",
                                      aes(label=paste0("Pearson's r=",signif(concordance,2))))+
    xlab("rep1 (log2(rpkm))") + ylab("rep2 (log2(rpkm)") + theme_bw()
  
})
plots <- plots.2 # or comment out 
plot_grid(plots[[1]],plots[[2]],labels = c("A","B"),ncol = 1)

ggsave(filename = paste0(fig_dir,"fig1s.pdf"),width = 5,height = 4,
       scale = 1.5)
ggsave(filename = paste0(fig_dir,"fig1s.2.pdf"),width = 5,height = 4,
       scale = 1.5)

# Fig1C: max rpkm scater ------------------------------------------------------------------
pd.max.avg <- pd.max %>% group_by(gene,geno,sti) %>% summarise(rpkm_avg = mean(rpkm))
pd.max.avg.ctrl <-  pd.max.avg %>% filter(geno=="Control")%>%
  spread(key = sti,value = rpkm_avg)
fig1b <- list()
if(T){
  fig1b[[1]]<- ggplot(pd.max.avg.ctrl,aes(TNF,LPS))+geom_point() + 
    geom_abline(slope = 1,intercept = 0,color='blue')+theme_bw()
  fig1b[[2]] <- ggplot(pd.max.avg.ctrl,aes(IL1,LPS))+geom_point() + 
    geom_abline(slope = 1,intercept = 0,color="blue")+theme_bw()
  fig1b[[3]]<- ggplot(pd.max.avg.ctrl,aes(IL1,TNF))+geom_point() + 
    geom_abline(slope = 1,intercept = 0,color='blue')+theme_bw()

}
plot_grid(fig1b[[1]],fig1b[[2]],fig1b[[3]],nrow = 1)
for(i in 1:3){
  ggsave(filename = paste0(fig_dir,"/subfigs/subfig1b_",i,".png"),
         plot = fig1b[[i]] ,  
         width = 2,height = 2)
  ggsave(filename = paste0(fig_dir,"/subfigs/subfig1b_",i,".eps"),
         plot = fig1b[[i]] + theme(text = element_blank()),  
         width = 2,height = 2)
  
}

# Fig1F: hm------------------------------------------------------------------
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)

# check sp thresholding 
pd.cate <- as.numeric(cut(sp.mat$b1.LT.ctrl,c(-Inf,-0.5,0.5,Inf)))
table(pd.cate)
# check data and sp 
all.equal(rownames(sp.mat),rownames(pd)) #TRUE 

pd.scale <- cbind(t(scale(t(pd[,1:21]))),
                  t(scale(t(pd[,43:63]))))

pd.scale[pd.scale>3] <- 3; pd.scale[pd.scale< -3] <- -3
pd.ord <-apply(pd.scale[,1:7],1, which.max)
pd.ord <- (data.frame(pd.cate,pd.ord,id=seq(1:nrow(pd))) %>% arrange(pd.cate,pd.ord))
seps<- sapply(1:3, function(i) max(which(pd.ord$pd.cate==i)))
cols <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
bks <- seq(-max(pd.scale),max(pd.scale),length.out = 21)

pheatmap(pd.scale[pd.ord$id,],scale = "none",cluster_cols = F,cluster_rows = F,color = cols,
         gaps_row  = seps,breaks = bks,gaps_col = grep("8",colnames(pd.scale)))

seps <- seps+1

sapply(1:2, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(subfig_dir,"z_rpkm_genoScale_",i,".eps"),onefile = F,width = 2,height = 6)
         dat <- pd.scale[pd.ord$id,((i-1)*21+1):(i*21)]
         pheatmap.2(pd=dat,bks=bks,cols = cols,cwid = 5,
                    legend=F,gaps_col = grep('8',colnames(dat)))
         dev.off()
       })

plotLegend(cols = cols,
           bks = bks,
           fnames = paste0(subfig_dir,"z_rpkm_batchScale_lg.eps"))

# Fig1C: sp scatter --------------------------------------------------------

### sp plot 
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
sp.mat <- sp.mat[,grepl('ctrl',colnames(sp.mat))]
sp.mat.log2 <- sp.mat
pd<- data.frame(LvT.sp = sp.mat$b1.LT.ctrl,
                LvI.sp = sp.mat$b1.LI.ctrl)

fig1c<- ggplot(pd,aes(LvT.sp,LvI.sp)) + geom_point() +
   geom_abline(slope = 1,intercept = 0,color="blue") + 
  theme_bw()
fig1c +  annotate("text",0,7.5, 
            label=paste0("concordance=",signif(ccc_fun(pd[,1],pd[,2]),2)),colour='blue')
ggsave(paste0(subfig_dir,'fig1c_scatter.png'),width = 2.5,height = 2.5)
fig1c +  theme(text = element_blank())
ggsave(paste0(subfig_dir,'fig1c_scatter.eps'),width =2,height = 2)


# Fig1F: sp heatmap  -----------------------------------------------------------------
sp.mat.log2[sp.mat.log2>3] <- 3 ; sp.mat.log2[sp.mat.log2 < -3] <- -3 
all.equal(rownames(sp.mat.log2),row.names(pd.scale))# True 

cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
          colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3)

sapply(1:2,
       function(i){
         setEPS()
         postscript(paste0(subfig_dir,"sp_maxFrac_genoScale2_",i,".eps"),
                    onefile = F,width = 3,height = 6)
         pheatmap.2(pd=sp.mat.log2[pd.ord$id,(i*3-2):(i*3)],bks=bks,cols=cols,legend=F,
                    cwid = 10,gaps_col=1:2)  
         dev.off()
       })


plotLegend(cols = cols,bks = bks,fnames = paste0(subfig_dir,"sp_zrpkm_genoScale_lg.eps"))
signif(bks,2)


# fig1.c Pie------------------------------------------------------------------
pd.venn.ctrl <- getVennData(simplify = T,sp.th.val = .5,sp.mat = sp.mat,gtype = "ctrl")
setEPS()
postscript(paste0(fig_folder,"subfig1c_pie_th1_",type,"_gtype_","ctrl.eps"),
           onefile = F,width = 3,height = 3)
pie(x = pd.venn.ctrl$Counts,labels = "",col = c(gray(.8),col.map['tnf'],col.map['lps']))
dev.off()


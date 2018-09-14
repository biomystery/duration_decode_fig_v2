source('../auxilary_functions.R')
fig_folder <- '~/Dropbox/grant/2017_RO1_GeneExpression/figures/Fig.2.1/'
# caRNA:heatmap by group  -------------------------------------------------
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)
#sp.mRNA <- read.csv(file="../data/mRNA.sp.avg.csv",header = T,
 #                   row.names = 1,stringsAsFactors = F)

cateogry = ifelse(sp.mRNA$b1.LT.ctrl<0.5,'non-sp.',ifelse(sp.mRNA$b1.LT.ko-sp.mRNA$b1.LT.ctrl<=-0.5,
                                                          "dy","LPS-sp-no-dy")) #no-specific, LPS-sp but not dy, LPS-sp &dy 
rev(table(cateogry))
names(cateogry )<- rownames(sp.mRNA)


cpm.avg.gene <- read.csv(file='../data/caRNA.cpm.max.geno.scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)
sum(rownames(cpm.avg.gene) %in% k_ord$gene.2)
rownames(cpm.avg.gene) <- k_ord$gene
all.equal(rownames(cpm.avg.gene),rownames(sp.mRNA))

pd <- cpm.avg.gene
#pd.new<- rbind(pd[cateogry=='non-sp.',],
#                     pd[cateogry=='LPS-sp-no-dy',],
#                     pd[cateogry=='dy',])

pd.new<- pd[cateogry=='dy',]
seps <- 0
#seps <- cumsum(rev(table(cateogry)))+1 # order: non-sp, lps-sp-no-dy,dy


bks <- seq(0,1,by = 0.1)
cols <- colorRampPalette(brewer.pal(9,"Blues"))(10)

require(pheatmap)
pheatmap(pd.new,cluster_cols = F,cluster_rows = F,scale = 'none',show_rownames = F,color = cols,
         breaks = bks,border_color = NA)

sapply(1:2, 
       function(i){
         #i <-2
         setEPS()
         postscript(paste0(fig_folder,"caRNA_max_bycat",i,".eps"),onefile = F,width = 2,height = 6)
         pheatmap.2(pd=pd.new[,((i-1)*12+1):(i*12)], bks = bks,gaps_col = c(6),cwid = 7,
                    legend = F,cols = cols,border_color = NA)
         dev.off()
       })

plotLegend(cols=cols,bks=bks,
           fnames = paste0(fig_folder,"caRNA_max_scale.eps"))

# sp. caRNA 

sp.caRNA <- read.csv('../data/caRNA.sp.peak.csv',stringsAsFactors = F,row.names = 1)
data.frame(rownames(sp.caRNA),rownames(pd),names(cateogry))

sp.caRNA.log2.new <-                            sp.caRNA[cateogry=='dy',]


cols <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))
bks <- c(-3,-2,-1,-.5,0,.5,1,2,3.001)

sp.caRNA.log2.new[sp.caRNA.log2.new>3] <- 3
gaps <- 40
sapply(1:2,
       function(i){
         setEPS()
         postscript(paste0(fig_folder,"sp.avg.caRNA.cat",i,".eps"),
                    onefile = F,width = 1.5,height = 6)
         pheatmap.2(pd=sp.caRNA.log2.new[,i],bks=bks,cols=cols,legend=F,
                    cwid = 10,border_color = NA)  
         dev.off()
       })


# mRNA heatmap ------------------------------------------------------------
pd <- read.csv(file='../data/mRNA.nfkbgene.rpkm.csv',header = T,row.names = 1)
pd <- pd[,1:42]
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)
all.equal(rownames(pd),k_ord$gene)

# max scale 
pd.scale <- cbind(t(apply(pd[,1:21],1,function(x) x/max(x))),
                  t(apply(pd[,22:42],1,function(x) x/max(x))))


# category 
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)

all.equal(rownames(sp.mRNA),rownames(pd.scale))
cateogry = ifelse(sp.mRNA$b1.LT.ctrl<0.5,'non-sp.',ifelse(sp.mRNA$b1.LT.ko-sp.mRNA$b1.LT.ctrl<=-0.5,
                                                          "dy","LPS-sp-no-dy")) #no-specific, LPS-sp but not dy, LPS-sp &dy 
rev(table(cateogry))
names(cateogry )<- rownames(sp.mRNA)

pd.scale.new<-  pd.scale[cateogry=='dy',]
pd.scale.new <- pd.scale.new[,-c(15:21,36:42)]

cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)
bks <- seq(0,1.0001,by = .1)
seps <- 0
sapply(1:2, 
       function(i){
         #i <-1
         setEPS()
         postscript(paste0(fig_folder,"fig2s_mRNA_by_cat",i,".eps"),onefile = F,width = 2,height = 6)
         dat <- pd.scale.new[,((i-1)*14+1):(i*14)]
         pheatmap.2(pd=dat,bks=bks,cols = cols,cwid = 5,
                    legend=F,gaps_col = grep('8',colnames(dat)),border_color = NA)
         dev.off()
       })

### sp plot 
all.equal(rownames(pd.scale),rownames(sp.mRNA))
all.equal(rownames(pd.scale),names(cateogry))
sp.mRNA.new<-    sp.mRNA[cateogry=='dy',]


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
                    cwid = 10,border_color = NA)#,gaps_col=c(1,2))  
         dev.off()
       })

signif(bks,2)

cat(rownames(pd.scale.new),sep = '\t')


# sp caRNA vs cytoRNA -----------------------------------------------------

all.equal(rownames(sp.caRNA.log2.new),rownames(sp.mRNA.new))
data.frame(rownames(sp.caRNA.log2.new),rownames(sp.mRNA.new))

pd <- data.frame(sp.ca = sp.caRNA.log2.new$LT.ctrl,sp.mRNA = sp.mRNA.new$b1.LT.ctrl,
                 gene=rownames(sp.caRNA.log2.new))
require(plotly)
ggplotly(ggplot(pd,aes(sp.ca,sp.mRNA))+ geom_point(aes(text=gene)),
         tooltip=c("sp.ca","sp.mRNA",'text'))

plot(sp.caRNA.log2.new$LT.ctrl,sp.mRNA.new$b1.LT.ctrl,xlim=c(-1,3),ylim=c(-1,3))
abline(a=c(-1,1),b=1)

plot(sp.caRNA.log2.new$LT.mt,sp.mRNA.new$b1.LT.ko,xlim=c(-1,3),ylim=c(-1,3))
abline(a=1,b=1)




# mRNA half-life

genes.hf <- read.csv(file="../data/v4-hf-final.csv",row.names = 1,stringsAsFactors = F)
all(rownames(sp.mRNA.new) %in% rownames(genes.hf))
pd.hf <- data.frame(hf = genes.hf[rownames(sp.mRNA.new),],
                    group = c(rep(1,40),rep(2,nrow(sp.mRNA.new)-40)))

tmp <- sp.mRNA.new$b1.LT.ctrl - sp.caRNA.log2.new$LT.ctrl
pd.hf <- data.frame(hf = genes.hf[rownames(sp.mRNA.new),],
                    group = cut(tmp,c(min(tmp,na.rm = T)-.001,0.5,max(tmp,na.rm = T)+0.001)))

pd.hf <- pd.hf%>%gather(key = sti_time,hf,-ncol(pd.hf)) %>% 
  filter(sti_time %in% c("hf.LPS0","hf.LPS30","hf.TNF30"))
pd.hf$group<- droplevels(pd.hf$group)
pd.hf <- pd.hf[-which(is.na(pd.hf$group)),]
ggplot(pd.hf,aes(group,hf)) + geom_boxplot(aes(fill=sti_time))

ggplot(pd.hf,aes(group,hf)) + geom_boxplot(aes(fill=sti_time),outlier.shape=NA) +  
  scale_y_continuous(limits = quantile(na.omit(pd.hf$hf), c(0.1, 0.9)))

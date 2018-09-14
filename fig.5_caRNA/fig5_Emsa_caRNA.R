rm(list=ls())
source(file='../auxilary_functions.R')
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.5_caRNA/subfigs/'


# cmp v1 & v2:nRMSD -------------------------------------------------------------
ev <- new.env(); load('../fig.4_modelv1/clustering_cateIII_clean.Rdata',ev)
nrmsd.v1 <- ev$pd %>% group_by(gene)%>%
  spread(type,normCnt.frac) %>% mutate(res=Sim.-Exp.)%>%
  summarise(nrmsd=sqrt(sum(res^2)/length(res))/(max(Exp.)-min(Exp.)))
nrmsd.v2 <- read.csv(file="caRNAFit-avg-sp-v6b/caRNAFit-avg-sp-funs-v6-r2.csv",
                     stringsAsFactors = F,row.names = 1)
dic.gene<- getSymbol(grep('ENSM*',unique(nrmsd.v2$gene),value = T)) # deal with the name 
nrmsd.v2$gene[grep('ENSM*',nrmsd.v2$gene)]<- dic.gene[nrmsd.v2$gene[grep('ENSM*',nrmsd.v2$gene)]]
rownames(nrmsd.v2)[grep('ENSM*',rownames(nrmsd.v2))]<- dic.gene[rownames(nrmsd.v2)[grep('ENSM*',rownames(nrmsd.v2))]]
sum(nrmsd.v1$gene %in% nrmsd.v2$gene) - nrow(nrmsd.v1)

nrmsd.v12 <- data.frame(nrmsd_v1=nrmsd.v1$nrmsd,
                        nrmsd_v2=nrmsd.v2[nrmsd.v1$gene,]$nrmsd,
                        gene = nrmsd.v1$gene)

p <- ggplot(nrmsd.v12,aes(nrmsd_v1,nrmsd_v2,label=gene)) + geom_point()+ 
  geom_abline(slope = 1,intercept = 0)+geom_hline(yintercept = .13)+geom_vline(xintercept = .13)
require(plotly); ggplotly(p)

# v2c (use exon reads for caRNA )
nrmsd.v2c <- read.csv(file="caRNAFit-avg-sp-v6c/caRNAFit-avg-sp-funs-v6-r2.csv",
                     stringsAsFactors = F,row.names = 1)
nrmsd.v2 <- nrmsd.v2c
nrmsd.v2$gene[grep('ENSM*',nrmsd.v2$gene)]<- dic.gene[nrmsd.v2$gene[grep('ENSM*',nrmsd.v2$gene)]]
rownames(nrmsd.v2)[grep('ENSM*',rownames(nrmsd.v2))]<- dic.gene[rownames(nrmsd.v2)[grep('ENSM*',rownames(nrmsd.v2))]]
nrmsd.v12$nrmsd_v2c <- nrmsd.v2[nrmsd.v1$gene,]$nrmsd
p <- ggplot(nrmsd.v12,aes(nrmsd_v2,nrmsd_v2c,label=gene)) + geom_point()+ 
  geom_abline(slope = 1,intercept = 0)+geom_hline(yintercept = .13)+geom_vline(xintercept = .13)

# save caRNA exon value into csv file
ev <- new.env();load(file="./caRNAFit-avg-sp-v6c/caRNAFit-avg-sp-pre-v6.Rdata",ev)
dic.gene<- getSymbol(grep('ENSM*',rownames(ev$caRNA.frac.avg),value = T)) # deal with the name 
rownames(ev$caRNA.frac.avg)[grep('ENSM*',rownames(ev$caRNA.frac.avg))]<- dic.gene[rownames(ev$caRNA.frac.avg)[grep('ENSM*',rownames(ev$caRNA.frac.avg))]]

write.csv(file='./caRNAFit-avg-sp-v6c/caRNA.cpm.max.geno.scale.csv',ev$caRNA.frac.avg)


# cmp v1 & v2:trajectory -------------------------------------------------------------
v1.simData <- read.csv(file='../fig.4_modelv1//mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
v1.simData.new <- read.csv(file="../fig.4_modelv1//mRNA-Fit-avg-sp-v1e/bestFit.tc.csv",
                           stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
names(v1.simData.new) <- sub("mRNA","normCnt.frac",names(v1.simData.new)) 
v1.simData <- rbind(v1.simData,v1.simData.new) ; rm(v1.simData.new)
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 


# v1 plot 
ev.1 <- ev
rownames(nrmsd.v12) <- nrmsd.v12$gene; nrmsd.v12$gene <- NULL 
if(T){
  p.genes <- c('Fpr1','Lcn2','Tnfsf11','Ccl1')
  p.genes <- unique(ev.1$pd$gene[!ev.1$pd$v1 & ev.1$pd$clust==2])
  v1.simData$gene[grep('ENSM*',v1.simData$gene)] <- dic.gene[v1.simData$gene[grep('ENSM*',v1.simData$gene)]]
  maxScale.mRNA$gene[grep('ENSM*',maxScale.mRNA$gene)] <- dic.gene[maxScale.mRNA$gene[grep('ENSM*',maxScale.mRNA$gene)]]
  
  names(col.map ) <- tolower(names(col.map))
  p.v1<- plotEgFit(eg.genes = p.genes,simData = v1.simData,
                   expData = maxScale.mRNA,isFinal = F,wide = T)
  
  # v2 data 
  v2.simData <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6b/bestFit.tc.csv',
                         stringsAsFactors = F)
  names(v2.simData) <- sub("mRNA","normCnt.frac",names(v2.simData)) 
  v2.simData$gene[grep('ENSM*',v2.simData$gene)] <- dic.gene[v2.simData$gene[grep('ENSM*',v2.simData$gene)]]
  p.2 <- p.v1 %+% rbind(p.v1$data,data.frame(v2.simData%>%filter(gene %in%p.genes),type='v2_Sim.'))
  
  
  # v2 exp data
  ev <- new.env()
  load("../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-pre-v6.Rdata",ev)
  rownames(ev$caRNA.frac.avg)[rownames(ev$caRNA.frac.avg)%in%names(dic.gene)] <- dic.gene[rownames(ev$caRNA.frac.avg)[rownames(ev$caRNA.frac.avg)%in%names(dic.gene)]]
  pd.ca.exp <- ev$caRNA.frac.avg[p.genes,]
  pd.ca.exp <- pd.ca.exp%>% mutate(gene=rownames(pd.ca.exp))%>%
    gather(key = sample,value = normCnt.frac,1:24)%>%
    mutate(time=sapply(sample,function(x) as.numeric(unlist(strsplit(x,split = "_"))[3])),
           sti = sapply(sample,function(x) unlist(strsplit(x,split = "_"))[2]),
           geno= sapply(sample,function(x) unlist(strsplit(x,split = "_"))[1]),
           type="v2_input")
  geno.dic <- c(I="wt",D="mt");sti.dic <- c(L='lps',T='tnf')
  
  pd.ca.exp$sti <-sti.dic[pd.ca.exp$sti]
  pd.ca.exp$geno <-geno.dic[pd.ca.exp$geno]
  
  pd.ca.exp$gene <- factor(pd.ca.exp$gene,levels = p.genes)
  
  p.3 <- p.2 + geom_point(data = pd.ca.exp,
                   aes(time,normCnt.frac,shape=geno,colour=sti)) +
    geom_line(data = pd.ca.exp,
              aes(linetype=geno))
  
  # v2c, do it again as v2
  v2.simData <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6c/bestFit.tc.csv',
                         stringsAsFactors = F)
  names(v2.simData) <- sub("mRNA","normCnt.frac",names(v2.simData)) 
  v2.simData$gene[grep('ENSM*',v2.simData$gene)] <- dic.gene[v2.simData$gene[grep('ENSM*',v2.simData$gene)]]
  #p.3 <- p.3 %+% rbind(p.3$data,data.frame(v2.simData%>%filter(gene %in%p.genes),type='v2c_Sim.'))
  p.3 <- p.2 %+% rbind(p.2$data,data.frame(v2.simData%>%filter(gene %in%p.genes),type='v2c_Sim.'))
  
  ev <- new.env()
  load("../fig.5_caRNA/caRNAFit-avg-sp-v6c/caRNAFit-avg-sp-pre-v6.Rdata",ev)
  rownames(ev$caRNA.frac.avg)[rownames(ev$caRNA.frac.avg)%in%names(dic.gene)] <- dic.gene[rownames(ev$caRNA.frac.avg)[rownames(ev$caRNA.frac.avg)%in%names(dic.gene)]]
  pd.ca.exp <- ev$caRNA.frac.avg[p.genes,]; pd.ca.exp <- as.data.frame(pd.ca.exp)
  pd.ca.exp <- pd.ca.exp%>% mutate(gene=rownames(pd.ca.exp))%>%
    gather(key = sample,value = normCnt.frac,1:24)%>%
    mutate(time=sapply(sample,function(x) as.numeric(unlist(strsplit(x,split = "_"))[3])),
           sti = sapply(sample,function(x) unlist(strsplit(x,split = "_"))[2]),
           geno= sapply(sample,function(x) unlist(strsplit(x,split = "_"))[1]),
           type="v2c_input")
  geno.dic <- c(I="wt",D="mt");sti.dic <- c(L='lps',T='tnf')
  
  pd.ca.exp$sti <-sti.dic[pd.ca.exp$sti]
  pd.ca.exp$geno <-geno.dic[pd.ca.exp$geno]
  
  pd.ca.exp$gene <- factor(pd.ca.exp$gene,levels = p.genes)
  
  p.4 <- print(p.3 + geom_point(data = pd.ca.exp,
                          aes(time,normCnt.frac,shape=geno,colour=sti)) +
    geom_line(data = pd.ca.exp,
              aes(linetype=geno)))
  
  
  
  # print 
  print(nrmsd.v12[p.genes,],digits = 2)
}
p.3
p.4 %+% (p.4$data%>% filter(type %in% c("Exp.","Sim.","v2_Sim.","v2c_Sim.")))

# updated to new ensemble exon reads  -------------------------------------
# genes in catIII 
gene.dic <- read.csv(file='../data/mRNA.cluster_old.csv',stringsAsFactors = F)
ev.1 <- new.env(); load('../fig.4_modelv1/clustering_cateIII_clean.Rdata',ev.1)
gene.dic <- subset(gene.dic,gene.2 %in% ev.1$pd$gene)
all(unique(ev.1$pd$gene %in% gene.dic$gene.2))

# exon reads rpkm 
require(edgeR) 
reads.exon <- read.delim(file='./data/new_ensembled_reads/counts-exon.txt',skip = 1,
                         stringsAsFactors = F)
length.exon <- reads.exon$Length
reads.exon.rpkm <- rpkm(reads.exon[,7:ncol(reads.exon)],gene.length = reads.exon$Length)
rownames(reads.exon.rpkm) <- substr(reads.exon$Geneid,1,18); rm(reads.exon)
reads.exon.rpkm <- reads.exon.rpkm[gene.dic$X,]


# whole genes reads rpkm 
reads.gene <- read.delim(file='./data/new_ensembled_reads/counts-whole-gene.txt',skip = 1,
                         stringsAsFactors = F)
length.gene <- reads.gene$Length
reads.gene.rpkm <- rpkm(reads.gene[,7:ncol(reads.gene)],gene.length = reads.gene$Length)
rownames(reads.gene.rpkm) <- substr(reads.gene$Geneid,1,18); rm(reads.gene)
reads.gene.rpkm <- reads.gene.rpkm[gene.dic$X,]

# compare length 
plot(log10(length.gene),log10(length.exon)); abline(a=0,b=1)

# comapre rpkm 
ggplot(data.frame(exon=matrix(reads.exon.rpkm,ncol = 1),
                  gene= matrix(reads.gene.rpkm,ncol = 1)),
       aes(log2(exon+1),log2(gene+1)))+ geom_point()+ 
  geom_abline(slope = 1,intercept = 0,linetype=2) + ggtitle("rpkm compare")

# exon row max rpkm 
par(mfrow=c(1,2))
pd <- log2(rowMax(reads.exon.rpkm)+1)
hist(pd,xlab = "log2 rpkm+1",main='exon')
abline(v=quantile(pd),lty=2)
pd <- log2(rowMax(reads.gene.rpkm)+1)
hist(pd,xlab = "log2 rpkm+1",main='gene')
abline(v=quantile(pd),lty=2)
par(mfrow=c(1,1))

# save the data for model v6c fitting 
gene.dic <- read.csv(file='../data/mRNA.cluster_old.csv',stringsAsFactors = F) # all genes 
ev <- new.env();load(file="./caRNAFit-avg-sp-v6b/caRNAFit-avg-sp-pre-v6.Rdata",ev)
all(rownames(ev$caRNA.frac.avg) %in% gene.dic$gene) #TRUE 
cord <- c(12:17,12,18:22,1:6,1,7:11);
data.frame(colnames(reads.exon.rpkm)[cord] ,
           colnames(ev$caRNA.frac.avg))
reads.exon.rpkm <-reads.exon.rpkm[,cord] 
all(rownames(reads.exon.rpkm)==gene.dic$X) # TRUE 
rownames(reads.exon.rpkm) <- gene.dic$gene # same rowname set
colnames(reads.exon.rpkm) <- colnames(ev$caRNA.frac.avg) # same colname
reads.exon.rpkm <- reads.exon.rpkm[rownames(ev$caRNA.frac.avg),] # same rowname order 
colnames(reads.exon.rpkm[,1:12])
reads.exon.rpkm[,1:12] <- reads.exon.rpkm[,1:12]/apply(reads.exon.rpkm[,1:12],1,max)
reads.exon.rpkm[,13:24] <- reads.exon.rpkm[,13:24]/apply(reads.exon.rpkm[,13:24],1,max)
sum(is.na(reads.exon.rpkm)) # 0 good 
range(reads.exon.rpkm) # 0 - 1 
ev$caRNA.frac.avg <- reads.exon.rpkm # replace old matrix 
save(file = "caRNAFit-avg-sp-pre-v6.Rdata",list = ls(ev),envir = ev)

pheatmap(cbind(ev$caRNA.frac.avg,reads.exon.rpkm),show_rownames = F,
               scale = 'none',cluster_cols = T,gaps_col = seq(0,48,by = 6))

pd <- data.frame(gene=matrix(unlist(ev$caRNA.frac.avg),ncol = 1),
                 exon=matrix(reads.exon.rpkm,ncol = 1))
ggplot(pd,aes(gene,exon)) + 
  geom_smooth(method = "lm") +  stat_density2d(aes(fill = ..density..), 
                                               geom = "tile", contour = FALSE, n = 100) +
  scale_fill_continuous(low = "white", high = "red")+geom_point(alpha=0.2)+ 
  theme_bw()+ 
  coord_cartesian(xlim = c(-0.01,1.01),ylim = c(-0.01,1.01),expand = F)+
  ggtitle(paste('cor=',signif(cor(pd)[2],2)))

# cmp v1 & v2 -------------------------------------------------------------
r2.v3 <- read.csv(file='./emsa-caRNA-v1/v1-cmp-N.csv',stringsAsFactors = F,row.names = 1)
r2.v2 <-  read.csv(file='./caRNAFit-avg-sp-v6/caRNAFit-avg-sp-funs-v6-r2.csv',
                   row.names = 1,stringsAsFactors = F)
r2.v1 <- read.csv(file = '../fig.4_modelv1/mRNA-Fit-avg-sp-v1c/result.csv',
                   row.names = 1,stringsAsFactors = F)
all.equal(rownames(r2.v1),rownames(r2.v3),rownames(r2.v2))

r2.cmp <- data.frame(r2.v1=r2.v1$r2,
                     r2.v2=r2.v2$r2,
                     r2.v3=r2.v3$r2,
                     gene = rownames(r2.v1),stringsAsFactors = T)

# v1 & v3
r2.th = 0.7 
require(venn)
venn(x = list(v1.fit = which(r2.cmp[,1]>=r2.th),
              v3.fit = which(r2.cmp[,3]>=r2.th)),
     transparency = 1,
     cexil=2,cexsn = 4)

# cmp v1 & v2 & v3 -------------------------------------------------------------


r2.cmp.rescued <- r2.cmp %>% filter(r2.v3>=r2.th & r2.v1< r2.th) %>% arrange(desc(r2.v3))
r2.cmp.rescued %>% filter(r2.v2>=r2.th) # 10  genes 

r2.cmp.notrescued <- r2.cmp %>% filter(r2.v3<r2.th & r2.v1< r2.th)
cat( as.character(r2.cmp.rescued$gene),sep = '\n')
cat( as.character(r2.cmp.notrescued$gene),sep = '\t')

tmp <- r2.cmp %>% filter(r2.v3<r2.th & r2.v1>= r2.th)
cat(as.character(tmp$gene),sep='\t')
r2.cmp %>% filter(r2.v3>=r2.th,r2.v2<r2.th,r2.v1<r2.th) %>% arrange(desc(r2.v2))
r2.cmp.rescued<- r2.cmp.rescued %>% arrange(desc(r2.v3))
require(tidyr)
require(plotly)
 p <- ggplot(r2.cmp %>% gather(key=model,value=r2,1:3),aes(model,r2,group=gene)) +
   geom_point(aes(colour=model))+
   geom_line()+
  geom_hline(yintercept = .7,colour='gold') 
 
p$data <- r2.cmp %>% gather(key=model,value=r2,c(1,3))
p %+% geom_line(data=subset(r2.cmp,r2.v1<r2.th & r2.v3>=r2.th)%>% gather(key=model,value=r2,c(1,3)),
                colour='blue')
p$data <- subset(r2.cmp,r2.v1<r2.th & r2.v3>=r2.th)%>% gather(key=model,value=r2,c(1,2,3))
p 
ggplotly()

# plot the tc for the rescued genes ---------------------------------------
v3.simData <- read.csv(file='./emsa-caRNA-v1/bestFit.tc.csv',stringsAsFactors = F)
names(v3.simData) <- sub("caRNA","normCnt.frac",names(v3.simData)) 
v1.simData <- read.csv(file='../fig.4_modelfit/mRNA-Fit-avg-sp-v1c/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 

maxScale.caRNA <- plotExp(dtype = 'caRNA',scale = "geno",savetofile = F)
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <-maxScale.caRNA$cond <- NULL
maxScale.mRNA$sti <- toupper(maxScale.mRNA$sti )


# plot
p.genes <- as.character(r2.cmp.rescued$gene[1:5+30])

p.v3<- plotEgFit(eg.genes = p.genes,simData = v3.simData,
          expData = maxScale.caRNA,isFinal = F)
pd_lab <- data.frame(r2=r2.v3[p.genes,'r2'],
                     gene=p.genes,
                     type=rep('Sim.',length(p.genes)),
                     stringsAsFactors = F)
pd_lab$gene <- factor(pd_lab$gene,levels = p.genes)
p.lab.v3 <- p.v3 + geom_text(data= pd_lab,
                   aes(x = 6,y=1,label=paste0('R2=',signif(r2,2))),colour='red',
                   hjust=0,vjust=1) 
p.lab.v3

# v1 
p.v1<- plotEgFit(eg.genes = p.genes,simData = v1.simData,
              expData = maxScale.mRNA,isFinal = F)

pd_lab <- data.frame(r2=r2.v1[p.genes,'r2'],
                     gene=p.genes,
                     type=rep('Sim.',length(p.genes)),
                     stringsAsFactors = F)
pd_lab$gene <- factor(pd_lab$gene,levels = p.genes)
p.lab.v1 <- p.v1 + geom_text(data= pd_lab,
                       aes(x = 6,y=1,label=paste0('R2=',signif(r2,2))),colour='red',
                       hjust=0,vjust=1) 
p.lab.v1




# fig.5A - sp comparason (old, use mean) -------------------------------------------------

sp.caRNA <- read.csv(file='./data/caRNA-max-avg-sp.csv',header = T,
                     row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)

sum(rownames(sp.caRNA) %in% k_ord$gene.2)
sp.caRNA <- sp.caRNA[k_ord$gene.2,]
rownames(sp.caRNA) <- k_ord$gene

caRNA.sum<- data.frame(ctrl.L = rowSums(cpm.avg.gene[,1:6]),
                       ctrl.T = rowSums(cpm.avg.gene[,7:12]),
                       mt.L = rowSums(cpm.avg.gene[,13:18]),
                       mt.T = rowSums(cpm.avg.gene[,19:24]))

sp.caRNA.2 <- data.frame(LT.ctrl=caRNA.sum$ctrl.L/caRNA.sum$ctrl.T,
                         LT.mt = caRNA.sum$mt.L/caRNA.sum$mt.T)
sp.caRNA.2[apply(sp.caRNA.2,2,is.infinite)] <- NA
sp.caRNA.2[sp.caRNA.2==0] <- NA
sp.caRNA.log2 <- log2(sp.caRNA.2)
rownames(sp.caRNA.log2) <- rownames(caRNA.sum)
all.equal(rownames(caRNA.sum),rownames(sp.mRNA))
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
all.equal(rownames(sp.mRNA.log2),rownames(sp.caRNA.log2))

sp.mRNA.log2 <- log2(sp.mRNA)
sp.mRNA.log2[sp.mRNA.log2>3] <- 3
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


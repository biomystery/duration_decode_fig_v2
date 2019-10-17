source("../auxilary_functions.R")

# export parameters  ------------------------------------------------------
load('../fig.4_modelv1/data/clustering_cateIII_clean.Rdata')
pars.v1 <- read.csv(file='./models/mRNA-Fit-avg-sp-v1e/result-RMSD.csv',stringsAsFactors = F)
pars.v1.2 <- read.csv(file='./models/mRNA-Fit-avg-sp-v1d/result-RMSD.csv',stringsAsFactors = F)

pars.v1.all <- rbind(pars.v1[,-c(1,6)],pars.v1.2[-(76:88),-1])

dic.gene<- read.csv('../data/mRNA.cluster_old.csv',stringsAsFactors = F)[,c("gene","gene.2")]%>%
  column_to_rownames('gene')
pars.v1.all$gene =dic.gene[pars.v1.all$gene,'gene.2']

# R2 histgram -------------------------------------------------------------

## R2 histgram
ggplot(anno.row%>% mutate(r2=ifelse(r2<0,0,r2)),aes(r2,..count..))+ geom_histogram(binwidth = .1,colour='white')+
  theme_bw()
ggsave(filename = paste0(fig_folder,'subfig.r2hist.eps'),width = 3,height = 3)

## plot sim vs. exp : half-life
pd.hf <- data.frame(hf=log(2)/fit.r2[ph.pd$gene,"k_deg"]/60,
                    r2=fit.r2[ph.pd$gene,'r2'])
bks.r2 <- c(min(pd.hf)-0.001,0.5,5,33)
bks.r2 <- c(.15,.5,1,3,5,16,max(pd.hf)+0.01)
cols.r2 <- colorRampPalette(c("azure","chartreuse4"))(length(bks.r2)-1)
postscript(paste0(fig_folder,"fitHf.eps"),
           onefile = F,width = 1,height = 6)
pheatmap.2(pd=pd.hf[ord,1],bks=bks.r2,cols=cols.r2,legend=F,
           cwid = 10)  
dev.off()

plotLegend(cols = cols.r2,bks = 1:length(bks.r2),fnames = paste0(fig_folder,"fitHf_hm_scale.eps")) 


# half-life compare
ggplot(pd.hf%>%mutate(isGoodFit=r2>=.7),aes(isGoodFit,hf,fill=isGoodFit)) + geom_boxplot() + scale_fill_brewer(palette = "Set1")+
  ylab('Half-life (hr)')


anno.row$r2<- cut(anno.row$r2,breaks = bks.r2)
anno.row.col <- list(r2=c("(-0.753,0]"=cols.r2[1],
                          "(0,0.5]" =cols.r2[2],
                          "(0.5,0.6]" =cols.r2[3],
                          "(0.6,0.7]" =cols.r2[4],
                          "(0.7,0.8]"=cols.r2[5],
                          "(0.8,0.9]" =cols.r2[6]))
rownames(anno.row) <- rownames(ph.pd)  
pheatmap(ph.pd[ord,-1],scale = 'none',cluster_cols = F,cluster_rows = F,
         color = ,gaps_col = seq(0,ncol(ph.pd)-1,by=7),
         annotation_row = anno.row,annotation_colors = anno.row.col)

# fig4 cmp sp  ------------------------------------------------------------

# calc sp 
pd.sp <- pd %>% group_by(geno,gene,type,sti) %>% summarise(max.frac = max(normCnt.frac))
pd.sp<- pd.sp %>% group_by(gene,geno,type) %>% 
  summarise(sp = log2(first(max.frac)/last(max.frac)))
pd.sp <- ungroup(pd.sp)

# draw box 
all.equal(rownames(dat)[ord],unique(pd.sp$gene)[ord])
n.gene <- length(unique(pd.sp$gene))

for(g in c('wt','mt')){
  setEPS()
  postscript(file=paste0(fig_folder,'subfig4.spcmp.',g,'.eps'),
             height = 8,width = 2)
  plotSp(genotype = g)
  dev.off()
}

# old plot-heatmap --------------------------------------------------------

## plot sim vs. exp 
cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)
bks <- seq(0,1.0001,by = .1)

for(n in c("mt_lps","wt_lps",'mt_tnf','wt_tnf')){
  lswitch <- ifelse(n =='mt_tnf',T,F)
  setEPS()
  postscript(paste0(fig_folder,n,".eps"),onefile = F,width = 2,height = 6)
  tmp<- unlist(strsplit(n,'_'))
  dat <- subset(pd,geno==tmp[1] & sti==tmp[2])  %>% 
    unite(key,geno,type,sti,time) %>%
    spread(key,normCnt.frac)
  rownames(dat ) <- dat[,1]; dat <- dat[,-1]
  #seps <- NULL
  seps <- cumsum(table(cut(anno.row$r2,bks.r2[-2])))+1
  pheatmap.2(pd=dat[ord,],bks=bks,cols = cols,cwid = 5,
             legend=F,gaps_col =7)
  dev.off()
}

plotLegend(cols = cols,
           bks = bks,
           fnames = paste0(fig_folder,"sim_hm_scale.eps"))

cat(rownames(ph.pd)[ord],sep = '\n')
gene.dic <- read.csv(file='../data/mRNA.cluster_old.csv',stringsAsFactors = F)
all(rownames(ph.pd) %in% gene.dic$gene)
rownames(gene.dic) <- gene.dic$gene
cat(gene.dic[rownames(ph.pd)[ord],]$gene.2,sep = '\n')
data.frame(gene.dic[rownames(ph.pd)[ord],]$gene.2,
           rownames(ph.pd)[ord])

## plot sim vs. exp : R2 
bks.r2 <- c(min(anno.row)-0.001,0,.5,.6,.7,.8,.9,1.0)

bks.r2 <- c(-0.001,.5,.6,.7,.8,.9,1.0)
cols.r2 <- colorRampPalette(c("black","yellow"))(length(bks.r2)-1)
anno.row.2 <- anno.row; anno.row.2[anno.row<0] <- 0 
seps <- cumsum(table(cut(anno.row.2$r2,bks.r2)))+1
postscript(paste0(fig_folder,"fitR2.eps"),
           onefile = F,width = 1,height = 6)
pheatmap.2(pd=anno.row.2[ord,],bks=bks.r2,cols=cols.r2,legend=F,
           cwid = 10,border_color = grey(.5))  
dev.off()

plotLegend(cols = cols.r2,bks = bks.r2,fnames = paste0(fig_folder,"fitR2_hm_scale.eps"))

# clust cate III data  ----------------------------------------------------
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 
sp.th <- 2^0.5
sp.mat <- read.csv(file='../data/mRNA.sp.peak.csv',row.names = 1,stringsAsFactors = F)
sp.pd <- data.frame(ctrl.sp= sp.mat$b1.LT.ctrl,
                 mt.sp = sp.mat$b1.LT.ko,
                 gene = rownames(sp.mat),
                 stringsAsFactors = T) %>% 
  mutate(cate=ifelse(ctrl.sp<0.5,'I',
                     ifelse(mt.sp>ctrl.sp*2^(-.5),'II','III')))


pd <- subset(maxScale.mRNA,maxScale.mRNA$gene %in%
               unique(sp.pd$gene[sp.pd$cate=="III"]))
ph.pd<-pd %>% unite(key,geno,sti,time) %>%
  spread(key,normCnt.frac)
rownames(ph.pd)  <- ph.pd$gene; ph.pd$gene<-NULL 

## kmeans
ks.res <- sapply(c(2:10,15,20,25,87),function(x){
  set.seed(20)
  ks <- kmeans(ph.pd,centers = x)
  ks$tot.withinss
  })
plot(c(2:10,15,20,25,87),ks.res)
set.seed(20)
ks <- kmeans(ph.pd,centers = 10)

## plot sim vs. exp 
cols <- colorRampPalette(brewer.pal(9,"Reds"))(10)
bks <- seq(0,1.0001,by = .1)
cl_ord <- c(4,2,7,10,3,5,8,1,6,9)
ord <- unlist(sapply(cl_ord,function(x) which(ks$cluster==x)))
cord <- c(15:28,1:14)
seps <- cumsum(unlist(sapply(cl_ord,function(x) sum(ks$cluster==x))))

anno_row <- data.frame(cluster=factor(as.numeric(factor(ks$cluster,
                                                        levels = cl_ord))))
rownames(anno_row) <- rownames(ph.pd)
write.csv(file='kclusters.csv',anno_row,quote = F)
anno_row.2 <- data.frame(nrmsd= as.factor(pd.nRMSD$nrmsd<0.13 ))
rownames(anno_row.2) <- pd.nRMSD$gene
pheatmap(ph.pd[ord,cord],scale = "none",cluster_rows = F,
         cluster_cols = F,show_rownames = T,
         color = cols,breaks = bks,gaps_row = seps,
         #fontsize_row = 4,
         gaps_col = c(7,14,21),
         annotation_row = anno_row,
         annotation_legend = T)

## average profile 
pd <- subset(maxScale.mRNA,maxScale.mRNA$gene %in%
               unique(sp.pd$gene[sp.pd$cate=="III"]))
pd <- pd %>% mutate(clust=anno_row[gene,"cluster"])


save(file = 'clustering_cateIII.Rdata',list = c('pd','pd.sum'))

# clusting III data 2 -----------------------------------------------------
env=new.env();load('./data/clustering_cateIII.Rdata',env)

# filter specific genes 
pd <- pd%>%filter(gene!=f.genes)

# update few gene name
dic.gene<- getSymbol(grep('ENSM*',unique(pd$gene),value = T))
pd$gene[grep('ENSM*',pd$gene)]<- dic.gene[pd$gene[grep('ENSM*',pd$gene)]]

# transform to heatmap format
ph.pd<-pd %>% unite(key,geno,sti,time) %>%
  spread(key,normCnt.frac)
rownames(ph.pd)  <- ph.pd$gene; ph.pd$gene<-NULL 

# plot heatmap
cl_ord <- c(1,2,4,5,7,8,9)
ord <- unlist(sapply(cl_ord,function(x) which(ph.pd$clust==x)))
cord <- c(15:28,1:14)+1
seps <- cumsum(unlist(sapply(cl_ord,function(x) sum(ph.pd$clust==x))))
bks <- seq(0,1.0001,by = .1)

pheatmap(ph.pd[ord,cord],scale = "none",cluster_rows = F,
         cluster_cols = F,show_rownames = T,
         color = cols,breaks = bks,gaps_row = seps,
         fontsize_row = 5,
         gaps_col = c(7,14,21))

save(file = 'clustering_cateIII_clean.Rdata',list = c('pd','pd.sum'))
# v1 fitting for each clusters  -----------------------------------------
v1.simData <- read.csv(file='./mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 
v1.weight <- calc_scores.old(v1.simData)
all.equal(v1.weight$gene,v1$pd$gene)
v1$pd$nrmsd_r_w <- v1.weight$nrmsd_r
v1$pd$nrmsd_m_w <- v1.weight$nrmsd_m

for (i in 1:9){
  selc.genes<- rownames(anno_row)[anno_row==i]
  selc.genes<- selc.genes[selc.genes%in%v1.simData$gene]
  #plotEgFit(eg.genes = selc.genes,wide = T)
  tmp <-(v1$pd %>% filter(gene %in% selc.genes))
  tmp[,2:10] <- signif(tmp[,2:10],3)
  write.csv(file=paste0('tmp',i,'.csv'),
                        t(tmp))
}

v1.old<-calc_scores.old(v1.simData,weigh.pt = 0)
sapply(c("Ripk2","Tnfaip2","Cebpb","Adtrp","Kcnn3"),
       function(x) v1.old%>%filter(gene==x) %>% select(nrmsd_r))
sapply(c("Trim47","Angptl4"),
       function(x) v1.old%>%filter(gene==x)%>%select(nrmsd_r))
       
# explore adjR2-Pilra  ----------------------------------------------------
fit.r2 <- read.csv(file = '../fig.4_modelv1/mRNA-Fit-avg-sp-v1d/result.csv',row.names = 1)
fit.r2.ca <- read.csv(file = '../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-funs-v6-r2.csv',row.names = 1)
head(fit.r2);head(fit.r2.ca)
fit.r2 <- fit.r2[rownames(fit.r2.ca),]
pd.r2 <- rbind(data.frame(fit.r2.ca[,-(3:5)],type="ca"),
               data.frame(fit.r2[,c("r2","gene")],type="emsa"))
pd.r2$type <- relevel(pd.r2$type,'emsa')
pd.r2.2<- pd.r2 %>% spread(type,r2) %>% 
  mutate(delta=ca-emsa) %>% arrange(desc(ca))

# selected genes 
selc.genes <- as.character((pd.r2.2 %>% filter(ca>=0.7))$gene)
pd.r2 <- pd.r2 %>% filter(gene %in% selc.genes)
# plot changes 
ggplot(pd.r2 %>% filter(gene %in% selc.genes),aes(type,r2,group=gene)) + geom_line() + geom_point() + 
  ylim(0,1) +   geom_hline(yintercept = .7,linetype=2,colour=grey(.3)) +
  geom_text(data = (pd.r2.2%>% filter(gene %in% selc.genes))[1:5,],colour='red',
            aes(x=2.2,y=ca,label=gene)) + 
  geom_text(data = (pd.r2.2%>% filter(gene %in% selc.genes))[6:9,],colour='red',
            aes(x=0.8,y=emsa,label=gene)) #,position = position_dodge(width = .5)

# load best fit data emsa
v1.simData <- read.csv(file='./mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 

p <- plotEgFit(eg.genes = selc.genes,wide = T)


# load caRNA fit data 
v2.simData <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6a/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v2.simData) <- sub("mRNA","normCnt.frac",names(v2.simData)) 
v2.simData$sti <- toupper(v2.simData$sti)
head(p$data);head(v2.simData)
p.2 <- p %+% rbind(p$data,data.frame(v2.simData%>%filter(gene %in%selc.genes),type='caSim.'))


# load caRNA exp data 
ev <- new.env()
load("../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-pre-v6.Rdata",ev)
pd.ca.exp <- ev$caRNA.frac.avg[selc.genes,]
pd.ca.exp <- pd.ca.exp%>% mutate(gene=rownames(pd.ca.exp))%>%
  gather(key = sample,value = normCnt.frac,1:24)%>%
  mutate(time=sapply(sample,function(x) as.numeric(unlist(strsplit(x,split = "_"))[3])),
         sti = sapply(sample,function(x) unlist(strsplit(x,split = "_"))[2]),
         geno= sapply(sample,function(x) unlist(strsplit(x,split = "_"))[1]),
         type="caExp.")
geno.dic <- c(I="wt",D="mt");sti.dic <- c(L='LPS',T='TNF')
pd.ca.exp$sti <-sti.dic[pd.ca.exp$sti]
pd.ca.exp$geno <-geno.dic[pd.ca.exp$geno]
pd.ca.exp$gene <- factor(pd.ca.exp$gene,levels = selc.genes)

p.2 + geom_point(data = pd.ca.exp,
                 aes(time,normCnt.frac,shape=geno,colour=sti)) +
  geom_line(data = pd.ca.exp,
            aes(linetype=geno))

# r2 plot 
pd.r2$gene <- droplevels(pd.r2$gene);
pd.r2$gene <- factor(pd.r2$gene,levels = levels(p.2$data$gene))
p.r2 <- ggplot(pd.r2,aes(type,r2)) + geom_bar(stat = 'identity') + facet_grid(.~gene)
p.r2 + theme_bw() + theme(axis.title = element_blank(),
                          legend.position = "none")+ scale_y_continuous(position = 'right')


# overlay with half-life 
load("./mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata",ev)
pd.hfs <- round(log(2)/ev$genes.kdeg)
pd.hfs[selc.genes]



# explore Ripk2 -----------------------------------------------------------
g <- 'Ripk2'
g.v1 <- v1.simData%>% filter(gene ==g)
g.v2 <- v2.simData %>% filter(gene==g)
g.exp <- maxScale.mRNA %>% filter(gene==g)
g.v1.exp <- g.v1 %>% filter(time %in% g.exp$time)
g.v2.exp <- g.v2 %>% filter(time %in% g.exp$time)
if(T){
  par(mfrow=c(2,1))
  plot(g.exp$normCnt.frac,ylab='Expression',xlab='data points')
  lines(g.v1.exp$normCnt.frac)
  lines(g.v2.exp$normCnt.frac,col=2)
  g.v1.res <- g.v1.exp$normCnt.frac - g.exp$normCnt.frac
  g.v2.res <- g.v2.exp$normCnt.frac - g.exp$normCnt.frac
  plot(g.v1.res,pch=16,ylim=range(c(g.v1.res,g.v2.res)),
       ylab='residuals',xlab='data.points')
  points(g.v2.res,pch=16,col=2)
  abline(h=0)
  par(mfrow=c(1,1))
  
}

g.v1.rmsd <- sqrt(sum(g.v1.res^2)/length(g.v1.res))
g.v1.r2 <- 1- var(g.v1.res)/var(g.v1.exp$normCnt.frac)
ss.tot <- sum((g.exp$normCnt.frac-mean(g.exp$normCnt.frac))^2)
g.v1.r2.new <- 1- sum(g.v1.res^2)/ss.tot 
g.v2.rmsd <- sqrt(sum(g.v2.res^2)/length(g.v1.res))
g.v2.r2 <- 1- var(g.v2.res)/var(g.v1.exp$normCnt.frac)
g.v2.r2.new <- 1- sum(g.v2.res^2)/ss.tot

data.frame(g.v1.rmsd,g.v2.rmsd)
data.frame(g.v1.r2,g.v2.r2,g.v1.r2.new,g.v2.r2.new)

# half-life explore -------------------------------------------------------

normCnt <- read.csv(file="../../half-life/Supriya/ActDBrowser/allNormCount.csv",
                    row.names = 1)
gene.dic  <- read.csv(file="../../half-life/Supriya/ActDBrowser/gene.dic.csv",
                      stringsAsFactors = F,row.names = 1)
idx <- rownames(gene.dic)[which(gene.dic=="Ccl1")]
pd.normCnt <- normCnt[idx,]
pd.normCnt <- data.frame(time = c(0,.5,1,2,4,6),
                         val= as.numeric(pd.normCnt[grep("b2_3.0_LPS",names(pd.normCnt))]))
pd.normCnt$val <- round(pd.normCnt$val)
glm.fit <- glm(val~time, data = pd.normCnt,family = "poisson"(link = "log"))
plot(pd.normCnt)
lines(0:6,exp(predict(glm.fit,data.frame(time=0:6))))
1-glm.fit$deviance/glm.fit$null.deviance # 0.821542 

pd.normCnt <- normCnt[idx,]
pd.normCnt <- data.frame(time = c(0,.5,1,2,4,6),
                         val= as.numeric(pd.normCnt[grep("b2_3.0_LPS",names(pd.normCnt))]))
pd.normCnt <- data.frame(time = c(0,.5,1,2,4,6),
                         val= as.numeric(pd.normCnt[grep("(b1_3.0_LPS)|(b1_0.0_U_0.0)",names(pd.normCnt))]))

pd.normCnt$val <- log2(pd.normCnt$val+1)
lm.fit.1 <- lm(val~time,pd.normCnt)
summary(lm.fit.1)
lm.fit.2 <- lm(val~time,pd.normCnt[4:6,])
summary(lm.fit.2) 
plot(pd.normCnt); 
abline(lm.fit.1); abline(lm.fit.2,col=2)
car::Anova(lm.fit.1,lm.fit.2)

# cmp new/old R2 vs. RMSD -------------------------------------------------
v1.simData <- read.csv(file='./mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 
cmp_R2 <- function(v1.simData){
  v1.simData.exp <- v1.simData%>%filter(time %in% maxScale.mRNA$time)
  maxScale.mRNA <- maxScale.mRNA%>% filter(gene %in% v1.simData$gene)
  maxScale.mRNA <- maxScale.mRNA[,colnames(v1.simData.exp)]
  maxScale.mRNA$sti <- as.character(maxScale.mRNA$sti)
  maxScale.mRNA$geno <- as.character(maxScale.mRNA$geno)
  all.equal(v1.simData.exp[,-2],maxScale.mRNA[,-2])
  maxScale.mRNA$res <- v1.simData.exp$normCnt.frac - maxScale.mRNA$normCnt.frac
  pd<- maxScale.mRNA%>% group_by(gene) %>% 
    summarise(rmsd = sqrt(sum(res^2)/length(res)),
              mean_data= mean(normCnt.frac),
              range_data=max(normCnt.frac)-min(normCnt.frac),
              r2.old = 1-var(res)/var(normCnt.frac),
              r2.new = 1- sum(res^2)/sum((normCnt.frac - mean(normCnt.frac))^2)) %>%
    mutate(nrmsd_m=rmsd/mean_data,
           nrmsd_r=rmsd/range_data)
  pd %>% gather(r2.type,r2,2:3)
  print(signif(quantile(pd$r2.new-pd$r2.old,c(0,.25,.5,.75,1)),2))
  p<-ggplot(pd,aes(rmsd,r2.old)) + geom_point()+ 
    geom_point(aes(rmsd,r2.new),colour=2,shape=2)+
    ylab("R2") + ggtitle("R2 new (red) vs. R2 old (black)")
  list(p=p,pd=pd)
}
v1<-cmp_R2(v1.simData)
v2 <- cmp_R2(v2.simData)
write.csv(file = './tmp.csv',v2$pd%>%filter(r2.old>0.7) %>% arrange(desc(r2.old)),
          quote = F)

# adjust R2 by weighing peaks  --------------------------------------------
# replicate peak rows 
cmp_adjR2 <- function(v1.simData,rep.num=5){
  v1.simData.exp <- v1.simData%>%filter(time %in% maxScale.mRNA$time)
  maxScale.mRNA <- maxScale.mRNA%>% filter(gene %in% v1.simData$gene)
  maxScale.mRNA <- maxScale.mRNA[,colnames(v1.simData.exp)]
  maxScale.mRNA$sti <- as.character(maxScale.mRNA$sti)
  maxScale.mRNA$geno <- as.character(maxScale.mRNA$geno)
  all.equal(v1.simData.exp[,-2],maxScale.mRNA[,-2])
  maxScale.mRNA$res <- v1.simData.exp$normCnt.frac - maxScale.mRNA$normCnt.frac
  
  # adding peak replicates 
  maxScale.mRNA <-maxScale.mRNA %>% 
    group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%
    dplyr::select(1:6) %>% 
    dplyr::slice(rep(1:n(),each=rep.num)) %>% #5 
    rbind(maxScale.mRNA)
  
  
  pd<- maxScale.mRNA%>% group_by(gene) %>% 
    summarise(rmsd = sqrt(sum(res^2)/length(res)),
              mean_data= mean(normCnt.frac),
              range_data=max(normCnt.frac)-min(normCnt.frac),
              r2.old = 1-var(res)/var(normCnt.frac),
              r2.new = 1- sum(res^2)/sum((normCnt.frac - mean(normCnt.frac))^2))%>%
    mutate(nrmsd_m = rmsd/mean_data,
           nrmsd_r = rmsd/range_data)
  pd %>% gather(r2.type,r2,2:3)
  print(signif(quantile(pd$r2.new-pd$r2.old,c(0,.25,.5,.75,1)),2))
  p<-ggplot(pd,aes(rmsd,r2.old)) + geom_point()+ 
    geom_point(aes(rmsd,r2.new),colour=2,shape=2)+
    ylab("R2") + ggtitle("R2 new (red) vs. R2 old (black)")
  list(p=p,pd=pd)
}

v1.Radj <- cmp_adjR2(v1.simData)

v2.Radj <- cmp_adjR2(v2.simData)

v1.Radj$pd%>%filter(gene=="Ripk2" )

pd.cmp <- rbind(data.frame(v1$pd,model="v1",weigh=F),
                data.frame(v2$pd,model="v2",weigh=F),
                data.frame(v1.Radj$pd,model="v1",weigh=T),
                data.frame(v2.Radj$pd,model="v2",weigh=T))

#rmsd v1 vs v2 
ggplot(pd.cmp %>% filter(!weigh)%>%
  dplyr::select(gene,rmsd,model,weigh) %>%
  spread(model,rmsd), aes(v1,v2)) + geom_point() + 
  geom_abline(slope = 1,intercept = 0)


if(T){
  x <- (pd.cmp%>%filter(!weigh & model=="v1"))$rmsd
  y <- (pd.cmp%>%filter(weigh & model=="v1"))$rmsd
  plot(x,y,xlab="rmsd before weigh",ylab="rmsd after weigh")
  abline(a=0,b=1)
  sapply(c(0,.25,.5,.75,1),function(n) {
    m <- quantile(x,n)
    abline(v=m,lty=2)
    text(m,.3,labels = paste(names(m),signif(m,2)),srt=90)})
  sapply(c(0,.25,.5,.75,1),function(n) {
    m <- quantile(y,n)
    abline(h=m,lty=2)
    text(.3,m,labels = paste(names(m),signif(m,2)))})
}

id <- identify(x,y)
data.frame((pd.cmp%>%filter(!weigh & model=="v1"))[id,],
           rmsd.w = (pd.cmp%>%filter(weigh & model=="v1"))[id,"rmsd"])
write.csv(file='tmp.csv',data.frame((pd.cmp%>%filter(!weigh & model=="v1"))[id,c(1,2)],
           rmsd.w = (pd.cmp%>%filter(weigh & model=="v1"))[id,"rmsd"]))
selc.genes <- (pd.cmp%>%filter(!weigh & model=="v1"))$gene[id]
plotEgFit(selc.genes,wide = T)

x<-(pd.cmp%>%filter(!weigh & model=="v1"))$r2.new
y<- (pd.cmp%>%filter(weigh & model=="v1"))$r2.new
plot(x, y,
     xlab="r2.new before weigh",ylab="r2.new after weigh",
     xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1);abline(v=.7,lty=2);abline(h = .7,lty=2)

table(x>.7,y>.7)
table((pd.cmp%>%filter(!weigh & model=="v1"))$r2.old>=.7,x>=.7)

# weigh gained genes 
t.genes <- (pd.cmp%>%filter(!weigh & model=="v1"))$gene[y>=0.7 & x <0.7 ]

# plot these genes results
pd.cmp.v1 <-pd.cmp%>% filter(gene %in% t.genes & model=="v1") %>%
  dplyr::select(gene,r2.new,weigh) %>%
  mutate(weigh=paste0("w_",weigh))%>%
  spread(weigh,r2.new) %>% 
  arrange(w_TRUE)
pd.cmp.v1.rmsd <-pd.cmp%>% filter(gene %in% t.genes & model=="v1") %>%
  dplyr::select(gene,rmsd,weigh) %>%
  mutate(weigh=paste0("w_",weigh))%>%
  spread(weigh,rmsd) 

tmp <- pd.cmp.v1.rmsd[match(pd.cmp.v1$gene,pd.cmp.v1.rmsd$gene),]

plotEgFit(eg.genes = pd.cmp.v1$gene,wide = T)  
write.csv(file="tmp.csv",t(tmp),quote = F)



# NRMSD thresholding  ------------------------------------------------------
# plot as from bad to good using nrmsd_m value 

v1.Radj <- cmp_adjR2.2(v1.simData)

selc.genes <- (v1.Radj$pd %>% arrange(nrmsd_m))
selc.genes <- (v1.Radj$pd %>% arrange(desc(r2.new)))
ords<- c(seq(1,88,by=9),89) 

pdf(file = 'adjR2_check_v1.pdf',height = 3,width = 9)
sapply(2:length(ords),
       function(i) {
         print(plotEgFit(eg.genes = selc.genes$gene[ords[i-1]:(ords[i]-1)],
                                 wide = T))
         })
dev.off()
for(i in 2:length(ords))
  print( selc.genes$r2.new[ords[i-1]:(ords[i]-1)])



# fake weigh both pt and pv ----------------------------------------------------


calc_scores.old <- function(v1.simData,weight.pv=5,
                        weigh.pt=2){ # weigth on both peak and peak location
  
  # norm residual 
  v1.simData.exp <- v1.simData%>%filter(time %in% maxScale.mRNA$time)
  maxScale.mRNA <- maxScale.mRNA%>% filter(gene %in% v1.simData$gene)
  maxScale.mRNA <- maxScale.mRNA[,colnames(v1.simData.exp)]
  maxScale.mRNA$sti <- as.character(maxScale.mRNA$sti)
  maxScale.mRNA$geno <- as.character(maxScale.mRNA$geno)
  norm.res <- maxScale.mRNA
  norm.res$res <- v1.simData.exp$normCnt.frac - maxScale.mRNA$normCnt.frac
  
  # peak value residual
  peak.sim <- v1.simData%>% group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%dplyr::select(1:5)
  
  peak.exp <-   maxScale.mRNA%>% group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%dplyr::select(1:5) 
  
  all.equal(with(peak.sim,paste(sti,geno,gene,sep="-")),
            with(peak.exp,paste(sti,geno,gene,sep="-")))
  
  peak.res <- peak.sim; peak.res$res <- peak.sim$normCnt.frac - peak.exp$normCnt.frac
  peak.res <- peak.res%>% dplyr::slice(rep(1:n(),each=weight.pv))#weight it 
  
  # add peak time res
  peak.sim <- peak.sim %>% mutate(normCnt.frac = time/8 )
  peak.exp <- peak.exp %>% mutate(normCnt.frac = time/8 ) #scale time from 0-1 
  peak.res<- peak.res%>% rbind(peak.sim %>% 
                      mutate(res=peak.sim$normCnt.frac-peak.exp$normCnt.frac)%>%
                      dplyr::slice(rep(1:n(),each=weigh.pt)))
  
  # merge peak and norm res
  all.res <- rbind(peak.res,norm.res)
  
  
  pd<- all.res%>% group_by(gene) %>% 
    summarise(rmsd = sqrt(sum(res^2)/length(res)),
              mean_data= mean(normCnt.frac),
              range_data=max(normCnt.frac)-min(normCnt.frac),
              r2.new = 1- sum(res^2)/sum((normCnt.frac - mean(normCnt.frac))^2))%>%
    mutate(nrmsd_m = rmsd/mean_data,
           nrmsd_r = rmsd/range_data)
}

# calc new nrmsd and plot it 
selc.genes <- calc_scores(v1.simData) %>% arrange(nrmsd_m)
selc.genes <- selc.genes %>% arrange(desc(r2.new))
table(selc.genes$nrmsd_m <0.3)
table(selc.genes$r2.new >0.7)
ords<- c(seq(1,88,by=9),89) 

pdf(file = 'weighted_pv_5_pt3_ordR.pdf',height = 3,width = 9)
sapply(2:length(ords),
       function(i) print(plotEgFit(eg.genes = selc.genes$gene[ords[i-1]:(ords[i]-1)],
                                   wide = T)))
dev.off()

# Real weight  ------------------------------------------------------------
calc_score <- function(v1.simData,weight.pv=2,
                        weigh.pt=1){ # weigth on both peak and peak location
  
  # norm residual 
  v1.simData.exp <- v1.simData%>%filter(time %in% maxScale.mRNA$time)
  maxScale.mRNA <- maxScale.mRNA%>% filter(gene %in% v1.simData$gene)
  maxScale.mRNA <- maxScale.mRNA[,colnames(v1.simData.exp)]
  maxScale.mRNA$sti <- as.character(maxScale.mRNA$sti)
  maxScale.mRNA$geno <- as.character(maxScale.mRNA$geno)
  norm.res <- maxScale.mRNA
  norm.res$res <- v1.simData.exp$normCnt.frac - maxScale.mRNA$normCnt.frac
  
  # peak value residual
  peak.sim <- v1.simData%>% group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%dplyr::select(1:5)
  
  peak.exp <-   maxScale.mRNA%>% group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%dplyr::select(1:5) 
  
  all.equal(with(peak.sim,paste(sti,geno,gene,sep="-")),
            with(peak.exp,paste(sti,geno,gene,sep="-")))
  
  peak.res <- peak.sim; peak.res$res <- peak.sim$normCnt.frac - peak.exp$normCnt.frac
  peak.res <- peak.res%>% mutate(res=weight.pv*res)#weight it 
  
  # add peak time res
  peak.sim.t <- peak.sim %>% mutate(normCnt.frac = time/8 )
  peak.exp.t <- peak.exp %>% mutate(normCnt.frac = time/8 ) #scale time from 0-1 
  peak.res<- peak.res%>% rbind(peak.sim %>% 
                                 mutate(res=(peak.sim.t$normCnt.frac-peak.exp.t$normCnt.frac)*weigh.pt))
  
  # merge peak and norm res
  all.res <- rbind(peak.res,norm.res)
  
  
  pd<- all.res%>% group_by(gene) %>% 
    summarise(rmsd = sqrt(sum(res^2)/length(res)),
              mean_data= mean(normCnt.frac),
              range_data=max(normCnt.frac)-min(normCnt.frac),
              r2.new = 1- sum(res^2)/sum((normCnt.frac - mean(normCnt.frac))^2))%>%
    mutate(nrmsd_m = rmsd/mean_data,
           nrmsd_r = rmsd/range_data)
}

v1.new <- calc_score(v1.simData)
v1.new%>%filter(gene=="Acpp")
v1.old <- calc_scores.old(v1.simData)
v1.old%>%filter(gene=="Acpp")

Trim47.rank <- sapply(1:20,function(x)
  {v1.new <- calc_score(v1.simData,weight.pv = x) %>% arrange(nrmsd_m)
  which(v1.new$gene=="Ripk2")} )

v1.new <- calc_score(v1.simData,weight.pv = 20) %>% arrange(nrmsd_m)
data.frame(1:20,Trim47.rank)



pdf(file = 'weighted_pv_20_ord.pdf',height = 3,width = 9)
for (i in 2:length(ords))
  print(plotEgFit(eg.genes = v1.new$gene[ords[i-1]:(ords[i]-1)],
                                   wide = T));
dev.off()

table(v1.new$r2.new>.7)
table(v1.new$nrmsd_m<.3)
v1.new %>% filter(gene=="Trim47")

# cmp r2 and caRNA sp -----------------------------------------------------
sp.caRNA <- read.csv(file='../fig.5_caRNA//data/caRNA-max-avg-sp.csv',header = T,
                     row.names = 1,stringsAsFactors = F)
sp.mRNA <- read.csv(file="../data/mRNA.sp.peak.csv",header = T,
                    row.names = 1,stringsAsFactors = F)

fit.r2 <- read.csv(file = '../fig.4_modelv1/mRNA-Fit-avg-sp-v1d/result.csv',row.names = 1)
all(rownames(fit.r2)%in% rownames(sp.caRNA))
all(rownames(fit.r2)%in% rownames(sp.mRNA))

pd.sp <- data.frame(sp.caRNA[rownames(fit.r2),],
                    r2=fit.r2$r2,
                    sp.mRNA[rownames(fit.r2),c("b1.LT.ctrl","b1.LT.ko")])
require(ggplot2)
bks.r2 <- c(min(pd.sp$r2)-0.001,0,.5,.6,.7,.8,.9,1.0)
cols.r2 <- colorRampPalette(c("black","yellow"))(length(bks.r2))
sp.th <- 2^0.5
p<- ggplot(data.frame(ctrl.sp=c(0.5,0.5,8,8,0.5,0.5,8,8,-3,-3,0.5,0.5),
                      mt.sp=c(0.5/sp.th,8,8,8/sp.th, 
                              -3,0.5/sp.th,8/sp.th,-3,
                              -3,8,8,-3),
                      gp=factor(rep(c(2,3,1),each=4))),
           aes(ctrl.sp,mt.sp))+   coord_cartesian(xlim=range(pd.sp[,1:2]),ylim=range(pd.sp[,1:2])) +
  geom_polygon(aes(fill=gp)) + scale_fill_manual(values = alpha(c("blue"), c(.05,.1,.2)))+
  theme_bw()+
  geom_abline(slope = 1,intercept = 0,colour=grey(.8)) +
  geom_hline(yintercept = c(-.5,.5),colour=grey(.8),linetype=2)+
  geom_vline(xintercept = c(-.5,.5),colour=grey(.8),linetype=2) +
  geom_point(data = pd.sp,aes(LT.ctrl,LT.mt,colour=cut(r2,bks.r2)),size=2)+
  geom_abline(slope = 1/sp.th,intercept = 0,colour="blue",linetype=2) +
  scale_colour_manual(values = cols.r2)

ggsave(file=paste0(fig_folder,'subfig4_caRNA.eps'),p+
         theme(legend.position = 'none',axis.title.x = element_blank(),axis.title.y = element_blank(),
               axis.text.x = element_blank(),axis.text.y = element_blank()),
       width = 2.5,height = 2.5,device=cairo_ps)


pd.sp <- pd.sp %>% mutate(cate=ifelse(LT.ctrl<0.5,'I',
                                      ifelse(LT.mt>LT.ctrl*2^(-.5),'II','III')))
pd.sp <- pd.sp %>% mutate(isGoodFit=ifelse(r2>=0.7,T,F),gene=rownames(fit.r2))

(pd.sp %>% unite(col = cb,cate,isGoodFit)%>% dplyr::filter(cb=='I_FALSE'))$gene
(pd.sp %>% unite(col = cb,cate,isGoodFit)%>% dplyr::filter(cb=='I_TRUE'))$gene

# delta 

pd.sp <- pd.sp %>% mutate(wt.delta= b1.LT.ctrl - LT.ctrl,
                          mt.delta=b1.LT.ko - LT.mt)

require(ggforce)
ggplot(pd.sp,aes(wt.delta,mt.delta)) + geom_point(aes(colour=cut(r2,bks.r2)),size=2) +
  scale_colour_manual(values = cols.r2)+theme_bw()+theme(legend.position = "none") + 
  geom_vline(xintercept = 0.5,linetype=2)
ggsave(filename = paste0(fig_folder,'subfig.delta_sp.eps'),width = 2.5,height = 2.5)

pd.2 <- table(cut((pd.sp %>% filter(wt.delta<0.5))$r2,breaks = c(-1,.7,1)))
pd.2 <- rbind(pd.2,table(cut((pd.sp %>% filter(wt.delta>=0.5))$r2,breaks = c(-1,.7,1))))
pd.2 <- as.data.frame(pd.2)
pd.2 <- tidyr::gather(pd.2,key = r2,value = num, c(1,2)) %>% mutate(consistent=rep(c(T,F),2))
pd.2$consistent <- factor(pd.2$consistent,levels = c(T,F))

ggplot(pd.2,aes(consistent,num,fill=r2)) + 
  geom_bar(stat='identity') + scale_fill_manual(values = cols.r2[c(1,length(cols.r2))])+ 
  theme_bw()+theme(legend.position = "bottom")
ggsave(filename = paste0(fig_folder,'subfig.r2vSpdelta.eps'),width = 2,height = 3)

# explore caFit -----------------------------------------------------------
fit.r2 <- read.csv(file = '../fig.4_modelv1/mRNA-Fit-avg-sp-v1d/result.csv',row.names = 1)
fit.r2.ca <- read.csv(file = '../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-funs-v6-r2.csv',row.names = 1)
head(fit.r2);head(fit.r2.ca)
fit.r2 <- fit.r2[rownames(fit.r2.ca),]
pd.r2 <- rbind(data.frame(fit.r2.ca[,-(3:5)],type="ca"),
               data.frame(fit.r2[,c("r2","gene")],type="emsa"))
pd.r2$type <- relevel(pd.r2$type,'emsa')
pd.r2.2<- pd.r2 %>% spread(type,r2) %>% mutate(delta=ca-emsa) %>% 
  arrange(desc(ca))

# selected genes 
selc.genes <- as.character((pd.r2.2 %>% filter(ca>=0.7))$gene)
pd.r2 <- pd.r2 %>% filter(gene %in% selc.genes)
# plot changes 
ggplot(pd.r2 %>% filter(gene %in% selc.genes),aes(type,r2,group=gene)) + geom_line() + geom_point() + 
  ylim(0,1) +   geom_hline(yintercept = .7,linetype=2,colour=grey(.3)) +
  geom_text(data = (pd.r2.2%>% filter(gene %in% selc.genes))[1:5,],colour='red',
            aes(x=2.2,y=ca,label=gene)) + 
  geom_text(data = (pd.r2.2%>% filter(gene %in% selc.genes))[6:9,],colour='red',
            aes(x=0.8,y=emsa,label=gene)) #,position = position_dodge(width = .5)
ggsave(filename = paste0(fig_folder,'r2_cmp_sub.pdf'),width = 4,height = 8)


# load best fit data emsa
v1.simData <- read.csv(file='./mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 

p <- plotEgFit(eg.genes = selc.genes,wide = T)


# load caRNA fit data 
v2.simData <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6a/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v2.simData) <- sub("mRNA","normCnt.frac",names(v2.simData)) 
v2.simData$sti <- toupper(v2.simData$sti)
head(p$data);head(v2.simData)
p.2 <- p %+% rbind(p$data,data.frame(v2.simData%>%filter(gene %in%selc.genes),type='caSim.'))

#p + theme_bw() + theme(axis.title = element_blank(),axis.text = element_blank(),
#                       strip.text = element_blank(),legend.position = "none")
ggsave(p.2,filename = paste0(fig_folder,'v2_improve.png'),height = 4,width = 9)


# load caRNA exp data 
ev <- new.env()
load("../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-pre-v6.Rdata",ev)
pd.ca.exp <- ev$caRNA.frac.avg[selc.genes,]
pd.ca.exp <- pd.ca.exp%>% mutate(gene=rownames(pd.ca.exp))%>%
  gather(key = sample,value = normCnt.frac,1:24)%>%
  mutate(time=sapply(sample,function(x) as.numeric(unlist(strsplit(x,split = "_"))[3])),
         sti = sapply(sample,function(x) unlist(strsplit(x,split = "_"))[2]),
         geno= sapply(sample,function(x) unlist(strsplit(x,split = "_"))[1]),
         type="caExp.")
geno.dic <- c(I="wt",D="mt");sti.dic <- c(L='LPS',T='TNF')
pd.ca.exp$sti <-sti.dic[pd.ca.exp$sti]
pd.ca.exp$geno <-geno.dic[pd.ca.exp$geno]
pd.ca.exp$gene <- factor(pd.ca.exp$gene,levels = selc.genes)

p.2 + geom_point(data = pd.ca.exp,
                 aes(time,normCnt.frac,shape=geno,colour=sti)) +
  geom_line(data = pd.ca.exp,
            aes(linetype=geno))

# r2 plot 
pd.r2$gene <- droplevels(pd.r2$gene);
pd.r2$gene <- factor(pd.r2$gene,levels = levels(p.2$data$gene))
p.r2 <- ggplot(pd.r2,aes(type,r2)) + geom_bar(stat = 'identity') + facet_grid(.~gene)
p.r2 + theme_bw() + theme(axis.title = element_blank(),
                          legend.position = "none")+ scale_y_continuous(position = 'right')



# overlay with half-life 

load("./mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata",ev)
pd.hfs <- round(log(2)/ev$genes.kdeg)
pd.hfs[selc.genes]




# explot unfit - remaining ------------------------------------------------
fit.r2 <- read.csv(file = '../fig.4_modelv1/mRNA-Fit-avg-sp-v1d/result.csv',row.names = 1)
fit.r2.ca <- read.csv(file = '../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-funs-v6-r2.csv',row.names = 1)
head(fit.r2);head(fit.r2.ca)
fit.r2 <- fit.r2[rownames(fit.r2.ca),]
pd.r2 <- rbind(data.frame(fit.r2.ca[,-(3:5)],type="ca"),
               data.frame(fit.r2[,c("r2","gene")],type="emsa"))
pd.r2$type <- relevel(pd.r2$type,'emsa')
pd.r2.2<- pd.r2 %>% spread(type,r2) %>% mutate(delta=ca-emsa) %>% arrange(desc(emsa))

# selected genes 
selc.genes <- as.character((pd.r2.2 %>% filter(ca<0.7))$gene)
pd.r2 <- pd.r2 %>% filter(gene %in% selc.genes)

# load best fit data emsa
v1.simData <- read.csv(file='./mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 
p <- plotEgFit(eg.genes = selc.genes,wide = T)



# load caRNA fit data 
v2.simData <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6a/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v2.simData) <- sub("mRNA","normCnt.frac",names(v2.simData)) 
v2.simData$sti <- toupper(v2.simData$sti)
head(p$data);head(v2.simData)
p.2 <- p %+% rbind(p$data,data.frame(v2.simData%>%filter(gene %in%selc.genes),type='caSim.'))


# load the caRNA data
ev <- new.env()
load("../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-pre-v6.Rdata",ev)
pd.ca.exp <- ev$caRNA.frac.avg[selc.genes,]
pd.ca.exp <- pd.ca.exp%>% mutate(gene=rownames(pd.ca.exp))%>%
  gather(key = sample,value = normCnt.frac,1:24)%>%
  mutate(time=sapply(sample,function(x) as.numeric(unlist(strsplit(x,split = "_"))[3])),
         sti = sapply(sample,function(x) unlist(strsplit(x,split = "_"))[2]),
         geno= sapply(sample,function(x) unlist(strsplit(x,split = "_"))[1]),
         type="caExp.")
geno.dic <- c(I="wt",D="mt");sti.dic <- c(L='LPS',T='TNF')
pd.ca.exp$sti <-sti.dic[pd.ca.exp$sti]
pd.ca.exp$geno <-geno.dic[pd.ca.exp$geno]
pd.ca.exp$gene <- factor(pd.ca.exp$gene,levels = selc.genes)

p.2 + geom_point(data = pd.ca.exp,
                 aes(time,normCnt.frac,shape=geno,colour=sti)) +
  geom_line(data = pd.ca.exp,
            aes(linetype=geno))+ theme_bw() + theme(axis.title = element_blank(),axis.text = element_blank(),
                                                    legend.position = "none")


#ggsave(p.2,filename = paste0(fig_folder,'v2_improve.png'),height = 4,width = 9)

# r2 plot 
pd.r2$gene <- droplevels(pd.r2$gene);
pd.r2$gene <- factor(pd.r2$gene,levels = levels(p.2$data$gene))
p.r2 <- ggplot(pd.r2,aes(type,r2)) + geom_bar(stat = 'identity') + facet_grid(.~gene)
p.r2 + theme_bw() + theme(axis.title = element_blank(),
                          legend.position = "none")+ scale_y_continuous(position = 'right')

# overlay with half-life 
ev <- new.env()
load("./mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata",ev)
pd.hfs <- round(log(2)/ev$genes.kdeg)
pd.hfs[selc.genes]

# overlay with changing half-life 
ev <- new.env()
load("../fig.5_caRNA/mRNA-Fit-avg-sp-v2a/mRNA-Fit-avg-sp-v2-pre.Rdata",ev)
pd.hfs <- cbind(ev$genes.kdeg.lps,ev$genes.kdeg.tnf)
pd.hfs <- data.frame(pd.hfs[selc.genes,],gene=selc.genes)
pd.hfs <- pd.hfs %>% gather(sample,hfs,1:10) %>% mutate(hfs=round(log(2)/hfs))

pd.hfs$time  <- sapply(pd.hfs$sample,function(x) 
  as.numeric(unlist(strsplit(x,split="_"))[2]))
pd.hfs$sti  <- sapply(pd.hfs$sample,function(x) unlist(strsplit(x,split="_"))[1])
pd.hfs$time[is.na(pd.hfs$time)] <- 0 
pd.hfs$sti[pd.hfs$sti=="U"] <- "LPS"; pd.hfs$sti[pd.hfs$sti=="U.1"] <- "TNF"
pd.hfs$hfs[is.na(pd.hfs$hfs)] <- min(pd.hfs$hfs,na.rm = T)
pd.hfs$gene <- factor(pd.hfs$gene,levels = selc.genes)
p.hf <- ggplot(pd.hfs,aes(time,hfs,colour=sti,group=sti)) + 
  geom_line()

p.hf + theme_bw() + facet_grid(.~gene)+
  theme(axis.title = element_blank(),
        legend.position = "none",axis.text.x = element_blank())+ 
  scale_y_continuous(position = 'right') + 
  scale_colour_manual(values = col.map)

## load v2a sim
v1a.simData <- read.csv(file='../fig.5_caRNA/mRNA-Fit-avg-sp-v2a/bestFit.tc.csv',
                        stringsAsFactors = F)
names(v1a.simData) <- sub("mRNA","normCnt.frac",names(v1a.simData)) 
v1a.simData$sti <- toupper(v1a.simData$sti)

p.3 <- p.2 %+% rbind(p.2$data,data.frame(v1a.simData%>%filter(gene %in%selc.genes),type='V1aSim.'))

v2a.simData <- read.csv(file='../fig.4_modelv1/caRNA-Fit-avg-sp-v7/bestFit.tc.csv',
                        stringsAsFactors = F)
names(v2a.simData) <- sub("mRNA","normCnt.frac",names(v2a.simData)) 
v2a.simData$sti <- toupper(v2a.simData$sti)
p.4 <- p.3 %+% rbind(p.3$data,data.frame(v2a.simData%>%filter(gene %in%selc.genes),type='V2aSim.'))


## load 
# nfkbie  -----------------------------------------------------------------
tnf3hr.hf.data <- data.frame( time=c(0,.5,1,3,6),
                              val=c(11.27,11.1,10.92,9.62,8.15))                          
tnf3hr.hf.data <- data.frame( time=c(0,.5,1,3,6),
                              val=c(6.37,7.93,7.78,6.84,4.72))                          
tnf3hr.hf.data <- data.frame( time=c(0,.5,1,3,6),
                              val=c(10.73,11.06,10.78,8.67,7.25))                          

tnf3hr.hf.fit <- lm(val~time,tnf3hr.hf.data)
par(mfrow=c(2,2))
plot(tnf3hr.hf.fit.2)
par(mfrow=c(1,1))
tnf3hr.hf.fit.2 <- lm(val~time,tnf3hr.hf.data[1:3,])
plot(tnf3hr.hf.data)
abline(tnf3hr.hf.fit)
abline(tnf3hr.hf.fit.2,col=2)
summary(tnf3hr.hf.fit); summary(tnf3hr.hf.fit.2)
car::Anova(tnf3hr.hf.fit,tnf3hr.hf.fit.2)

tnf3hr.hf.data$val <- round(2^tnf3hr.hf.data$val -1) 
glm.fit<- glm(val ~ time, data=tnf3hr.hf.data,  family="poisson"(link="log"))
tnf3hr.hf.data$val <- log(tnf3hr.hf.data$val)
glm.fit.2<- lm(val ~ time, data=tnf3hr.hf.data)

plot(tnf3hr.hf.data)

lines(0:6,exp(predict(glm.fit,newdata = data.frame(time=0:6))))
#http://stats.stackexchange.com/questions/160253/compare-two-lm-where-one-is-calculated-on-a-subset-of-the-data
tnf3hr.hf.data.2 <- rbind(data.frame(tnf3hr.hf.data,group="A"),
                          data.frame(tnf3hr.hf.data[1:3,],group="B"))

lmtest <- lm(val ~ time*group,data = tnf3hr.hf.data.2)
summary(lmtest)
log(2)/(-tnf3hr.hf.fit$coefficients[2])
summary(tnf3hr.hf.fit)
summary(tnf3hr.hf.fit.2)
log(2)/(-tnf3hr.hf.fit.2$coefficients[2])*60

anova(c(tnf3hr.hf.fit,tnf3hr.hf.fit.2))


# inconsistency found between mannual computing and the browser 


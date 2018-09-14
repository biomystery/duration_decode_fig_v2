source(file="../auxilary_functions.R")
fig_folder <- '~/Dropbox/Projects/DurationDecoding/figure/Fig.5_caRNA/subfigs/'
# R2 compare --------------------------------------------------------------
r2.v1 <- read.csv(file = '../fig.4_modelfit/mRNA-Fit-avg-sp-v1c/result.csv',row.names = 1)
r2.v2b <- read.csv(file = '../fig.5_caRNA/caRNAFit-avg-sp-v6/caRNAFit-avg-sp-funs-v6-r2.csv',row.names = 1)
r2.v2a <- read.csv(file = '../fig.5_caRNA/emsa-caRNA-v1/v1-cmp-N.csv',row.names = 1)
all.equal(rownames(r2.v1),rownames(r2.v2a),rownames(r2.v2b))
r2.all <- data.frame(v1 = r2.v1$r2,
                     v2a=r2.v2a$r2,
                     v2b=r2.v2b$r2)
rownames(r2.all) <- rownames(r2.v1)
bks <- c(min(r2.all)-.0001,0,.5,.6,.7,max(r2.all)+.0001)
bks <- c(min(r2.all)-.0001,.7,max(r2.all)+.0001)
pheatmap(subset(r2.all,v1<0.7),color = colorRampPalette(c('black','yellow'))(length(bks)-1),
         breaks = bks,
         cluster_rows = T,cluster_cols = F,
         legend = F,treeheight_row = 0)


# derive time-change half-life  -------------------------------------------

## error model and final half-life
pd.hf.all <- read.csv(file='../../half-life/Supriya/ActDBrowser/final.genome.kdeg.csv',stringsAsFactors = F)
pd.hf.all <- pd.hf.all %>% mutate(hf=1/kdeg*60) %>% filter(hf>=0)

require(LSD)

heatscatter(x=log10(pd.hf.all$hf),y=pd.hf.all$adjR2,xlab = "hf(log10, mins)",ylab='adjR2',
            main = 'Whole genome scale half-life')
abline(v = log10(as.numeric(quantile(pd.hf.all$hf, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))),
       lty=2,col="chartreuse4")
sapply(c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
       function(x) text(log10(as.numeric(quantile(pd.hf.all$hf,x))),
                        -0.4,paste0(x*100,'%, ',signif(as.numeric(quantile(pd.hf.all$hf,x)),2)),
                        pos = 4,srt=90,col = "chartreuse4"))
abline(h = (as.numeric(quantile(pd.hf.all$adjR2, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
sapply(c(0., 0.25, 0.5, 0.75, 0.99),
       function(x) text(4,(as.numeric(quantile(pd.hf.all$adjR2,x))),
                        paste0(x*100,'%, ',signif(as.numeric(quantile(pd.hf.all$adjR2,x)),2)),
                        pos=4,col = "#3366ff"))

## filter 0.86 adjR2 
require(tidyr)
pd.hf.all <- pd.hf.all %>% filter(adjR2>=0.86)
pd.hf.pair <- pd.hf.all%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
na.fil <- complete.cases(pd.hf.pair)

x <- log10(pd.hf.pair$b1_3.0_LPS[na.fil]); y <- log10(pd.hf.pair$b2_3.0_LPS[na.fil])
heatscatter(x,y,
            xlab = 'half-life (log10, mins,rep1)',
            ylab = 'half-life (log10, mins,rep2)',
            main = 'Half-life with adjR2>=.86')
lines(loess.smooth(x,y),col=2,lwd=2)
abline(b=1,a=log10(2),lty=2)
abline(b=1,a=-log10(2),lty=2)
abline(b=1,a=0,lty=2)
lm.fit <-lm(y~x)
abline(lm.fit)
require(MASS)
abline(rlm(y~x),col=2)

abline(h = (as.numeric(quantile(y, c(0.,.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
abline(v = (as.numeric(quantile(x, c(0.,.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm=T))),
       lty=2,col="#3366ff")
sapply(c(0.,.05, 0.25, 0.5, 0.75, 0.95,1.0),
       function(x1) text(1, as.numeric(quantile(y,x1,na.rm = T)),
                         paste0(x1*100,'%, ',signif(10^as.numeric(quantile(y,x1,na.rm = T)),2)),
                         pos=4,col = "#3366ff"))
sapply(c(0.,.05, 0.25, 0.5, 0.75, 0.95,1.0),
       function(x1) text(as.numeric(quantile(x,x1,na.rm = T)),1.2,
                         paste0(x1*100,'%, ',signif(10^as.numeric(quantile(x,x1,na.rm = T)),2)),
                         pos=4,srt=90,col = "#3366ff"))
legend(x=2.2,y=3.2,legend = paste0("y=",
                                   signif(as.numeric(lm.fit$coefficients[2]),3),'x+',
                                   signif(as.numeric(lm.fit$coefficients[1]),2)),
       lty = 1,col=1)



## load error model 
require(MASS)
load(file='../fig.3_hf_hypothesis/whole.genome.error.model.hf.RData') #rm.fit 

## Predicted final half-life  
pd.hf.pair <- pd.hf.all%>% dplyr::select(ensembleID,sample,hf) %>% spread(sample,hf)
pd.hf.pair$U <- apply(pd.hf.pair[,c("b1_0.0_U","b2_0.0_U")],1,mean,na.rm=T)
pd.hf.pair$LPS_3 <- apply(pd.hf.pair[,c("b1_3.0_LPS","b2_3.0_LPS")],1,mean,na.rm=T)
pd.hf.pair$LPS_3[is.nan(pd.hf.pair$LPS_3)] <- NA 
pd.hf.pair$U[is.nan(pd.hf.pair$U)] <- NA 
pd.hf.pair <- subset(pd.hf.pair, select = -c(b1_0.0_U,b2_0.0_U,b1_3.0_LPS,b2_3.0_LPS))

colnames(pd.hf.pair)[c(1,7,2,5,8,3,6,4)]
pd.hf.pair <- pd.hf.pair[,c(1,7,2,5,8,3,6,4)]
colnames(pd.hf.pair) <- c("ensembleID","U","LPS_0.5","LPS_1","LPS_3","TNF_0.5","TNF_1","TNF_3")
pd.hf.pair$Symbol <- getSymbol(pd.hf.pair$ensembleID)

pd.hf.lb <- pd.hf.ub <- pd.hf.pair

for(i in 1:(ncol(pd.hf.pair)-2)){
  new <- data.frame(x=log10(pd.hf.pair[,i+1]))  
  sd.pred <- 10^as.numeric(predict.rlm(rlm.fit,newdata = new))
  pd.hf.lb[,i+1] <- pd.hf.pair[,i+1]-sd.pred
  pd.hf.ub[,i+1] <- pd.hf.pair[,i+1]+sd.pred
}

pd.hf.pair[,2:8] <- round(pd.hf.pair[,2:8])
pd.hf.lb[,2:8] <- round(pd.hf.lb[,2:8])
pd.hf.ub[,2:8] <- round(pd.hf.ub[,2:8])


write.csv(file=paste0(fig_folder,'Table.S.genome.hf.val.csv'),
          pd.hf.pair,quote = F,row.names = F)
write.csv(file=paste0(fig_folder,'Table.S.genome.hf.lb.csv'),
          pd.hf.lb,quote = F,row.names = F)
write.csv(file=paste0(fig_folder,'Table.S.genome.hf.ub.csv'),
          pd.hf.ub,quote = F,row.names = F)



# fig.5A - heatmap gene body quant -------------------------------------------------
cpm.peak.gene <- read.csv(file='../data/caRNA.cpm.max.geno.scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 
k_ord<-read.csv(file = '../data/mRNA.cluster.csv',stringsAsFactors = F)
sum(rownames(cpm.peak.gene) %in% k_ord$X)
cpm.peak.gene <- cpm.peak.gene[k_ord$X,]
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])

pd <- cpm.peak.gene


bks <- seq(0,1,by = 0.2)
cols <-colorRampPalette(c("azure","blue"))(6)

require(pheatmap)
pheatmap(pd,cluster_cols = F,cluster_rows = F,scale = 'none',show_rownames = F,color = cols,
         breaks = bks)


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


# plot eg profiles --------------------------------------------------------
v3.simData <- read.csv(file = '../data/bestFit.tc-caRNA.csv',stringsAsFactors = F)
names(v3.simData) <- sub("mRNA","normCnt.frac",names(v3.simData)) 
maxScale.caRNA <- plotExp(dtype = 'caRNA',scale = "geno",savetofile = F)
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- maxScale.caRNA$cond <- NULL 

p <-plotEgFit(s.genes[1:6],simData = v3.simData,expData = maxScale.mRNA,isFinal = F)
p + theme_bw() + theme(axis.title = element_blank(),axis.text = element_blank(),
                       strip.text = element_blank(),legend.position = "none")
ggsave(filename = paste0(fig_folder,'v3_eg.eps'),width = 2,height = 6)
pdf(file = 'v3.cmp.pdf')

for(i in 1:8){
  p <-plotEgFit(s.genes[(6*(i-1)+1):(6*i)],simData = v3.simData,expData = maxScale.mRNA,isFinal = F)
  pd_lab <- cbind(p$data %>% group_by(gene) %>% dplyr::summarise(r2=unique(r2)),
                  sti='LPS',type="Sim.")
  pd_lab$gene <- as.character(pd_lab$gene)
  pd_lab$hf <- round(log(2)/fit.r2[pd_lab$gene,'k_deg']);pd_lab$gene <- factor(pd_lab$gene,levels = pd_lab$gene)
  p <- p + geom_text(data= pd_lab ,
                     aes(x = 0,y=1,label=paste0('R2=',signif(r2,2),'\n hf=',hf)),colour='red',
                     hjust=0,vjust=1) 
  
  print(p)
}
dev.off()



# fig.5C - simulation  R2----------------------------------------------------

v1.simData <- read.csv(file='../fig.4_modelfit/mRNA-Fit-avg-sp-v1c/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
v1.simData.sub <- subset(v1.simData,time %in% unique(maxScale.mRNA$time))
pd <- rbind(cbind(v1.simData.sub,type='Sim.'),
            cbind(maxScale.mRNA[,colnames(v1.simData.sub)],type='Exp.'))
ph.pd<-(pd %>% unite(key,geno,type,sti,time)
        %>% spread(key,normCnt.frac))
fit.r2 <- read.csv(file = '../fig.4_modelfit/mRNA-Fit-avg-sp-v1c/result.csv',
                   row.names = 1,stringsAsFactors = F)
fit.r2.caRNA <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6/caRNAFit-avg-sp-funs-v6-r2.csv',
                         row.names = 1,stringsAsFactors = F)
filter <- fit.r2$r2<0.7
# plot violin
pd <- data.frame(r2 = c(fit.r2$r2,fit.r2.caRNA$r2),
                 type=rep(c("mRNAFit",'caRNAFit'),each=nrow(fit.r2)))
ggplot(pd,aes(factor(type),r2)) + geom_violin(aes(fill=factor(type)))
last_plot() + geom_jitter(height = 0)

pd <- data.frame(r2 = c(fit.r2$r2[filter],fit.r2.caRNA$r2[filter])
                 ,type=rep(c("mRNAFit",'caRNAFit'),each=nrow(fit.r2[filter,])),
                 gene= rep(rownames(fit.r2)[filter]),2)
s.genes <- as.character(pd$gene[(pd$r2>=0.7)])
pd <- subset(pd,gene %in% s.genes)
pd$type <- relevel(pd$type,'mRNAFit')
ggplot(pd,aes(type,r2,group=gene)) + geom_line() + 
  geom_hline(yintercept = 0.7,colour='red') + 
  geom_point(shape=1)+
  theme_bw()
ggsave(filename = paste0(fig_folder,'subfig5b_fitR2.eps'),width = 3,height = 2.5)

#tmp <- subset(pd,type=='mRNAFit' & gene %in% s.genes)
tmp <- subset(pd,type=='caRNAFit' & gene %in% s.genes)
cat(as.character(tmp[order(tmp$r2,decreasing = F),]$gene),sep = '\t')
s.genes <- as.character(tmp[order(tmp$r2,decreasing = T),]$gene)
s.genes <- as.character(unique(pd$gene))

pd %>% spread(type,r2) %>% arrange(desc(caRNAFit)) %>% select(-X2) %>%
  mutate(delta_R2=caRNAFit - mRNAFit) 
pd.2 <- pd %>% spread(type,r2) %>% arrange(desc(caRNAFit)) %>% select(gene,mRNAFit,caRNAFit)
rownames(pd.2) <- pd.2$gene; pd.2$gene <- NULL 
pheatmap(pd.2,color = c('black','yellow'),
         breaks = c(min(pd.2)-.0001,.7,max(pd.2)+.0001),
         cluster_cols = F,cluster_rows = F )


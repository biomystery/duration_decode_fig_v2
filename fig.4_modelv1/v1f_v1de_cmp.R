# preprocessing 
ev <- new.env();load('./data/clustering_cateIII_clean.Rdata',ev)
ls(ev)
v1.simData <- read.csv(file='./models/mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
v1.simData.new <- read.csv(file="./models/mRNA-Fit-avg-sp-v1e/bestFit.tc.csv",
                           stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
names(v1.simData.new) <- sub("mRNA","normCnt.frac",names(v1.simData.new)) 
v1.simData <- rbind(v1.simData,v1.simData.new) ; rm(v1.simData.new)

# subset simlutation same time points as experiments 
v1.simData.sub <- subset(v1.simData,time %in% unique(ev$pd$time))
sum(unique(ev$pd$gene) %in% unique(v1.simData.sub$gene))


pd <- rbind(cbind(subset(v1.simData.sub,gene %in% unique(ev$pd$gene)),
                  type='Sim.'),
            cbind(ev$pd[,colnames(v1.simData.sub)],
                  type='Exp.'))

#dictionary for the cluster info
dic.clust <- unique.data.frame(ev$pd %>% dplyr::select(gene,clust))
rownames(dic.clust)<- dic.clust$gene
pd <- pd%>%mutate(clust=dic.clust[gene,'clust'])


# dictionary for v1 fitting results
load('../fig.4_modelv1/data/clustering_cateIII_clean.Rdata') #new

dic.v1<- pd %>% group_by(gene)%>% 
  spread(type,normCnt.frac) %>%mutate(res=Sim.-Exp.)%>%
  dplyr::summarise(nrmsd=sqrt(sum(res^2)/length(res))/(max(Exp.)-min(Exp.)))
dic.v1<-as.data.frame(dic.v1)
rownames(dic.v1) <- dic.v1$gene
pd <- pd%>%mutate(v1=dic.v1[gene,"nrmsd"]<=0.13)

write.csv(file="./table_model.csv",dic.v1)

# export data s
tmp.var <- read.csv(file = './mRNA-Fit-avg-sp-v1e/result.csv',stringsAsFactors = F,
                    row.names = 1)
glist.v1 <- unique(pd$gene[pd$v1])
glist.v1.left <- unique(pd$gene[!pd$v1])
all(glist.v1.left %in% gene.dic$gene.2)
all(glist.v1 %in% gene.dic[rownames(tmp.var),]$gene.2) 
rownames(tmp.var) <- gene.dic[rownames(tmp.var),]$gene.2

write.csv(subset(gene.dic,gene.2%in% glist.v1.left),file='glist.v1.left.csv')
write.csv(tmp.var[unique(pd$gene[pd$v1]),],file='v1_model_good_genes.csv') 




# compare v1de vs. v1f ----------------------------------------------------
require(dplyr); require(tidyr);require(ggplot2)
load('../fig.4_modelv1/clustering_cateIII_clean.Rdata') #new

dic.v1<- pd %>% group_by(gene)%>% 
  spread(type,normCnt.frac) %>%mutate(res=Sim.-Exp.)%>%
  dplyr::summarise(nrmsd=sqrt(sum(res^2)/length(res))/(max(Exp.)-min(Exp.)))
dic.v1<-as.data.frame(dic.v1)
rownames(dic.v1) <- dic.v1$gene

nrmsd.v2 <- read.csv(file = './mRNA-Fit-avg-sp-v1f/result.csv',
                     stringsAsFactors = F)
colnames(nrmsd.v2)[1:2] <- colnames(dic.v1)
gene.dic <- read.csv(file='../data/mRNA.cluster_old.csv',stringsAsFactors = F)
gene.dic.sub <- gene.dic[,c('gene','gene.2')]; rownames(gene.dic.sub) <- gene.dic.sub$gene
nrmsd.v2$gene <- gene.dic.sub[nrmsd.v2$gene,'gene.2']

nrmsd.cmp <- rbind(data.frame(dic.v1,model="v1"),
                   data.frame(nrmsd.v2[,1:2],model='v2'))

# results.2 
nrmsd.v2 <- read.csv(file = './mRNA-Fit-avg-sp-v1f/result_2.csv',
                     stringsAsFactors = F)
colnames(nrmsd.v2)[1:2] <- colnames(dic.v1)
nrmsd.v2$gene <- gene.dic.sub[nrmsd.v2$gene,'gene.2']




# plot 
require(ggplot2)
p<-ggplot(nrmsd.cmp%>%filter(!(model=='v1'& nrmsd<.13)),
       aes(model,nrmsd,group=gene))+ geom_point()+geom_line()+
  geom_hline(yintercept = .13,colour='red')
p
table(p$data%>%spread(model,nrmsd)%>%
  mutate(improved=v2<v1)%>%select(improved))

# merged results 
res <- p$data%>%spread(model,nrmsd)
res %>% right_join(y=nrmsd.v2[,1:2],by = "gene")

write.csv(file='v1f_v1de_cmp.csv',res)

  
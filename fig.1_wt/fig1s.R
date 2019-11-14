ev <- new.env()
load(file = '../data/lrpkmD_clean_data_all.Rdata',ev)
pd <- ev$lrpkm_all
k_ord<-read.csv(file = '../data/cluster.csv',stringsAsFactors = F)

all.equal(rownames(pd),k_ord$gene)
seps <- sapply(1:8, function(x) which(k_ord$cluster==x)[1])
source(file='../auxilary_functions.R')

# fig1s-all logRpkm ------------------------------------------------------------------

pdf(file='fig1s_all_rpkm.pdf',height = 5)
pd <- pd[order(rownames(pd)),]
for(i in 1:177){
  p <- plotExpSingleGene(i,maxFrac = F)  
  p2 <- p %+% (p$data %>% group_by(gene,genotype,batch) %>% mutate(value=value/max(value))) + 
    ylab('maxFrac') + ggtitle(label = NULL)
  multiplot(p,p2,cols = 1)
  cat(i)
}
dev.off()
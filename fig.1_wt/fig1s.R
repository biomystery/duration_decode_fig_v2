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


# cmp mylist with kim -----------------------------------------------------
kim_rela <- read.csv(file='../data/Kims.RelA.genes.lfc2.txt',header = F,stringsAsFactors = F)
table(kim_rela$V1 %in% k_ord$X)

#FALSE  TRUE 
#463   100 



# old fig1.C 　donut plot ------------------------------------------------------------------
source(file='../auxilary_functions.R')
sp.th <- 2^0.5
pd.venn.ctrl <- getVennData(simplify = T)
pd.venn.ctrl.th <- getVennData(sp.th.val = 2,simplify = T)

p.1 <- donutPlot(pd.venn.ctrl,inside.v = "LT.nosp")
p.2 <- donutPlot(pd.venn.ctrl.th,inside.v = "LT.nosp")

pd <- rbind(p.1$data,p.2$data)
pd$th <- c(rep(0.5,nrow(p.1$data)),rep(1,nrow(p.2$data)))
pd <- donutPlot(pd,trans = F,group.v = 'th') 
pd <- pd + theme(strip.text.x=element_blank())
#pd + geom_text(aes(x=))

ggsave(filename = '../../../DurationDecoding/figure/Fig.1/subfigs/subfig.1C.eps',pd,width = 4,height = 2)



# batch 2
pd<- data.frame(LvT.sp = sp.mat.log2$b2.LT.ctrl,
                LvI.sp = sp.mat.log2$b2.LI.ctrl)
last_plot() %+% pd

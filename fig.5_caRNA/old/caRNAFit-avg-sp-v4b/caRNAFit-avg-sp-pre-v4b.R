# v2b- range of kdeg 
rm(list = ls())
# load data  --------------------------------------------------------------
# caRNA reads
caRNA.time <- c(0,.5,1,2,4,8)*60
caRNA.frac.avg <- read.csv(file = '../data/caRNA-cpm-max-scale.csv',header = T,row.names = 1) 

# mRNA Reads & name covert 
load(file = '../caRNAFit-avg-sp-v2/caRNAFit-avg-sp-pre-v2.Rdata')
all.equal(names(genes.kdeg),rownames(genes.normCount))

colnames(caRNA.frac.avg)
colnames(genes.normCount)


## update half-life time 
genes.kdeg.new <- read.csv(file='../../../half-life/Supriya/final.basal.hf.csv',row.names = 1)
sum(rownames(genes.kdeg.new) %in% names(genes.kdeg)) #70
g.list <- rownames(genes.kdeg.new); g.list <- g.list[g.list!='Ifnb1']

genes.normCount <- genes.normCount[g.list,]
caRNA.frac.avg <- caRNA.frac.avg[g.list,]
genes.kdeg.new <- genes.kdeg.new[g.list,]

genes.kdeg.new$kdeg <- log(2)/genes.kdeg.new$hf.mins
genes.kdeg.new$ub <- log(2)/genes.kdeg.new$ub
genes.kdeg.new$lb <- log(2)/genes.kdeg.new$lb
genes.kdeg.new <- genes.kdeg.new[,c("kdeg","ub","lb")]
colnames(genes.kdeg.new) <- colnames(genes.kdeg.new)[c(1,3,2)]
genes.kdeg.new [apply(genes.kdeg.new, 1, min)<0 ,]

tmp.val = c(median(genes.kdeg.new$kdeg),min(genes.kdeg.new$lb[genes.kdeg.new$lb>0]),max(genes.kdeg.new$ub))
(genes.kdeg.new [apply(genes.kdeg.new, 1, min)<0 ,] <- tmp.val[col(genes.kdeg.new [apply(genes.kdeg.new, 1, min)<0 ,] )])
genes.kdeg.new[c('Ccl1','Acpp',"Lcn2"),]
genes.kdeg <- genes.kdeg.new

## save data  
save(file='caRNAFit-avg-sp-pre-v4b.Rdata',list = c("caRNA.frac.avg","genes.normCount","caRNA.time","genes.timepoints","genes.kdeg"))


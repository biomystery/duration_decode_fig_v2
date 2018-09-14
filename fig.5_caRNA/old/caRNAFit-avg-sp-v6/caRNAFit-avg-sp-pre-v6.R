# v2b- range of kdeg 
rm(list = ls())
# load data  --------------------------------------------------------------
# caRNA reads
caRNA.time <- c(0,.5,1,2,4,8)*60
caRNA.frac.avg <- read.csv(file = '../data/caRNA-cpm-max-scale.csv',header = T,row.names = 1) 

# mRNA Reads & name covert 
load(file = '../old/caRNAFit-avg-sp-v2/caRNAFit-avg-sp-pre-v2.Rdata')
all.equal(names(genes.kdeg),rownames(genes.normCount))

colnames(caRNA.frac.avg)
colnames(genes.normCount)


## update half-life time 
genes.kdeg.new <- read.csv(file='../../../half-life/Supriya/final.nfkb.genes.kdeg.csv',row.names = 1)
genes.kdeg.se.new <- read.csv(file='../../../half-life/Supriya/final.nfkb.genes.kdeg.se.csv',row.names = 1)
sum(rownames(genes.kdeg.new) %in% names(genes.kdeg)) #70
g.list <- rownames(genes.kdeg.new); g.list <- g.list[g.list!='Ifnb1']

genes.normCount <- genes.normCount[g.list,]
caRNA.frac.avg <- caRNA.frac.avg[g.list,]
genes.kdeg.new <- genes.kdeg.new[g.list,]
genes.kdeg.se.new <- genes.kdeg.se.new[g.list,]


genes.kdeg.ub <- log(2)*(genes.kdeg.new+ genes.kdeg.se.new)/60 # 1/kdeg' = log(2)/kdeg ; 1/h -> 1/mins
genes.kdeg.lb <- log(2)*(genes.kdeg.new-genes.kdeg.se.new)/60
genes.kdeg <- log(2)*genes.kdeg.new/60; #in log2 scale fitted 
genes.kdeg.isNA <- is.na(genes.kdeg)

### set default values for NAs 
df <- c(median(genes.kdeg[!genes.kdeg.isNA]),range(genes.kdeg,na.rm = T))
genes.kdeg[genes.kdeg.isNA] <- df[1]
genes.kdeg.lb[genes.kdeg.isNA] <- df[2]
genes.kdeg.ub[genes.kdeg.isNA] <- df[3]

### Assign the NA to the maximum ranges 
#tmp.val = c(median(genes.kdeg.new$kdeg),min(genes.kdeg.new$lb[genes.kdeg.new$lb>0]),max(genes.kdeg.new$ub))
#(genes.kdeg.new [apply(genes.kdeg.new, 1, min)<0 ,] <- tmp.val[col(genes.kdeg.new [apply(genes.kdeg.new, 1, min)<0 ,] )])
#genes.kdeg.new[c('Ccl1','Acpp',"Lcn2"),]
#genes.kdeg <- genes.kdeg.new

## save data  
save(file='nfkb-kdeg-final.Rdata',list=c("genes.kdeg","genes.kdeg.lb","genes.kdeg.ub",'genes.kdeg.isNA'))


genes.kdeg <- genes.kdeg$U; genes.kdeg.lb <- genes.kdeg.lb$U ; genes.kdeg.ub <- genes.kdeg.ub$U
names(genes.kdeg.ub) <- names(genes.kdeg.lb) <- names(genes.kdeg) <- rownames(genes.kdeg.new)
save(file='caRNAFit-avg-sp-pre-v6.Rdata',list = c("caRNA.frac.avg","genes.normCount","caRNA.time","genes.timepoints","genes.kdeg"
                                                  ,"genes.kdeg.lb","genes.kdeg.ub"))

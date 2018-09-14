
# v2b- range of kdeg 
rm(list = ls())
# v1d processing  --------------------------------------------------------------
if(F){# caRNA reads
caRNA.time <- c(0,.5,1,2,4,8)*60
caRNA.frac.avg <- read.csv(file = './caRNA.cpm.max.geno.scale.csv',header = T,row.names = 1) 

# mRNA Reads & name covert 
ev <- new.env()
load(file = '../mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata',ev)

all.equal(names(ev$genes.kdeg),rownames(ev$genes.normCount))
all.equal(names(ev$genes.kdeg),rownames(caRNA.frac.avg))
all.equal(ev$genes.cat$gene.2,rownames(caRNA.frac.avg)) #TRUE 
rownames(caRNA.frac.avg)  <- names(ev$genes.kdeg)


fit.r2 <- read.csv(file = '../mRNA-Fit-avg-sp-v1d/result.csv',row.names = 1)
genes.selected <- rownames(fit.r2)[fit.r2$r2<0.7]
all(genes.selected %in% names(ev$genes.kdeg)) #TRUE 

ev$genes.kdeg <- ev$genes.kdeg[genes.selected]
}

# add new category  -------------------------------------------------------
load("./caRNAFit-avg-sp-pre-v6.Rdata")
genes.cat <- read.csv(file = "../mRNA-Fit-avg-sp-v1e/mRNA.cat.csv",stringsAsFactors = F)
table(genes.cat$cate)
all(genes.cat$gene%in% rownames(caRNA.frac.avg))
all(genes.cat$gene%in% rownames(caRNA.frac.avg))
ev <- new.env()
load(file = '../mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata',ev)
all(genes.cat$gene%in% names(ev$genes.kdeg))
genes.kdeg <- ev$genes.kdeg

# save data  --------------------------------------------------------------
save(file='caRNAFit-avg-sp-pre-v6.Rdata',list = c("caRNA.frac.avg","genes.normCount","caRNA.time","genes.timepoints",
                                                  "genes.kdeg","genes.cat",
                                                  "genes.kdeg.lb","genes.kdeg.ub"))




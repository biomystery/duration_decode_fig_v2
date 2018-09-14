
rm(list = ls())
# load data  --------------------------------------------------------------
# caRNA reads
caRNA.time <- c(0,.5,1,2,4,8)*60
caRNA.frac.avg <- read.csv(file = '../../data/caRNA-cpm-max-scale.csv',header = T,row.names = 1) 

# mRNA Reads & name covert 
load(file = '../../../fig.4_modelfit/mRNA-Fit-avg-sp-v1/mRNA-Fit-avg-sp-v1-pre.Rdata')
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
tmp.2 <- rownames(caRNA.frac.avg)
tmp <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
             values = grep('ENSMUSG',tmp.2,value = T), mart = ensembl)
rownames(tmp) <- tmp$ensembl_gene_id
tmp.2[grep('ENSMUSG',tmp.2)] <- tmp[tmp.2[grep('ENSMUSG',tmp.2)],"mgi_symbol"]
rownames(caRNA.frac.avg) <-  tmp.2 

all.equal(names(genes.kdeg),rownames(genes.normCount))
caRNA.frac.avg <- caRNA.frac.avg[names(genes.kdeg),]

colnames(caRNA.frac.avg)
colnames(genes.normCount)

## save data  
save(file='caRNAFit-avg-sp-pre-v2.Rdata',list = c("caRNA.frac.avg","genes.normCount","caRNA.time","genes.timepoints","genes.kdeg"))

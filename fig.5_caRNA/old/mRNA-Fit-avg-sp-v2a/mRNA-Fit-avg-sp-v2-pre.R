rm(list=ls())
load(file='../mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata')

# select genes  -----------------------------------------------------------
fit.r2 <- read.csv(file = '../mRNA-Fit-avg-sp-v1d/result.csv',row.names = 1)
genes.selected <- rownames(fit.r2)[fit.r2$r2<0.7]
genes <- read.csv(file='../mRNA-Fit-avg-sp-v1d/mRNA.cluster.csv',row.names = 2,stringsAsFactors = F)
all(genes.selected %in% genes$gene) #TRUE 
genes <- subset(genes,gene %in% genes.selected)

# half-life for selected genes  -------------------------------------------
hfs <- read.csv(file = "./Table.S.genome.hf.val.csv",stringsAsFactors = F,row.names = 1)
hfs <- hfs[rownames(genes),]
genes[grep('NA',rownames(hfs)),] # 3 genes that 
hfs <- hfs[-grep('NA',rownames(hfs)),]

## now boundary 
hfs.lb <- read.csv(file = "./Table.S.genome.hf.lb.csv",stringsAsFactors = F,row.names = 1)
hfs.lb <- hfs.lb[rownames(hfs),]
hfs.ub <- read.csv(file = "./Table.S.genome.hf.ub.csv",stringsAsFactors = F,row.names = 1)
hfs.ub <- hfs.ub[rownames(hfs),]

## check previous data 
ev <- new.env(); load('./mRNA-Fit-avg-sp-v2-pre.Rdata',ev)
head(ev$genes.kdeg.lps)
head(ev$genes.kdeg.tnf)





# current: only focus genes with half-life all available  -----------------
genes.all.available <- complete.cases(hfs[,1:7])
sum(genes.all.available) #17 genes 
nrow(hfs) # 35 
hfs <- hfs[genes.all.available,]


## processing half-life 
genes.kdeg.lps <- log(2)/hfs[,1:4]
genes.kdeg.tnf <- log(2)/hfs[,c(1,5:7)]
genes.kdeg.times <- c(0,.5,1,3,8)*60 # in mins

# 2. 8hr half-life is the same as the basal (more)
genes.kdeg.lps$LPS_8.0 <- genes.kdeg.tnf$TNF_8.0 <- genes.kdeg.lps$U

# rename the half-life 
all.equal(rownames(genes.normCount),genes.cat$gene)
subset(genes.cat,X %in% rownames(genes.kdeg.lps))
rownames(genes.cat) <- genes.cat$X
rownames(genes.kdeg.lps) <- rownames(genes.kdeg.tnf) <- genes.cat[rownames(genes.kdeg.tnf),]$gene

# alternative process -----------------------------------------------------
### assumptions: 
if(F){
  # will have the all NA matrix for v2 
  # 1. All NA genes in TNF : use measured half-life in lps
  genes.kdeg.isAllNA.lps <- apply(is.na(genes.kdeg.lps), 1, all)
  genes.kdeg.isAllNA.tnf <- apply(is.na(genes.kdeg.tnf), 1, all)
  sum(genes.kdeg.isAllNA.lps)#0
  sum(genes.kdeg.isAllNA.tnf)#3
  sum(genes.kdeg.isAllNA.lps | genes.kdeg.isAllNA.tnf)
  
  genes.kdeg.tnf[genes.kdeg.isAllNA.tnf,] <- genes.kdeg.lps[genes.kdeg.isAllNA.tnf,]
  
  
  # 2. 8hr half-life is the same as the basal (more)
  genes.kdeg.lps$LPS_8.0 <- genes.kdeg.tnf$TNF_8.0 <- genes.kdeg.lps$U
  
  # 2. 3 NA genes: extend this half-life for all time points. 
  isNA.no <- apply(is.na(genes.kdeg.lps), 1 , sum)
  for(gene in which(isNA.no==4)){
    notNa <- !is.na(genes.kdeg.lps[gene,])
    genes.kdeg.lps[gene,!notNa] <- genes.kdeg.lps[gene,notNa]
  }
  
}


# save the results  -------------------------------------------------------

save(file='mRNA-Fit-avg-sp-v2-pre.Rdata',list=c("genes.kdeg.lps",
                                                "genes.kdeg.tnf",
                                                "genes.kdeg.times",
                                                 "genes.normCount",
                                            "input.emsa.ifnarikba",
                                            "input.emsa.ifnar", "genes.timepoints"))


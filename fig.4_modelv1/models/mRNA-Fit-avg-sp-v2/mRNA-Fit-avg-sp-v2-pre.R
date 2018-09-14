## load data 
load(file='../mRNA-Fit-avg-sp-v1c/mRNA-Fit-avg-sp-v1c-pre.Rdata')
load(file = '../caRNAFit-avg-sp-v6/nfkb-kdeg-final.Rdata')
sum(genes.kdeg.isNA) 
colnames(genes.kdeg)
genes.kdeg[genes.kdeg.isNA] <- NA

write.csv(file='nfkb.genes.kdeg.csv',genes.kdeg,row.names = T)



genes.kdeg <- read.csv(file='nfkb.genes.kdeg.csv',row.names = 1)
genes.kdeg.lps <- genes.kdeg[,1:4]
genes.kdeg.tnf <- genes.kdeg[,c(1,5:7)]
genes.kdeg.times <- c(0,.5,1,3,8)*60 # in mins
genes.kdeg.isNA <- is.na(genes.kdeg)

### assumptions: 
# will have the all NA matrix for v2 
# 1. All NA genes: not fits - or use the medeam half-life 
genes.kdeg.isAllNA.lps <- apply(is.na(genes.kdeg.lps), 1, all)
genes.kdeg.isAllNA.tnf <- apply(is.na(genes.kdeg.tnf), 1, all)
genes.kdeg.tnf[genes.kdeg.isAllNA.tnf,] <- median(genes.kdeg[!genes.kdeg.isNA])
genes.kdeg.lps[genes.kdeg.isAllNA.lps,] <- median(genes.kdeg[!genes.kdeg.isNA])

# 2. 8hr half-life is the same as the basal (more)
genes.kdeg.lps$LPS_8.0 <- genes.kdeg.tnf$TNF_8.0 <- genes.kdeg.lps$U

# 2. 3 NA genes: extend this half-life for all time points. 
isNA.no <- apply(is.na(genes.kdeg.lps), 1 , sum)
for(gene in which(isNA.no==4)){
  notNa <- !is.na(genes.kdeg.lps[gene,])
  genes.kdeg.lps[gene,!notNa] <- genes.kdeg.lps[gene,notNa]
}



## nfkb input function 
data_folder<-"../data/"
input.emsa.ifnar <- read.csv(file=paste0(data_folder,'Final_EMSA_ifnar.csv'))
input.emsa.ifnarikba <- read.csv(file=paste0(data_folder,'Final_EMSA_ifnarikba.csv'))
input.emsa.ifnar[,2:4] <- input.emsa.ifnar[,2:4]/max(input.emsa.ifnar[,2:4])
input.emsa.ifnarikba[,2:4] <- input.emsa.ifnarikba[,2:4]/max(input.emsa.ifnarikba[,2:4])



save(file='mRNA-Fit-avg-sp-v1c-pre.Rdata',list=c("genes.kdeg","genes.kdeg.lb","genes.kdeg.ub",
                                                 "genes.normCount",
                                            "input.emsa.ifnarikba",
                                            "input.emsa.ifnar", "genes.timepoints"
))


###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------
numWorkers <- 32
load('../caRNAFit-avg-sp-v2/caRNAFit-avg-sp-pre-v2.Rdata')
source(file='caRNAFit-avg-sp-funs-v4.R')
# run fit  ----------------------------------------------------------------
caFit <- mclapply(names(genes.kdeg),
                           FUN = runFitCaRNADDE,
                           mc.cores = numWorkers)
names(caFit) <- unlist(lapply(caFit, function(x) x$gene))
save(file='caRNAFit-avg-sp-run-v4.Rdata',caFit)


###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------
numWorkers <- 32
load('./caRNAFit-avg-sp-pre-v4b.Rdata')
source(file='caRNAFit-avg-sp-funs-v4b.R')
# run fit  ----------------------------------------------------------------
caFit <- mclapply(rownames(genes.kdeg),
                           FUN = runFitCaRNADDE,
                           mc.cores = numWorkers)
names(caFit) <- unlist(lapply(caFit, function(x) x$gene))
save(file='caRNAFit-avg-sp-run-v4.Rdata',caFit)


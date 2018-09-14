###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------
numWorkers <- 48
load('caRNAFit-avg-sp-pre-v2.Rdata')
source(file='caRNAFit-avg-sp-funs-v2.R')
# run fit  ----------------------------------------------------------------
caFit <- mclapply(names(genes.kdeg),
                           FUN = runFitCaRNA,
                           mc.cores = numWorkers)
save(file='caRNAFit-avg-sp-run-v2.Rdata',caFit)

names(caFit) <- unlist(lapply(caFit, function(x) x$gene))
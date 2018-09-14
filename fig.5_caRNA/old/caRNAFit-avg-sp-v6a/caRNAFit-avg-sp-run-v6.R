###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------

load('./caRNAFit-avg-sp-pre-v6.Rdata')
source(file='caRNAFit-avg-sp-funs-v6.R')
require(parallel)
numWorkers <- length(genes.kdeg)

# run fit  ----------------------------------------------------------------
caFit <- mclapply(names(genes.kdeg),
                           FUN = runFitCaRNADDE,
                           mc.cores = numWorkers)
names(caFit) <- unlist(lapply(caFit, function(x) x$gene))
save(file='caRNAFit-avg-sp-run-v6.Rdata',caFit)

source(file='./caRNAFit-avg-sp-post-v6.R')
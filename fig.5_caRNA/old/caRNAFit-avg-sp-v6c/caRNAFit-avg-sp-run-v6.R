###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------

load('./caRNAFit-avg-sp-pre-v6.Rdata')
source(file='caRNAFit-avg-sp-funs-v6.R')
require(parallel)
numWorkers <- 48

# run fit  ----------------------------------------------------------------
new.genes <- read.csv(file='../mRNA-Fit-avg-sp-v1e/mRNA.cluster.csv',stringsAsFactors = F)

all(new.genes$X%in% names(genes.kdeg))
caFit <- mclapply(new.genes$X,
                           FUN = runFitCaRNADDE,
                           mc.cores = numWorkers)
names(caFit) <- unlist(lapply(caFit, function(x) x$gene))
save(file='caRNAFit-avg-sp-run-v6.Rdata',caFit)

source(file='./caRNAFit-avg-sp-post-v6.R')

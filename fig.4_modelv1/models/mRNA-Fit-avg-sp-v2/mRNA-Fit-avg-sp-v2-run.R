## load data 
load(file='mRNA-Fit-avg-sp-v1c-pre.Rdata')#,list=c("genes.kdeg","genes.normCount",
                                    #        "input.emsa.ifnarikba",
                                    #        "input.emsa.ifnar", "genes.timepoints"
source(file = "mRNA-Fit-avg-sp-v1c-funs.R") #rowMax, runModel.v1, runOptim.v1

require(parallel)
numWorkers <- 32


# main 
# round 2: fit control + knockout -----------------------------------
for(N in 1:6){
  v1.fit.MaxScale <- mclapply(rownames(genes.normCount),
                              FUN = function(x) {
                                lb <- c(-3,-3,-3,log10(genes.kdeg.lb[x]));
                                ub<- c(3,3,3,log10(genes.kdeg.ub[x]));
                                initP <- c(-1,-1,-1,log10(genes.kdeg[x]))
                                
                                runOptim.v1(initP = initP,lb = lb,ub = ub,
                                            g.name = x,input.nfkb = input.emsa.ifnar,
                                            input.nfkb.mt = input.emsa.ifnarikba)},
                              mc.cores = numWorkers)
  
  save(file = paste0('mRNA-Fit-avg-sp-v1c-run-N',N,'.Rdata'),
       list=c("v1.fit.MaxScale","N"))
}


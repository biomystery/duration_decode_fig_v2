###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------
rm(list=ls());data_folder <- "./data/"
require(parallel)
require(deSolve)
numWorkers <- 48


rowMax <- function(x) apply(x, 1, max)

runModelFromCaRNA<- function (caRNA_input,pars,times=manytp){
  
  kt<- pars$kt; kdeg<-pars$kdeg;
  
  #mRNAini <- c(mRNA= caRNA_input$val[1]*kt/kdeg) # init,ss assumption
  mRNAini <- c(mRNA=0)
  
  ## subroutine: The model
  odeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      caRNA <-approxfun(caRNA_input$time,caRNA_input$val)(Time)
      dmRNA    <- kt*caRNA - kdeg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, odeModel, pars)
  out<-as.data.frame(out)
}

runOptimFromCaRNA <- function(initP,lb,ub,caRNA,gname){ 
   obj <- function(pars){
    pars <- 10^pars # log10 scale
    out<-  sapply(1:2, function(x)
      runModelFromCaRNA(caRNA_input = data.frame(time=caRNA.time, 
                                                 val =unlist(caRNA[,x])),
                        data.frame(kt=pars[1],kdeg=genes.kdeg[gname,1]),times),simplify = F)
    idx <- which(out[[1]]$time %in% mRNAprofile$time); idx[length(idx)] <- idx[length(idx)]-1
    out<- unlist(lapply(out, function(x) x[idx,"mRNA"]))
    out <- out /max(out) 
    rmsd <- sqrt(sum((out - mRNAprofile$val)^2)/length(mRNAprofile$time))
  }
  
  
  times <- seq(0,max(caRNA.time),by=0.1);
  mRNAprofile<-data.frame(time=genes.timepoints,
                          val = genes.normCount[gname,1:14])
    
  try(fit <- nlminb(start = initP,obj,lower=lb,upper =ub));fit$par <- 10^fit$par  
  out<-  sapply(1:2, function(x)
    runModelFromCaRNA(caRNA_input = data.frame(time=caRNA.time, 
                                               val =unlist(caRNA[,x])),
                      data.frame(kt=fit$par,kdeg=genes.kdeg[gname,1]),times),simplify = F)
  
  fit$expdata<-mRNAprofile
  fit$bestFit <- out
  idx <- which(out[[1]]$time %in% mRNAprofile$time); idx[length(idx)] <- idx[length(idx)]-1
  fit$simScale <- max(unlist(lapply(out, function(x) x[idx,"mRNA"])))
  fit$gene <- gname
  fit$residue <- unlist(lapply(out, function(x) x[idx,"mRNA"]))/fit$simScale - mRNAprofile$val
  return(fit)
  
}


runFitCaRNA <- function(gname){
  #gname <- rownames(genes.normCount)[1]
  caRNA_tmp <- matrix(cpmfrac.avg.gene[gname,],ncol = 4)
  runOptimFromCaRNA(initP = .1,lb = -3,ub = 3,caRNA = caRNA_tmp,gname = gname)
}


# load data  --------------------------------------------------------------
caRNA.time <- c(0,.5,1,2,4,8)*60

cpmfrac.avg.gene <- read.csv(file = './data/caRNA-cpm-max-scale.csv',header = T,row.names = 1) 
load(file='../fig.1_wt/lrpkmD_clean_data_all.Rdata')
ls()

ev <- new.env()
load(file='../fig.4_modelfit/data/v4-hf.Rdata',envir = ev)
genes.hf <- read.csv(file="../fig.4_modelfit/data/v4-hf-final.csv",stringsAsFactors = F,row.names = 1)
genes.kdeg <- log(2)/(10^genes.hf) ; 
all.equal(rownames(genes.kdeg),names(ev$nfkbgenes))
rm(ev)

## lrpkm 
load(file="../fig.4_modelfit/data/lrpkmD_clean_data_all.Rdata")
genes.normCount <- 2^lrpkm_all[1:28]; rm(lrpkm_all)
genes.normCount <- genes.normCount[,c(1:7,15:21,8:14,22:28)]
genes.normCount <- genes.normCount/rowMax(genes.normCount[,1:14]) #normnizsed by wt 
genes.timepoints <- c(0,.5,1,2,3,5,8)*60 
colnames(genes.normCount)



# run fit  ----------------------------------------------------------------

caFit <- mclapply(rownames(genes.kdeg),
                           FUN = runFitCaRNA,
                           mc.cores = numWorkers)
save(file='caFit.Rdata',caFit)

## load data 
data_folder<-"./data/"
library(gplots);library(RColorBrewer)
#library(pheatmap)
require(parallel)
require(deSolve)
numWorkers <- 48

## pre-defined functions 
rowMax <- function(x) apply(x, 1, max)
runModel.v1<- function (nfkb_input,pars,times=manytp){
  k0 <- pars$k0;kt<-pars$kt;k_deg<-pars$k_deg;kd <- pars$kd;
    mRNAini <- c(mRNA= (kt*nfkb_input$val[1]/(kd+nfkb_input$val[1])+k0)/k_deg) # init,ss assumption
  
  ## subroutine: The model
  odeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      nfkb <-approxfun(nfkb_input$time,nfkb_input$val)(Time)
      dmRNA    <- kt*nfkb/(kd+nfkb)+k0-k_deg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, odeModel, pars)
  out<-as.data.frame(out)
  out$mRNA[nrow(out)] <- out$mRNA[nrow(out)-1]
  out
}

runOptim.v1 <- function(initP,lb,ub,g.name,input.nfkb){ 
  
  obj <- function(pars){
    pars <- 10^pars # log10 scale
    out <- sapply(c("LPS","TNF"), function(x)
      runModel.v1(nfkb_input = data.frame(time= input.nfkb$Time,val = input.nfkb[,x]),
                  pars= data.frame(k0=pars[1],kt=pars[2],kd=pars[3],k_deg=genes.kdeg[g.name,1]),
                  times = times),simplify = F)
    out<- unlist(lapply(out, function(x) subset(x,subset =out$LPS$time %in% mRNAprofile$time)$mRNA))
    out <- out /max(out) 
    rmsd <- sqrt(sum((out - mRNAprofile$val)^2)/length(mRNAprofile$time))
  }
  
  times <- seq(0,max(input.nfkb$Time),by=0.1);
  mRNAprofile<-data.frame(time=genes.timepoints,
                          val = genes.normCount[g.name,1:14])
  
  try(fit <- nlminb(start = initP,obj,lower=lb,upper =ub));fit$par <- 10^fit$par
  out <- sapply(c("LPS","TNF"), function(x)
    runModel.v1(nfkb_input = data.frame(time= input.nfkb$Time,val = input.nfkb[,x]),
                pars= data.frame(k0=fit$par[1],kt=fit$par[2],kd=fit$par[3],k_deg=genes.kdeg[g.name,1]),
                times = times),simplify = F)

  fit$expdata<-mRNAprofile
  fit$bestFit <- out
  fit$simScale <- max(unlist(lapply(out, function(x) subset(x,subset =out$LPS$time %in% mRNAprofile$time)$mRNA)))
  fit$gene <- g.name
  return(fit)
}

# main --------------------------------------------------------------------

## read half-lifes 
ev <- new.env()
load(file=paste0(data_folder,'v4-hf.Rdata'),envir = ev)
genes.hf <- read.csv(file=paste0(data_folder,"v4-hf-final.csv"),stringsAsFactors = F,row.names = 1)
genes.kdeg <- log(2)/(10^genes.hf) ; 
rm(genes.hf); 
all.equal(rownames(genes.kdeg),names(ev$nfkbgenes))
rm(ev)

## nfkb input function 
input.emsa.ifnar <- read.csv(file=paste0(data_folder,'Final_EMSA_ifnar.csv'))
input.emsa.ifnarikba <- read.csv(file=paste0(data_folder,'Final_EMSA_ifnarikba.csv'))
input.emsa.ifnar[,2:4] <- input.emsa.ifnar[,2:4]/max(input.emsa.ifnar[,2:4])
input.emsa.ifnarikba[,2:4] <- input.emsa.ifnarikba[,2:4]/max(input.emsa.ifnarikba[,2:4])


## lrpkm 
load(file=paste0(data_folder,'lrpkmD_clean_data_all.Rdata'))
genes.normCount <- 2^lrpkm_all[1:28]; rm(lrpkm_all)
genes.normCount <- genes.normCount[,c(1:7,15:21,8:14,22:28)]
genes.normCount <- genes.normCount/rowMax(genes.normCount[,1:14]) #normnizsed by wt 
genes.timepoints <- c(0,.5,1,2,3,5,8)*60 
colnames(genes.normCount)

# round 1: fit control predict knockout -----------------------------------
lb <- c(-3,-3,-3);ub<- c(3,3,3);initP <- c(0.1,0.1,.1)

v1.fit.rescale <- mclapply(rownames(genes.kdeg),
                   FUN = function(x) runOptim.v1(initP = initP,lb = lb,ub = ub,
                                                 g.name = x,input.nfkb = input.emsa.ifnar),
                   mc.cores = numWorkers)


# round 2: fit control + knockouts ----------------------------------------

save(file = 'v1.fit.rescale.RData',v1.fit.rescale)


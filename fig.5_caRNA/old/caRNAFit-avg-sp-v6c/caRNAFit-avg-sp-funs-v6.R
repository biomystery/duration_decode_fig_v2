###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------
require(parallel)
require(deSolve)
require(pracma)
# require data: caRNA.frac.avg, genes.kdeg,genes.timepoints,caRNA.time.

rowMax <- function(x) apply(x, 1, max)

runModelFromCaRNADDE<- function (caRNA_input,pars,times=manytp){
  
  kt<- pars$kt; kdeg<-pars$kdeg;tau <- pars$tau;
  
  mRNAini <- c(mRNA= caRNA_input$val[1]*kt/kdeg) # init,ss assumption
#  mRNAini <- c(mRNA=0)
  
  ## subroutine: The model
  ddeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      #caRNA <-ifelse(Time<= tau,caRNA_input$val[1],approxfun(caRNA_input$time,caRNA_input$val)(Time-tau))
      caRNA <-ifelse(Time<= tau,caRNA_input$val[1],pchipfun(caRNA_input$time,caRNA_input$val)(Time-tau))
      
      dmRNA    <- kt*caRNA - kdeg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, ddeModel, pars)
  out<-as.data.frame(out)
}

runOptimFromCaRNADDE <- function(initP,lb,ub,caRNA,gname){ 
   obj <- function(pars){
    pars <- 10^pars # log10 scale
    out<-  sapply(1:4, function(x)
      runModelFromCaRNADDE(caRNA_input = data.frame(time=caRNA.time, 
                                                 val =unlist(caRNA[,x])),
                        data.frame(kt=pars[1],kdeg=pars[2],tau=pars[3]),times),simplify = F)
    idx <- which(out[[1]]$time %in% mRNAprofile$time); idx[length(idx)] <- idx[length(idx)]-1
    out<- unlist(lapply(out, function(x) x[idx,"mRNA"]))
    out[1:14] <- out[1:14] /max(out[1:14]) 
    out[15:28] <- out[15:28] /max(out[15:28]) 
    nrmsd <- sqrt(sum((out - mRNAprofile$val)^2)/length(mRNAprofile$time))/(max(mRNAprofile$val)-min(mRNAprofile$val))
  }
  
  
  times <- seq(0,max(caRNA.time),by=0.1);
  mRNAprofile<-data.frame(time=genes.timepoints,
                          val = genes.normCount[gname,1:28])
    
  try(fit <- nlminb(start = initP,obj,lower=lb,upper =ub));fit$par <- 10^fit$par  
  out<-  sapply(1:4, function(x)
    runModelFromCaRNADDE(caRNA_input = data.frame(time=caRNA.time, 
                                               val =unlist(caRNA[,x])),
                      data.frame(kt=fit$par[1],kdeg=fit$par[2],tau=fit$par[3]),times),simplify = F)
  
  fit$expdata<-mRNAprofile
  fit$bestFit <- out
  idx <- which(out[[1]]$time %in% mRNAprofile$time); idx[length(idx)] <- idx[length(idx)]-1
  fit$simScale.wt <- max(unlist(lapply(out, function(x) x[idx,"mRNA"]))[1:14])
  fit$simScale.mt <- max(unlist(lapply(out, function(x) x[idx,"mRNA"]))[15:28])
  fit$gene <- gname
  tmp <- unlist(lapply(out, function(x) x[idx,"mRNA"])); tmp[1:14]<- tmp[1:14]/fit$simScale.wt
  tmp[15:28] <- tmp[15:28]/fit$simScale.mt
  fit$residue <- tmp - mRNAprofile$val
  return(fit)
  
}

runFitCaRNADDE <- function(gname){
  #gname <- rownames(genes.normCount)[1]
  caRNA_tmp <- matrix(caRNA.frac.avg[gname,],ncol = 4)
  runOptimFromCaRNADDE(initP = c(.1,log10(genes.kdeg[gname]),1),
                       lb = c(-3,log10(genes.kdeg.lb[gname]),-5),
                       ub = c(3,log10(genes.kdeg.ub[gname]),3),
                       caRNA = caRNA_tmp,gname = gname)
}


calR2 <- function(x){
  tmp <- caFit[[x]]
  1 - var(tmp$residue)/var(tmp$expdata$val)
  
}

plotFitInd <- function(fit){
  geno <- unique(substr(colnames(caRNA.frac.avg),1,3))
  tls <- list(paste0(geno[1],', ',fit$gene,',r2=',signif(caFit.R2[fit$gene],2)),
              paste0(geno[2],',nrmsd=',signif(fit$objective,2)),
              paste0(geno[3],',sim(l) exp(d)'),
              paste0(geno[4],'caRNA(b) mRNA(r)'))
  ylms <- max(sapply(1:4, function(i) max(c(fit$bestFit[[i]]$mRNA/fit$simScale,fit$expdata[(i*7-6):(i*7),2]),na.rm = T)))
  par(mfrow=c(2,2))
  for(i in 1:4){
    pd <- fit$bestFit[[i]]; pd$mRNA <- pd$mRNA/fit$simScale
    pd.expdata <- fit$expdata[(i*7-6):(i*7),]
    plot(pd,type='l',ylim=c(0,ylms),xlab='Time (h)',ylab='max frac',col=2,
         main=tls[[i]])
    points(pd.expdata,pch=16,col=2)
    lines(caRNA.time,caRNA.frac.avg[fit$gene,(i*6-5):(i*6)],col=1,type = "b",pch=16)
  }
  par(mfrow=c(1,1))
}
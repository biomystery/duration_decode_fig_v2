###### Fit1: use caRNAseq -> fit (RNAseq) -> half-lifes =------
require(parallel)
require(deSolve)
require(pracma)
# require data: caRNA.frac.avg, genes.kdeg,genes.timepoints,caRNA.time.

rowMax <- function(x) apply(x, 1, max)

runModelFromCaRNADDE<- function (caRNA_input,k_deg_input,
                                 pars,times=manytp){
  
  kt<- pars$kt; tau <- pars$tau;
  
  mRNAini <- c(mRNA= caRNA_input$val[1]*kt/k_deg_input$val[1]) # init,ss assumption
#  mRNAini <- c(mRNA=0)
  
  ## subroutine: The model
  ddeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      
      caRNA <-ifelse(Time<= tau,caRNA_input$val[1],pchipfun(caRNA_input$time,caRNA_input$val)(Time-tau))
      k_deg <- pchipfun(k_deg_input$time, k_deg_input$val)(Time)
      
      dmRNA    <- kt*caRNA - k_deg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, ddeModel, pars)
  out<-as.data.frame(out)
}

runOptimFromCaRNADDE <- function(initP,lb,ub,caRNA,kdegs,gname){ 
   obj <- function(pars){
    pars <- 10^pars # log10 scale
    out<-  sapply(1:4, function(x)
      runModelFromCaRNADDE(caRNA_input = data.frame(time=caRNA.time, 
                                                 val =unlist(caRNA[,x])),
                           k_deg_input = data.frame(time=caRNA.time, 
                                                    val =unlist(kdegs[,x])), 
                        data.frame(kt=pars[1],tau=pars[2]),times),simplify = F)
    idx <- which(out[[1]]$time %in% mRNAprofile$time); idx[length(idx)] <- idx[length(idx)]-1
    out<- unlist(lapply(out, function(x) x[idx,"mRNA"]))
    out[1:14] <- out[1:14] /max(out[1:14]) 
    out[15:28] <- out[15:28] /max(out[15:28]) 
    rmsd <- sqrt(sum((out - mRNAprofile$val)^2)/length(mRNAprofile$time))
  }
  
  
  times <- seq(0,max(caRNA.time),by=0.1);
  mRNAprofile<-data.frame(time=genes.timepoints,
                          val = genes.normCount[gname,1:28])
    
  try(fit <- nlminb(start = initP,obj,lower=lb,upper =ub));fit$par <- 10^fit$par  
  out<-  sapply(1:4, function(x)
    runModelFromCaRNADDE(caRNA_input = data.frame(time=caRNA.time, 
                                               val =unlist(caRNA[,x])),
                         k_deg_input = data.frame(time=caRNA.time, 
                                                  val =unlist(kdegs[,x])), 
                      data.frame(kt=fit$par[1],tau=fit$par[2]),times),simplify = F)
  
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
  kdeg_tmp <- data.frame(kdeg_lps_wt = genes.kdeg.lps[gname,],
                         kdeg_tnf_wt = genes.kdeg.tnf[gname,],
                         kdeg_lps_mt = genes.kdeg.lps[gname,],
                         kdeg_tnf_mt = genes.kdeg.tnf[gname,])
  runOptimFromCaRNADDE(initP = c(.1,1),
                       lb = c(-3,-5),
                       ub = c(3,3),
                       caRNA = caRNA_tmp,kdegs = kdeg_tmp,gname = gname)
}


calR2 <- function(x){
  tmp <- caFit[[x]]
  data.frame(r2=1 - var(tmp$residue)/var(tmp$expdata$val),
             kt=tmp$par[1],
             tau=tmp$par[2])
}

plotFitInd <- function(g='Jak2'){
  geno <- unique(substr(colnames(caRNA.frac.avg),1,3))
  tls <- list(paste0(geno[1],', ',g,',r2=',signif(caFit.R2[[g,'r2']],2)),
              paste0(geno[2],',rmsd=',signif(fit[[g]]$objective,2)),
              paste0(geno[3],',sim(l) exp(d)'),
              paste0(geno[4],'caRNA(b) mRNA(r)'))
  
  par(mfrow=c(2,2))
  for(i in 1:4){
    pd <- fit[[g]]$bestFit[[i]]; 
    if (i %in% 1:2)
      pd$mRNA <- pd$mRNA/fit[[g]]$simScale.wt
    else
      pd$mRNA<-pd$mRNA/fit[[g]]$simScale.mt
    pd.expdata <- fit[[g]]$expdata[(i*7-6):(i*7),]
    plot(pd,type='l',ylim=c(0,1.2),xlab='Time (h)',ylab='max frac',col=2,
         main=tls[[i]])
    points(pd.expdata,pch=16,col=2)
    lines(caRNA.time,caRNA.frac.avg[g,(i*6-5):(i*6)],col=1,type = "b",pch=16)
  }
  par(mfrow=c(1,1))
}

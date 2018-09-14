## load data 
rm(list = ls())
data_folder<-"./"
library(gplots);library(RColorBrewer)
#library(pheatmap)
require(parallel)
require(deSolve)


## pre-defined functions 
rowMax <- function(x) apply(x, 1, max)
runModel.v1.2<- function (nfkb_input,pars,times=manytp){
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

runOptim.v1.2 <- function(initP,lb,ub,g.name,input.nfkb.ctl,
                          input.nfkb.mt){ 
  
  obj <- function(pars){
    pars <- 10^pars # log10 scale
    out.ctl <- sapply(c("LPS","TNF"), function(x)
      runModel.v1.2(nfkb_input = data.frame(time= input.nfkb.ctl$Time,val = input.nfkb.ctl[,x]),
                    pars= data.frame(k0=pars[1],kt=pars[2],kd=pars[3],k_deg=genes.kdeg[g.name,1]),
                    times = times),simplify = F)
    out.mt <- sapply(c("LPS","TNF"), function(x)
      runModel.v1.2(nfkb_input = data.frame(time= input.nfkb.mt$time,val = input.nfkb.mt[,x]),
                    pars= data.frame(k0=pars[1],kt=pars[2],kd=pars[3],k_deg=genes.kdeg[g.name,1]),
                    times = times),simplify = F)
    
    out.ctl<- unlist(lapply(out.ctl, function(x) subset(x,subset =out.ctl$LPS$time %in% mRNAprofile$time)$mRNA))
    out.mt<- unlist(lapply(out.mt, function(x) subset(x,subset =out.mt$LPS$time %in% mRNAprofile$time)$mRNA))
    out <- c(out.ctl,out.mt)
    out <- out /max(out.ctl) 
    rmsd <- sqrt(sum((out - mRNAprofile$val)^2)/length(mRNAprofile$time))
  }
  
  times <- seq(0,max(input.nfkb.ctl$Time),by=0.1);
  mRNAprofile<-data.frame(time=genes.timepoints,
                          val = genes.normCount[g.name,1:28])
  
  try(fit <- nlminb(start = initP,obj,lower=lb,upper =ub));fit$par <- 10^fit$par
  out.ctl <- sapply(c("LPS","TNF"), function(x)
    runModel.v1.2(nfkb_input = data.frame(time= input.nfkb.ctl$Time,val = input.nfkb.ctl[,x]),
                  pars= data.frame(k0=fit$par[1],kt=fit$par[2],kd=fit$par[3],k_deg=genes.kdeg[g.name,1]),
                  times = times),simplify = F)
  out.mt <- sapply(c("LPS","TNF"), function(x)
    runModel.v1.2(nfkb_input = data.frame(time= input.nfkb.mt$time,val = input.nfkb.mt[,x]),
                  pars= data.frame(k0=fit$par[1],kt=fit$par[2],kd=fit$par[3],k_deg=genes.kdeg[g.name,1]),
                  times = times),simplify = F)
  out <- c(out.ctl,out.mt)  
  fit$expdata<-mRNAprofile
  fit$bestFit <- out
  fit$simScale <- max(unlist(lapply(out.ctl, function(x) subset(x,subset =out$LPS$time %in% mRNAprofile$time)$mRNA)))
  fit$gene <- g.name
  return(fit)
}


# main --------------------------------------------------------------------

## read half-lifes 
ev <- new.env()
load(file=paste0(data_folder,'v4-hf.Rdata'),envir = ev)
genes.hf <- read.csv(file=paste0(data_folder,"v4-hf-final.csv"),stringsAsFactors = F,row.names = 1)
genes.kdeg <- log(2)/(10^genes.hf) ; 
#rm(genes.hf); 
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

# model fit 
load(file = 'v1.fit.2.n1.RData')
ls()
fun.calR2 <- function(x){
  #tmp <- v1.fit.1[[1]]  
  res.idx <- which(x$bestFit$LPS$time %in% x$expdata$time) 
  res.lps <- (x$bestFit$LPS$mRNA/x$simScale)[res.idx] - x$expdata$val[1:length(res.idx)]
  res.tnf <- (x$bestFit$TNF$mRNA/x$simScale)[res.idx] - x$expdata$val[-(1:length(res.idx))]
  sqrt(sum((c(res.tnf,res.lps))^2)/length(x$expdata$time)) ==x$objective
  1 - var(c(res.lps,res.tnf))/var(x$expdata$val)
}


fun.plotFit <- function(x){
  #x <- v1.fit.2[[1]]
  plot(x$bestFit[[1]]$time,x$bestFit[[1]]$mRNA/x$simScale,
       type='l',xlab="time(mins)",ylab=x$gene,ylim=c(0,1),
       main = paste0("r2=",signif(x$r2,2),
                     ',obj=',signif(x$objective,2)))
  lines(x$bestFit[[2]]$time,x$bestFit[[2]]$mRNA/x$simScale,
        type='l',col=2)
  lines(x$bestFit[[3]]$time,x$bestFit[[3]]$mRNA/x$simScale,
        type='l',col=1,lty=2)
  lines(x$bestFit[[4]]$time,x$bestFit[[4]]$mRNA/x$simScale,
        type='l',col=2,lty =2)
  
  points(x$expdata,pch=rep(c(16,1),each=14),col=rep(c(1,2,1,2),each=7))
  
  text(300,.5,paste0("kdeg=",signif(genes.kdeg[x$gene,1],2)))
  #nrow(tmp$expdata)
}

v1.fit.2.r2 <- unlist(lapply(v1.fit.2, fun.calR2))
for(i in 1:177) v1.fit.2[[i]]$r2 <- v1.fit.2.r2[i]
v1.fit.2.obj <- unlist(lapply(v1.fit.2, function(x) x$obj))
  
# plot the results 
tmp.ord<- order(v1.fit.2.r2,decreasing = T)
tmp.ord.2<- order(v1.fit.2.obj)

pdf(file="v1.2.res.pdf")
par(mfrow=c(2,2))
hist(v1.fit.2.r2)
hist(v1.fit.2.obj)
par(mfrow=c(3,3))
sapply(1:9, function(x) fun.plotFit(v1.fit.2[[tmp.ord[x]]]))
sapply(169:177, function(x) fun.plotFit(v1.fit.2[[tmp.ord[x]]]))
sapply(1:9, function(x) fun.plotFit(v1.fit.2[[tmp.ord.2[x]]]))
sapply(169:177, function(x) fun.plotFit(v1.fit.2[[tmp.ord.2[x]]]))
par(mfrow=c(1,1))
dev.off()




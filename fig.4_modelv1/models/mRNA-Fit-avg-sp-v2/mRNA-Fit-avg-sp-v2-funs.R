require(deSolve)

## pre-defined functions 
rowMax <- function(x) apply(x, 1, max)
runModel.v1<- function (nfkb_input,pars,times=manytp){
  k0 <- pars$k0;kt<-pars$kt;k_deg<-pars$k_deg;kd <- pars$kd;
  mRNAini <- c(mRNA= (kt*nfkb_input$val[1]^N/(kd^N+nfkb_input$val[1]^N)+k0)/k_deg) # init,ss assumption
  
  ## subroutine: The model
  odeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      nfkb <-approxfun(nfkb_input$time,nfkb_input$val)(Time)
      dmRNA    <- kt*nfkb^N/(kd^N+nfkb^N)+k0-k_deg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, odeModel, pars)
  out<-as.data.frame(out)
  out$mRNA[nrow(out)] <- out$mRNA[nrow(out)-1]
  out
}
runOptim.v1 <- function(initP,lb,ub,g.name,input.nfkb,input.nfkb.mt){ 
  
  obj <- function(pars){
    pars <- 10^pars # log10 scale
    # sim contrl
    out <- sapply(c("LPS","TNF"), function(x)
      runModel.v1(nfkb_input = data.frame(time= input.nfkb$Time,val = input.nfkb[,x]),
                  pars= data.frame(k0=pars[1],kt=pars[2],kd=pars[3],k_deg=pars[4]),
                  times = times),simplify = F)
    out<- unlist(lapply(out, function(x) subset(x,subset =out$LPS$time %in% mRNAprofile$time)$mRNA))
    # sim mt 
    out.mt <- sapply(c("LPS","TNF"), function(x)
      runModel.v1(nfkb_input = data.frame(time= input.nfkb.mt$Time,val = input.nfkb.mt[,x]),
                  pars= data.frame(k0=pars[1],kt=pars[2],kd=pars[3],k_deg=pars[4]),
                  times = times),simplify = F)
    out.mt<- unlist(lapply(out.mt, function(x) subset(x,subset =out.mt$LPS$time %in% mRNAprofile$time)$mRNA))
    
    # combine ctrl. + mt 
    out <- c(out,out.mt) /max(out) # normalize to max in contrl 
    rmsd <- sqrt(sum((out - mRNAprofile$val)^2)/length(mRNAprofile$time))
  }
  
  times <- seq(0,max(input.nfkb$Time),by=0.1);
  mRNAprofile<-data.frame(time=rep(genes.timepoints,4), #4 time course
                          val = genes.normCount[g.name,1:28])
  
  try(fit <- nlminb(start = initP,obj,lower=lb,upper =ub));fit$par <- 10^fit$par
  out  <- sapply(c("LPS","TNF"), function(x)
    runModel.v1(nfkb_input = data.frame(time= input.nfkb$Time,val = input.nfkb[,x]),
                pars= data.frame(k0=fit$par[1],kt=fit$par[2],kd=fit$par[3],k_deg=fit$par[4]),
                times = times),simplify = F)

  out.mt  <- sapply(c("LPS","TNF"), function(x)
    runModel.v1(nfkb_input = data.frame(time= input.nfkb.mt$Time,val = input.nfkb.mt[,x]),
                pars= data.frame(k0=fit$par[1],kt=fit$par[2],kd=fit$par[3],k_deg=fit$par[4]),
                times = times),simplify = F)
  fit$bestFit.ctrl <- out ; fit$bestFit.mt <- out.mt  
  fit$expdata<-mRNAprofile
  fit$simScale <- max(unlist(lapply(out, function(x) subset(x,subset =out$LPS$time %in% mRNAprofile$time)$mRNA)))
  fit$gene <- g.name
  return(fit)
}

## post-process functions 
fun.calR2 <- function(x){
  #tmp <- v1.fit.1[[1]]  
  res.idx <- which(x$bestFit.ctrl$LPS$time %in% x$expdata$time) 
  res.all <- c(x$bestFit.ctrl$LPS$mRNA[res.idx],x$bestFit.ctrl$TNF$mRNA[res.idx],
               x$bestFit.mt$LPS$mRNA[res.idx],x$bestFit.mt$TNF$mRNA[res.idx])/x$simScale; 
  res.all <- res.all - x$expdata$val
  print(sqrt(sum(res.all^2)/length(x$expdata$time)) ==x$objective)
  1 - var(res.all)/var(x$expdata$val)
}

fun.plotFit <- function(x){
  plot(x$bestFit.ctrl$LPS$time,x$bestFit.ctrl$LPS$mRNA/x$simScale,
       type='l',xlab="time(mins)",ylab=x$gene,ylim=c(0,1.2),
       main = paste0("r2=",signif(x$r2,2),
                     ',obj=',signif(x$objective,2)))
  lines(x$bestFit.ctrl$TNF$time,x$bestFit.ctrl$TNF$mRNA/x$simScale,
        type='l',col=2)
  lines(x$bestFit.mt$LPS$time,x$bestFit.mt$LPS$mRNA/x$simScale,
        type='l',col=1,lty=2)
  lines(x$bestFit.mt$TNF$time,x$bestFit.mt$TNF$mRNA/x$simScale,
        type='l',col=2,lty=2)
  
  points(x$expdata,pch=rep(c(16,1),each=14),col=rep(rep(c(1,2),each=7),2))
  text(300,.5,paste0("kdeg=",signif(genes.kdeg[x$gene],2)))
  #nrow(tmp$expdata)
}

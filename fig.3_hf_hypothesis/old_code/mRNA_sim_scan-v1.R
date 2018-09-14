## load data 
rm(list = ls())
require(deSolve)    

library(gplots);library(RColorBrewer)
library(pheatmap)

rowMax <- function(x) apply(x, 1, max)

###################################
# Input curves
###################################
dataFolder <- './'
ifnar_input <- read.csv(file=paste0(dataFolder,'Final_EMSA_ifnar.csv'))
ifnarikba_input <- read.csv(file=paste0(dataFolder,'Final_EMSA_ifnarikba.csv'))

emsa_max <- max(ifnar_input[,2:4])
ifnar_input[,2:4] <- ifnar_input[,2:4]/emsa_max
#ifnarikba_input[,2:4] <- ifnarikba_input[,2:4]/max(ifnarikba_input[,2:4])
ifnarikba_input[,2:4] <- ifnarikba_input[,2:4]/emsa_max
ifnar_input
ifnarikba_input

###################################
# model simulation (wt)
###################################
runModel<- function (nfkb_input,pars,times=manytp){
  
  k0 <- pars$k0;kt<-pars$kt;k_deg<-pars$k_deg;kd <- pars$kd;n <- pars$n
  
  mRNAini <- c(mRNA= (k0+kt*(nfkb_input$val[1])^n/(kd^n+(nfkb_input$val[1])^n))/k_deg) # init,ss assumption
  
  ## subroutine: The model
  odeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      nfkb <-approxfun(nfkb_input$time,nfkb_input$val)(Time)
      dmRNA    <- k0 + kt*nfkb^n/(nfkb^n+kd^n) -k_deg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, odeModel, pars)
  out<-as.data.frame(out)
}




fun.runSim<- function(n=1,hf = c(60,60*2.5),kds = rep(0.5,2)){
  #hf <- c(10,300)
  manytp <- seq(0,480,by=0.1)
  times <- ifnar_input$Time
  idx <- which(manytp%in% times)
  idx[7]<- idx[7]-1
  N <- length(hf)
  kb = 0.001;kt = 0.5;#kd=0.5;  #n=1;
  pars <- data.frame(k0 = rep(kb,N), 
                     kt = rep(kt,N),
                     kd = kds,
                     n = rep(n,N),
                     k_deg = log(2)/hf)
  
  # run LPS 
  nfkb <-data.frame(time=ifnar_input$Time,
                    val =ifnar_input$LPS)
  
  sim_result <- sapply(1:length(hf),function(x) runModel(nfkb_input = nfkb,
                                                         pars = pars[x,],times = manytp)$mRNA)
  lps_sim_result <- matrix(unlist(sim_result),nrow = N,byrow = T)
  
  #lps_sim_result<- lps_sim_result[,idx]
  # run tnf 
  nfkb <-data.frame(time=ifnar_input$Time,
                    val =ifnar_input$TNF)
  
  sim_result <- sapply(1:length(hf),function(x) runModel(nfkb_input = nfkb,
                                                         pars = pars[x,],times = manytp)$mRNA)
  tnf_sim_result <- matrix(unlist(sim_result),nrow = N,byrow = T)
  
    sim.result <- rbind(lps_sim_result,tnf_sim_result)
  matplot(manytp,lps_sim_result[1,],type = "l",
          xlab="Time(mins)",ylab="expression",xaxt="n",
          main=paste0("n=",n,',hf=',hf[1],',kd=',kds[1]),
          lty=2,lwd=2)
  axis(1,at=nfkb$time);grid()
  lines(manytp,tnf_sim_result[1,],type = "l",lwd=2)
  plot(manytp,lps_sim_result[2,],type = "l",
       xaxt="n",xlab = "Time(mins)",ylab = "expression",
       main=paste0("n=",n,",hf=",hf[2],',kd=',kds[2]),lty=2,lwd=2)
  axis(1,at=nfkb$time);grid()
  lines(manytp,tnf_sim_result[2,],type = "l",lwd=2)
}

pdf(file='hf-theory-10vs300-kd.pdf',width =  7/3*5)
par(mfrow=c(3,6))
fun.runSim(1,hf = c(10,300),kds=rep(.2,2))
fun.runSim(1,hf = c(10,300),kds=rep(.8,2))
fun.runSim(1,hf = c(10,300),kds=rep(5,2))
fun.runSim(2,hf=c(10,300),kds=rep(.2,2))
fun.runSim(2,hf=c(10,300),kds=rep(.8,2))
fun.runSim(2,hf=c(10,300),kds=rep(5,2))
fun.runSim(3,hf=c(10,300),kds=rep(.2,2))
fun.runSim(3,hf=c(10,300),kds=rep(.8,2))
fun.runSim(3,hf=c(10,300),kds=rep(5,2))
par(mfrow=c(1,1))
dev.off()


# hf-check ----------------------------------------------------------------
hfs <- read.csv('./v4-hf-final.csv',header = T,row.names = 1,stringsAsFactors = F)
hfs['Fos',]
hfs["Ccl5",]

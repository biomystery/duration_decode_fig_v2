loadLib <- function (){
  require(deSolve)
  require(pracma)
  require(rootSolve)
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(gridExtra)
  
}
loadLib()
col.map <- c(LPS="#f8766d",TNF="#00ba38",IL1="#619cff")
# functions ---------------------------------------------------------------
## subroutine: The model
odeModel <- function (Time, State, Pars) {
  
  # ODES
  with(as.list(c(State, Pars)),{
    # aux
    O = 1-C-A  
    if(class(nfkb_input)=="data.frame"){
      nfkb <-pchipfun(nfkb_input$Time,nfkb_input$Value)(Time) # nodelay 
      fp <-   ifelse(Time < TS, kp *A , 0) 
    }else{
      nfkb <- nfkb_input
      fp <- kp * A
    }
    
    # fluxes 
    f1 <- k1*nfkb^n1/(nfkb^n1+Kd1^n1)*C
    f_1 <- k_1 * O 
    f2 <- k2*nfkb^n2/(nfkb^n2+Kd2^n2)*O+k2_0*O
    f_2 <- k_2 * A 
    f_3 <- k_3 * A^n3/(A^n3+Kd3^n3)
    
    fdeg <- kdeg*mRNA
    
    #rhe
    dC <- -f1 + f_1 + f_3 
    dA <- f2 - f_2 - f_3  # f1 -f2 + f_2 -f_1
    dmRNA    <- fp - fdeg
    
    # return 
    return(list(c(dC, dA, dmRNA)))
  })
}

runModel<- function (pars,times=manytp){
  ## init to steady state 
  pars.2 <- pars;pars.2$nfkb_input <- pars$nfkb_input$Value[1] 
  ss.init <- runsteady(y = c(C=1,A=0,mRNA=0), #init
            fun = odeModel, 
            parms = pars.2,time = c(0,1e5))$y
  
  ## solve ode by deSolve 
  out   <- ode(ss.init, times, odeModel, pars)
  out<-as.data.frame(out)
  out$mRNA[nrow(out)] <- out$mRNA[nrow(out)-1]
  out
}
# ss scan -----------------------------------------------------------------
## steady state check 
# with or wo shutdown 
pars<- list(k1=1,k2=1,k_1=log(2)/6,k_2=log(2)/480,k_3=0,
            Kd1=.5,Kd2=.5,Kd3=.5,n1=3,n2=3,n3=3,kp=1,kt=1,
            kdeg=log(2)/30,nfkb_input=0.02,k2_0=0,KdA=0.5,nA=3,
            TS=120)

# simulation --------------------------------------------------------------
## input 
#input.set  <- readRDS(file='../../Fig_code/fig.2_knockout/subfig2A.rds')
#ggplot_build(input.set)$data[[2]]
ev <- new.env()
load(file='../data/mRNA-Fit-avg-sp-v1c-pre.Rdata',ev)
nfkb.input.set <- rbind(data.frame(ev$input.emsa.ifnar%>% gather(key="Stimuli",value = "Value",2:4),
                                   Genotype='Control',stringsAsFactors = F),
                        data.frame(ev$input.emsa.ifnarikba%>% gather(key="Stimuli",value = "Value",2:4),
                                   Genotype='Mutant',stringsAsFactors = F))

## sample simulations
tps <- seq(0,480,by=0.1)
pars.b <- pars

runSim <- function(pars=pars.b) {
  res.sim<- lapply(list('Control','Mutant'), function(g)
    lapply(list("LPS","TNF"),function(s){
      pars$nfkb_input<- (nfkb.input.set%>% 
                           filter(Genotype==g&Stimuli==s) %>%
                           select(Time,Value))
      res.out <- runModel(pars,tps)[seq(1,4800,by = 20),]
      data.frame(res.out%>% gather(key = 'Species',value = 'level',2:4),
                 Genotype=g,
                 Stimuli = s)}))
  res.out <- do.call(rbind,list(res.sim[[1]][[1]],res.sim[[1]][[2]],res.sim[[2]][[1]],res.sim[[2]][[2]]))
  p.output <- ggplot(res.out,aes(time,level,colour=Stimuli)) +
    geom_line(aes(linetype=Genotype))

  tmp.pd <- res.out%>%filter(Species=='mRNA')%>%
    group_by(Genotype)%>% 
    mutate(level=level/max(level))
  tmp.pd$Species <- recode(tmp.pd$Species,mRNA='mRNA(n)')
  
  p.output<-p.output + geom_line(data = tmp.pd,
                       aes(time,level,colour=Stimuli,linetype=Genotype)) +
    facet_grid( Species~.,scales = 'free_y')  + 
    scale_color_manual(values = col.map)
  
  return(p.output)
}

if(F){
  p.input <- ggplot(nfkb.input.set%>%filter(Stimuli!='IL1')%>%mutate(Species="NFkBn"),
                    aes(Time,Value,colour=Stimuli))+
    geom_point(aes(shape=Genotype)) + geom_line(aes(linetype=Genotype))+ 
    facet_grid(Species~.)
  
  # wrapper 
  pars.b$k_3 <- .1
  pars.b$TS <- 150 ;pars.b$kdeg <- log(2)/120
  p.outputs <- lapply(c(6,60,480),function(x){
    pars.b$k_2 <- log(2)/x
    cat(x)
    runSim(pars = pars.b)
  })
  noLegend <- theme(legend.position = 'none')
  grid.arrange(p.input+noLegend,p.input+noLegend,p.input,
               p.outputs[[1]]+noLegend,
               p.outputs[[2]]+noLegend,p.outputs[[3]],
               nrow=2,heights=c(1,4),widths=c(3,3,3.5))
  
}

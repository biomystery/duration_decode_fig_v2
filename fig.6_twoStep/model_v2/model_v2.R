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
    }else{
      nfkb <- nfkb_input
    }

    # fluxes
    f1 <- k1*nfkb^n1/(nfkb^n1+Kd1^n1)*C
    f_1 <- k_1 * O
    f2 <- k2*nfkb^n2/(nfkb^n2+Kd2^n2)*O
    f_2 <- k_2 * A
    fdeg <- kdeg*mRNA
    fp <- kp * A

    #rhe
    dC <- -f1 + f_1
    dA <- f2 - f_2
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
pars<- list(k1=1,k2=1,k_1=log(2)/6,k_2=log(2)/480,
            Kd1=.5,Kd2=.5,n1=3,n2=3,kp=1,kt=1,
            kdeg=log(2)/30,nfkb_input=0.02
            )

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



# fitting -----------------------------------------------------------------
rpkm.all<- read.csv(file='../data/rpkm.all.csv',stringsAsFactors = F,
                  row.names = 1)

pd.maxRPKM <- rpkm.all %>% group_by(Genotype,gene,Species,Stimuli)%>%
  summarise(max.rpkm = max(rpkm))

fun.pltmaxRPKM <- function(pd=pd.maxRPKM,g){
  ggplot(pd%>% filter(gene==g),
         aes(Species,max.rpkm,fill=Species))+
    geom_bar(stat = 'identity')  +
    facet_wrap(~ Genotype+Stimuli,ncol = 4)}


noLegend <- theme(legend.position = 'none')

pd.tc <- rpkm.all %>% group_by(Genotype,gene,Species)%>%
  mutate(max.rpkm = max(rpkm),
         frac.exp = rpkm/max(rpkm))

pd.tc$Time <- as.numeric(pd.tc$Time)


fun.pltTC <- function(pd=pd.tc,g,dotOnly=F){
  p <- ggplot(pd.tc%>% filter(gene==g),aes(Time,frac.exp,colour=Species)) +
    geom_point()+ facet_wrap(~ Genotype+Stimuli,ncol = 4)
  if(!dotOnly)
    return(p+ geom_line())
  else
    return(p)
}

fun.obj <- function(pars.b,doTrace=F,showAll=F){
  ## need global variable - gt
  pd.sim <- runSim(pars = pars.b)$data  %>%
    filter(Species %in% c("A",'mRNA'))%>%
    group_by(Genotype,Species) %>%
    mutate(level=level/max(level))
  pd.sim$Species <- recode(pd.sim$Species,A='caRNA',mRNA='cytoRNA')
  pd.sim$Stimuli <- recode(pd.sim$Stimuli,LPS='Lps',TNF='Tnf')
  colnames(pd.sim)<-c('Time',"Species","frac.exp","Genotype","Stimuli")
  pd.sim$Time <- pd.sim$Time/60
  plt.tc <- fun.pltTC(g = gt,dotOnly = T)+noLegend
  plt.tc_sim <- plt.tc + geom_line(data = pd.sim)
  if(doTrace) print(plt.tc_sim)
  pd.tc_sim <- rbind(data.frame(pd.sim,type='Sim'),
                     data.frame(plt.tc$data[,colnames(pd.sim)],
                                type="Exp"))%>%
    spread(key=type,value=frac.exp)
  pd.tc_sim<- pd.tc_sim[complete.cases(pd.tc_sim),]%>%
    mutate(Residual= Sim-Exp)
  
  pd.weight <- pd.maxRPKM%>% filter(gene==gt)%>%
    group_by(Genotype,Species)%>%
    mutate(w=max(max.rpkm))%>%
    ungroup()%>%
    mutate(w=w/sum(w)*2)
  
  # weighted score
  if(showAll){
    return(list(score=(pd.tc_sim %>% 
                         left_join(pd.weight[,],
                                   by=c("Species","Genotype","Stimuli"))%>%
                         summarise(score= sum(w^2*Residual^2)))$score,
                plt = plt.tc_sim,
                simData = pd.sim))
  }else{
    return(list(score=(pd.tc_sim %>%
                         left_join(pd.weight[,],
                                   by=c("Species","Genotype","Stimuli"))%>%
                         summarise(score= sum(w^2*Residual^2)))$score))
    
  }
}



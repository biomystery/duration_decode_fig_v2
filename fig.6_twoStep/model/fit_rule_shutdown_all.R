rm(list=ls())
require(ggplot2)
require(dplyr)
# load all the data --------------------------------------------------------------

rpkm.all<- read.csv(file='../data/rpkm.all.csv',stringsAsFactors = F,
                    row.names = 1)
noLegend <- theme(legend.position = 'none')
# plot1_tc ----------------------------------------------------------------

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

# gene  -------------------------------------------------------------------
#gt <- 'Nfkbie'

#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

gt <- args[1]
#gt <- 'Acpp'


# plot2_maxRPKM -----------------------------------------------------------
pd.maxRPKM <- rpkm.all %>% group_by(Genotype,gene,Species,Stimuli)%>%
  summarise(max.rpkm = max(rpkm))
fun.pltmaxRPKM <- function(pd=pd.maxRPKM,g){
  ggplot(pd%>% filter(gene==g),
         aes(Species,max.rpkm,fill=Species))+
    geom_bar(stat = 'identity')  + 
    facet_wrap(~ Genotype+Stimuli,ncol = 4)}

# objective function -----------------------------------------------------------
source('./model_rule_shutdown.R')
fun.obj <- function(pars.b,doTrace=F){
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
  return(list(score=(pd.tc_sim %>% 
            left_join(pd.weight[,],
                      by=c("Species","Genotype","Stimuli"))%>%
            summarise(score= sum(w^2*Residual^2)))$score))
}
# test
if(T){ # TEST 
  system.time(fun.obj(pars.b = pars,doTrace = T))
  #   user  system elapsed 
  #6.020   0.078   6.174 
  system.time(a<-fun.obj(pars.b = pars,doTrace = F))
  #user  system elapsed 
  #5.412   0.045   5.496   
}

pars$TS <- 60*8 
fun.obj_wapper.1 <- function(pars.fit,plotFlag=F){
  pars.2 <- pars
  #pars<- list(k1=1,k2=1,k_1=log(2)/6,k_2=log(2)/480,k_3=0,
  #            Kd1=.5,Kd2=.5,Kd3=.5,n1=3,n2=3,n3=3,kp=1,kt=1,
  #            kdeg=log(2)/30,nfkb_input=0.02,k2_0=0,KdA=0.5,nA=3,
  #            TS=120)
  pars.2[c('k2','k_1','k_2','k_3',"Kd1","Kd2","Kd3",
           'kdeg')] <- 10^pars.fit
  
  
  a<- fun.obj(pars.b = pars.2)
  
  # update BEST_SCORE
  if(a$score<BEST_SCORE[length(BEST_SCORE)]) {
    BEST_SCORE <<- c(BEST_SCORE,a$score)}else{
    BEST_SCORE <<- c(BEST_SCORE,BEST_SCORE[length(BEST_SCORE)])}
  
  # plt.bestScore 
  NUM_ITER <<- NUM_ITER + 1
  cat('iter:',NUM_ITER,'pars:',unlist(pars.2),'\n')
  cat('score:',signif(a$score,2),'Best_score',signif(BEST_SCORE[length(BEST_SCORE)],2))
  
  # plots
  if(plotFlag){
    plt.bestScore <- ggplot(data.frame(Iter=1:NUM_ITER,
                                       BestScore = BEST_SCORE[-1]),
                            aes(Iter,BestScore))+
      geom_point()+geom_line()
    
    # plt.pars
    #require(ggradar)
    plt.pars <- ggplot(data.frame(par=names(pars.2)[-15],
                                  value=unlist(pars.2)[-15])%>%
                         filter(value!=0),
                       aes(par,log10(value)))+
      geom_bar(stat = 'identity')
    
    # plot
    grid.arrange(grobs=list(a$plt,plt.bestScore,plt.pars),
                 layout_matrix=rbind(c(1,1),
                                     c(2,3)))
  }
  
  # time
  cat(" t:",proc.time()[3]-PTM[3],'secs\n')
  return(a$score)
}

# optimization ------------------------------------------------------------
kdegs <- readRDS(file='../data/genes.kdeg.Rds')
resDir <- "./hill_1"
NUM_ITER <- 0;BEST_SCORE<- 9999

#pars.2[c('k2','k3','k_1','k_2','k_3',"Kd1","Kd2","Kd3",
#'kdeg','TS')] <- 10^pars.fit
pars.fit <- c(k2=0,k_1=-1,k_2=-1,k_3=-1,
              Kd1=-1,Kd2=-1,Kd3=-1,
              kdeg=log10(kdegs$val[[gt]]))
              #TS = log10(120));
lb <- c(rep(-3,7),log10(kdegs$lb[[gt]]))#,-3);
ub<-c(rep(3,7),log10(kdegs$ub[[gt]]))#,log10(481))
write.csv(file = paste0(resDir,'fit_setting.csv'),rbind(pars.fit,lb,ub))

# GSA
if(T){
  require(GenSA)
  set.seed(1234)
  PTM <- proc.time()
  try(fit <- GenSA(fn = fun.obj_wapper.1,
                   lower=lb,upper =ub,
                   control=list(verbose=T,
                                nb.stop.improvement=50,
                                max.call = 1000)))
  save(fit,file = paste0(resDir,gt,'.Rdata'))
}

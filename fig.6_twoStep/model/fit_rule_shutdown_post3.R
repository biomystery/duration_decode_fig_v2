rm(list=ls())
# load all the data --------------------------------------------------------------
require(ggplot2)
require(dplyr)
require(tidyr)
rpkm.all<- read.csv(file='./data/rpkm.all.csv',stringsAsFactors = F,
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

# plot2_maxRPKM -----------------------------------------------------------
pd.maxRPKM <- rpkm.all %>% group_by(Genotype,gene,Species,Stimuli)%>%
  summarise(max.rpkm = max(rpkm))
fun.pltmaxRPKM <- function(pd=pd.maxRPKM,g){
  ggplot(pd%>% filter(gene==g),
         aes(Species,max.rpkm,fill=Species))+
    geom_bar(stat = 'identity')  + 
    facet_wrap(~ Genotype+Stimuli,ncol = 4)}

# plot3_sp ----------------------------------------------------------------
glist <- read.csv(file='./data/glist_final.csv',header = F,stringsAsFactors = F)
pd.sp <- pd.tc %>%
  summarise(Sp = log2(max(frac.exp[Stimuli=="Lps"]/max(frac.exp[Stimuli=='Tnf']))))%>%
  spread(key=Genotype,value=Sp)%>%
  filter(gene %in% glist$V1)

sp.th <- 2^0.5 # or 2
pd.sp <- pd.sp%>%
  mutate(cate=ifelse(Control<0.5,"I",
                     ifelse(Mutant<Control/sp.th,"III","II")))
table((pd.sp%>%filter(Species=='cytoRNA'))$cate)

pd.sp$Species<- factor(pd.sp$Species,
                       levels = c('cytoRNA','caRNA'))
plt.sp.all<- ggplot(data.frame(Control=c(0.5,0.5,8,8,0.5,0.5,8,8,-3,-3,0.5,0.5),
                               Mutant=c(0.5/sp.th,8,8,8/sp.th, 
                                        -3,0.5/sp.th,8/sp.th,-3,
                                        -3,8,8,-3),
                               Class=factor(rep(c(2,3,1),each=4))),
                    aes(Control,Mutant))+  
  coord_cartesian(xlim=range(pd.sp[,3:4],na.rm = T),ylim=range(pd.sp[,3:4],na.rm = T)) +
  geom_polygon(aes(fill=Class)) + scale_fill_manual(values = alpha(c("blue"), c(.05,.1,.2)))+
  theme_bw()+
  geom_line(data = pd.sp,aes(group=gene,label=gene),colour=grey(.8))+
  geom_point(data = pd.sp,aes(Control,Mutant,colour=Species,label=gene),alpha=0.4)+
  geom_hline(yintercept = c(-.5,.5),colour=grey(.8),linetype=2)+
  geom_vline(xintercept = c(-.5,.5),colour=grey(.8),linetype=2)+
  geom_abline(slope = 1,intercept = 0,colour=grey(.8),linetype=2)+
  scale_color_brewer(palette = 'Set1')+
  theme(legend.position = 'none')

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
            summarise(score= sum(w^2*Residual^2)))$score,
            plt = plt.tc_sim,
            simData = pd.sim))
}

pars$TS <- 60*8 



# plot_best_fit -----------------------------------------------------------
resDir <- './3rd/'
res.files<- list.files(path = resDir,pattern = "*.Rdata")
gs <- sub(".Rdata","",res.files)
i <-0 
fun.pltBestFit <- function(res.file){
  fileLoc<- paste0(resDir,res.file)
  ev<- new.env();load(file = fileLoc,ev)  
  pars.2 <- pars
  pars.2[c('k2','k_1','k_2','k_3',"Kd1","Kd2","Kd3",
           'kdeg')] <- 10^ev$fit$par
  gt <<- sub(".Rdata","",res.file)
  a<- fun.obj(pars.b = pars.2)
  tls<-paste0(gt,":",signif(ev$fit$value,2)) 
  print(i<<- i+1)
  return(list(plt= a$plt+ggtitle(tls),
              pd = a$simData,
              fitRes = ev$fit))
}


plt.bestFits <- lapply(res.files,fun.pltBestFit)
names(plt.bestFits) <- gs
saveRDS(file=paste0(resDir,'allRes.Rds'),plt.bestFits)


# plot_wrapper ------------------------------------------------------------
noLegend <- theme(legend.position = 'none')
fun.pltWrapper<- function(gt){
  plt.tc <- plt.bestFits[[gt]]$plt+noLegend+
    scale_color_brewer(palette = 'Set1')
  pd.sp_sim<-plt.bestFits[[gt]]$pd %>% 
    group_by(Species,Genotype,Stimuli)%>%
    summarise(frac.max = max(frac.exp))%>%
    ungroup()%>%
    group_by(Genotype,Species)%>%
    summarise(sp=log2(frac.max[1]/frac.max[2]))%>%
    spread(key = Genotype,sp)
  
  plt.maxRPKM <-fun.pltmaxRPKM(g=gt)+noLegend+
    scale_fill_brewer(palette = 'Set1')
  
  plt.sp<- plt.sp.all + 
    geom_point(data = pd.sp%>%filter(gene==gt),
               aes(Control,Mutant,colour=Species),size=4)+
    geom_line(data = pd.sp%>%filter(gene==gt),aes(group=gene),colour='black')+
    ggtitle(gt)
  plt.sp<-plt.sp + geom_point(data=pd.sp_sim,
                      aes(colour=Species),size=4,shape=17)+
    geom_line(data=pd.sp_sim,colour='blue')
  grid.arrange(grobs=list(plt.tc,plt.maxRPKM,plt.sp),
               layout_matrix=rbind(c(1,3),
                                   c(2,3)))
}

if(T){
  i <- 0 
  pdf(file=paste0(resDir,'caRNAvcytoRNA_bestFit.pdf'),width = 12)
  sapply(gs,function(x) {
    print(fun.pltWrapper(x))
    i <<- i+1
    print(i)
  })
  dev.off()
}


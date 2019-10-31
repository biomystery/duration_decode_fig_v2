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
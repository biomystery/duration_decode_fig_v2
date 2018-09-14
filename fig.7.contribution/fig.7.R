require(tidyverse)
require(ggplot2)
fig_dir <- "~/Dropbox/Projects/DurationDecoding/figure/Fig.7/subfigs/"
# contribution (simulation)  -----------------------------------------------------------
# load simulation data 
sim.data <- readRDS(file = "../fig.6_twoStep/pd.fig6.Rds")
names(sim.data)
attach(sim.data)
dim(sim_ca_mat)
rownames(sim_ca_mat); colnames(sim_ca_mat)

calSp <- function(simData=sim_ca_mat){
  simData %>% add_column(gene=rownames(simData))%>%
    gather(key = "conditions",value = "expression",1:ncol(simData))%>%
    separate(conditions,into = c("geno","sti","time"),sep = "_")  %>%
    group_by(gene,geno,sti)%>% summarise(m=max(expression))%>%
    summarise(sp=m[1]/m[2])
}

sp.ca <- calSp()
sp.mRNA <- calSp(sim_cyto_mat)

pd <- left_join(sp.ca,sp.mRNA,by=c("gene","geno"),suffix=c(".ca",".cyto")) %>%
  summarise(ca.contri = (sp.ca[1]-sp.ca[2])/(sp.cyto[1]-sp.cyto[2]))%>%
  mutate(cyto.contri=1-ca.contri)%>% arrange(ca.contri)%>%
  gather(key = "Mechanism",value = "Contribution",2:3)

#pd$Contribution[pd$Contribution>1] = 1; pd$Contribution[pd$Contribution<0] = 0

pd$gene <- factor(pd$gene,levels = pd$gene[1:(nrow(pd)/2)])
pd$Mechanism <- factor(pd$Mechanism,levels = c("cyto.contri","ca.contri"))
recode_factor(pd$Mechanism,cyto.contri = "Chr.",ca.contri="T1/2")
ggplot(pd,aes(gene,Contribution))+
  geom_bar(aes(fill=Mechanism),stat = 'identity')+
  coord_flip()+scale_fill_brewer(palette = 'Set1')+
  theme(axis.title.y = element_blank())+ geom_hline(yintercept =  .5,colour="black",linetype=2)

ggsave(filename = paste0(fig_dir,"subfigA.eps"),width = 8,height = 12,scale = .75)
# contribution (data)  -----------------------------------------------------------
pd.sp.contri<- pd.sp %>% 
  group_by(gene)%>%
  summarize(ca.contri=(Control[1]-Mutant[1])/(Control[2]-Mutant[2]))
pd.sp.contri<-pd.sp.contri%>%arrange(ca.contri)
pd.sp.contri.2 <-pd.sp.contri%>%
  mutate(cyto.contri = 1-ca.contri)%>%
  gather(key=Species,Contribution,2:3)
pd.sp.contri.2$gene <- factor(pd.sp.contri.2$gene,
                              levels = pd.sp.contri$gene)
pd.sp.contri.2$Species<- factor(pd.sp.contri.2$Species,
                                levels = c("cyto.contri","ca.contri"))
ggplot(pd.sp.contri.2,aes(gene,Contribution))+
  geom_bar(aes(fill=Species),stat = 'identity')+
  coord_flip()+scale_fill_brewer(palette = 'Set1')
pd.sp%>%filter(gene=='Cx3cl1')
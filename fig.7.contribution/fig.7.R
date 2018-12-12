require(tidyverse)
require(ggplot2)
subfig_dir <- "../figures/Fig.7/subfigs/"
# contribution (simulation)  -----------------------------------------------------------
# load simulation data 
sim.data <- readRDS(file = "../fig.6_twoStep/pd.fig6.Rds")
scores.df.final <- read.csv("../fig.6_twoStep/data/final.score.csv",stringsAsFactors = F)
v2.genes <- scores.df.final%>% filter(minScore<0.13)%>%
  select(gene)
names(sim.data)
attach(sim.data)
dim(sim_ca_mat)
rownames(sim_ca_mat); colnames(sim_ca_mat)

calSp <- function(simData=sim_ca_mat){
  simData %>% add_column(gene=rownames(simData))%>%
    gather(key = "conditions",value = "expression",1:ncol(simData))%>%
    separate(conditions,into = c("geno","sti","time"),sep = "_")  %>%
    group_by(gene,geno,sti)%>% summarise(m=max(expression))%>%
    summarise(sp=log2(m[1]/m[2]))
}

sp.ca <- calSp()
sp.mRNA <- calSp(sim_cyto_mat)

pd <- left_join(sp.ca,sp.mRNA,by=c("gene","geno"),suffix=c(".ca",".cyto")) %>%
  summarise(ca.contri = (sp.ca[1]-sp.ca[2])/(sp.cyto[1]-sp.cyto[2]))%>%
  mutate(cyto.contri=1-ca.contri)%>% arrange(ca.contri)%>%
  gather(key = "Mechanism",value = "Contribution",2:3)

pd$Contribution[pd$Contribution>1] = 1; pd$Contribution[pd$Contribution<0] = 0

pd$gene <- factor(pd$gene,levels = pd$gene[1:(nrow(pd)/2)])
pd$Mechanism <- factor(pd$Mechanism,levels = c("cyto.contri","ca.contri"))
recode_factor(pd$Mechanism,cyto.contri = "Chr.",ca.contri="T1/2")
ggplot(pd%>%filter(gene %in% v2.genes$gene),aes(gene,Contribution))+
  geom_bar(aes(fill=Mechanism),stat = 'identity')+
  coord_flip()+scale_fill_brewer(palette = 'Set1')+
  theme(axis.title = element_blank(),legend.position = "none")+ 
  geom_hline(yintercept =  .5,colour="black",linetype=2)

ggsave(filename = paste0(subfig_dir,"subfigA.eps"),
       width = 8,height = 12,scale = .75)

# contribution (data)  -----------------------------------------------------------
pd.sp <- readRDS('../fig.5_caRNA/data/pd.merge.sp.Rds')
dat.sp.pd <-pd.sp%>%mutate(ca.delta=ca.ctrl.sp-ca.mt.sp,
               cyto.delta=cyto.ctrl.sp-cyto.mt.sp)%>%
  select(-contains("ctrl"))%>%
  gather(key = 'key',
         value="sp",c(1,2,4,5))%>%
  mutate(key=sub(".sp","",key))%>%
  separate(key,into = c("type","geno"),sep = "[.]")


# new version -------------------------------------------------------------
pd.new <- left_join(sp.ca,sp.mRNA,by=c("gene","geno"),suffix=c(".ca",".cyto"))
sp.mRNA.pd <- sp.mRNA%>%
  group_by(gene)%>%
  summarise(Delta=sp[1]-sp[2],
            Mutant=sp[2])%>%
  gather(key = "geno",value = "sp",2:3)
sp.ca.pd <- sp.ca%>%
  group_by(gene)%>%
  summarise(Delta=sp[1]-sp[2],
            Mutant=sp[2])%>%
  gather(key = "geno",value = "sp",2:3)

sp.pd <- rbind(data.frame(sp.mRNA.pd,type='cyto'),
               data.frame(sp.ca.pd,type='ca'))
  
ggplot(sp.pd,aes(gene,sp)) +
  geom_bar(stat = 'identity',aes(fill=geno))+
  coord_flip()+
  facet_wrap(~type)+
  geom_hline(yintercept = 0.5)

## merged sim and dat
pd.all <- bind_rows(dat.sp.pd%>%
  mutate(type=paste0(type,".dat"),
         geno= recode(geno,delta="Delta",mt="Mutant")),
  sp.pd%>%
    mutate(type=paste0(type,".sim")))%>%
  mutate(type=factor(type,levels = c("ca.dat","cyto.dat","ca.sim","cyto.sim")))

ggplot(pd.all,aes(gene,sp)) +
  geom_bar(stat = 'identity',aes(fill=geno))+
  coord_flip()+
  facet_wrap(~type,nrow = 1)+
  geom_hline(yintercept = 0.5)


# Merge the contribution  -------------------------------------------------
pd.all.2 <- pd.all%>%
  filter(type%in%c("ca.sim","cyto.sim"))%>%
  mutate(sp=ifelse(sp<0,0,sp))
  
  
ggplot(pd.all.2,aes(gene,sp)) +
  geom_bar(stat = 'identity',aes(fill=geno))+
  coord_flip()+
  facet_wrap(~type,nrow = 1)+
  geom_hline(yintercept = 0.5)

pd.all.3<- pd.all.2%>%
  group_by(gene)%>%
  unite(col = "geno_type",2:3)%>%
  spread(key = "geno_type",value = "sp")%>%
  mutate(cyto.sim_DegContri=cyto.sim_Delta-ca.sim_Delta)%>%
  select(-c(3,4))%>%
  mutate(cyto.sim_DegContri=ifelse(cyto.sim_DegContri<0,0,cyto.sim_DegContri))%>%
  gather(key = "Specificity",value = "SP",2:4)%>%
  mutate(Specificity=recode(Specificity,'ca.sim_Delta'="cyto.sim_CaContri",
                            'cyto.sim_Mutant'="cyto.sim_Remain"))%>%
  mutate(Specificity=factor(Specificity,levels = rev(c("cyto.sim_Remain",
                                                   "cyto.sim_CaContri",
                                                   "cyto.sim_DegContri"))))
  

ggplot(pd.all.3,aes(gene,SP)) +
  geom_bar(stat = 'identity',aes(fill=Specificity))+
  coord_flip(expand = F)+
  geom_hline(yintercept = 0.5)+
  scale_fill_manual(values = rev(c("grey",col.caRNA,col.cytoRNA)))
  
saveRDS(pd.all.3,file='fig7B.rds')
  
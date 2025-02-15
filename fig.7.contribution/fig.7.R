rm(list = ls(all.names = T))
source('../auxilary_functions.R')
require(data.table)
require(tidyverse)
require(ggplot2)
subfig_dir <- "../figures/Fig.7/subfigs/"
# contribution (simulation)  -----------------------------------------------------------
# load simulation data 
sim.data <- readRDS(file = "../fig.6_twoStep/pd.fig6.Rds")
v2.genes <- fread("../fig.6_twoStep/data/model_v2_pars.csv")%>%
  filter(score<0.13)%>%
  dplyr::select(gene)
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


# contribution (data)  -----------------------------------------------------------
pd.sp <- readRDS('../fig.5_caRNA/data/pd.merge.sp.Rds')
dat.sp.pd <-pd.sp%>%mutate(ca.delta=ca.ctrl.sp-ca.mt.sp,
               cyto.delta=cyto.ctrl.sp-cyto.mt.sp)%>%
  dplyr::select(-contains("ctrl"))%>%
  gather(key = 'key',
         value="sp",c(1,2,4,5))%>%
  mutate(key=sub(".sp","",key))%>%
  separate(key,into = c("type","geno"),sep = "[.]")


# new version -------------------------------------------------------------
sp.ca <- calSp()
sp.mRNA <- calSp(sim_cyto_mat)

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
if(T){
  pd.all.2 <- pd.all%>%
    filter(type%in%c("ca.sim","cyto.sim"),gene%in%v2.genes$gene)%>%
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
    dplyr::select(-c(3,4))%>%
    mutate(cyto.sim_DegContri=ifelse(cyto.sim_DegContri<0,0,cyto.sim_DegContri))%>%
    gather(key = "Specificity",value = "SP",2:4)%>%
    mutate(Specificity=recode(Specificity,'ca.sim_Delta'="cyto.sim_CaContri",
                              'cyto.sim_Mutant'="cyto.sim_Remain"))%>%
    mutate(Specificity=factor(Specificity,levels = rev(c("cyto.sim_Remain",
                                                         "cyto.sim_CaContri",
                                                         "cyto.sim_DegContri"))))
  
  
  pd.subfigB.ord <- pd.all.3%>%spread(key = Specificity,value = SP)%>%
    mutate(group=ifelse(cyto.sim_DegContri<0.15,1,
                        ifelse(cyto.sim_CaContri<.15,3,2)),
           sp=cyto.sim_DegContri+cyto.sim_CaContri+cyto.sim_Remain)%>%
    arrange(desc(group),(sp))%>%pull(gene)
  
  ggplot(pd.all.3%>%ungroup%>%mutate(gene=factor(gene,levels = pd.subfigB.ord)),aes(gene,SP)) +
    geom_bar(stat = 'identity',aes(fill=Specificity))+
    coord_flip(expand = F)+
    geom_hline(yintercept = 0.5)+  theme_bw()+theme(legend.title =element_blank(),
                                                    legend.position = c(.8,.2),legend.justification=c(1,0) )+
    scale_fill_manual(values = rev(c("grey",col.caRNA,col.cytoRNA)))
  
  
  ggsave(filename = paste0(subfig_dir,'Fig7B.pdf'),
         width = 4,height = 8,units = 'in',scale=1.5);system(paste0("open ",subfig_dir,'Fig7B.pdf'))
  

  ## subfigC 
  require(ggrepel)
  pd.subfigC <- pd.all.3%>%spread(key = Specificity,value = SP)%>%
    mutate(cate=ifelse(gene %in% c('Fpr1','Mmp3','Ccl5','Pilra'),'ca',
                       ifelse(gene %in% c('Nfkb2','Cebpb','Rab15','Rab20'),'cyto','none')))
  
  ggplot(pd.subfigC,aes( cyto.sim_DegContri,cyto.sim_CaContri)) +
    geom_point()+geom_point(data=pd.subfigC%>%filter(cate!='none'),aes(color=cate),size=4)+
    scale_color_manual(values = c('ca'=col.caRNA,'cyto'=col.cytoRNA))+
    geom_text_repel(data=pd.subfigC%>%filter(cate!='none'),aes(label=gene),point.padding = NA)+
    theme_bw()+theme(legend.position = 'none')+ coord_fixed()
  
  ggsave(filename = paste0(subfig_dir,'Fig7C.pdf'),useDingbats=F,
         width = 2,height = 4,units = 'in',scale = 2);system(paste0('open ',subfig_dir,'Fig7C.pdf'))
  
  ggplot(pd.subfigC,aes( cyto.sim_DegContri,cyto.sim_CaContri)) +
    geom_point()+
    scale_color_manual(values = c('ca'=col.caRNA,'cyto'=col.cytoRNA))+
    geom_text_repel(aes(label=gene),point.padding = NA)+
    theme_bw()+theme(legend.position = 'none')+ coord_fixed()
  ggsave(filename = paste0(subfig_dir,'Fig7C_label.pdf'),useDingbats=F,
         width = 2,height = 4,units = 'in',scale = 2);system(paste0('open ',subfig_dir,'Fig7C_label.pdf'))
  
}


# fig7BC_data -------------------------------------------------------------

if(T){
  pd.all.2 <- pd.all%>%
    filter(type%in%c("ca.dat","cyto.dat"),gene%in%v2.genes$gene)%>%
    mutate(sp=ifelse(sp<0,0,sp))
  
 
  pd.all.3<- pd.all.2%>%
    group_by(gene)%>%
    unite(col = "geno_type",2:3)%>%
    spread(key = "geno_type",value = "sp")%>%
    mutate(cyto.dat_DegContri=cyto.dat_Delta-ca.dat_Delta)%>%
    dplyr::select(-c(3,4))%>%
    mutate(cyto.dat_DegContri=ifelse(cyto.dat_DegContri<0,0,cyto.dat_DegContri))%>%
    gather(key = "Specificity",value = "SP",2:4)%>%
    mutate(Specificity=recode(Specificity,'ca.dat_Delta'="cyto.dat_CaContri",
                              'cyto.dat_Mutant'="cyto.dat_Remain"))%>%
    mutate(Specificity=factor(Specificity,levels = rev(c("cyto.dat_Remain",
                                                         "cyto.dat_CaContri",
                                                         "cyto.dat_DegContri"))))
  
  
  pd.subfigB.ord <- pd.all.3%>%spread(key = Specificity,value = SP)%>%
    mutate(group=ifelse(cyto.dat_DegContri<0.2,1,
                        ifelse(cyto.dat_CaContri<.15,3,2)),
           sp=cyto.dat_DegContri+cyto.dat_CaContri+cyto.dat_Remain)%>%
    arrange(desc(group),(sp))%>%pull(gene)
  

    
  ggplot(pd.all.3%>%ungroup%>%mutate(gene=factor(gene,levels = pd.subfigB.ord)),aes(gene,SP)) +
    geom_bar(stat = 'identity',aes(fill=Specificity))+
    coord_flip(expand = F)+
    geom_hline(yintercept = 0.5)+  theme_bw()+theme(legend.title =element_blank(),
                                                    legend.position = c(.8,.2),legend.justification=c(1,0) )+
    scale_fill_manual(values = rev(c("grey",col.caRNA,col.cytoRNA)))
  
  
  ggsave(filename = paste0(subfig_dir,'subfig7B_dat.pdf'),
         width = 4,height = 8,units = 'in',scale=1.5);  system(paste0('open ',subfig_dir,'subfig7B_dat.pdf'))
  
  ## subfigC 
  require(ggrepel)
  pd.subfigC <- pd.all.3%>%spread(key = Specificity,value = SP)%>%
    mutate(cate=ifelse(gene %in% c('Fpr1','Mmp3','Ccl5','Pilra'),'ca',
                       ifelse(gene %in% c('Nfkb2','Cebpb','Rab15','Rab20'),'cyto','none')))
  
  ggplot(pd.subfigC,aes( cyto.dat_DegContri,cyto.dat_CaContri)) +
    geom_point()+geom_point(data=pd.subfigC%>%filter(cate!='none'),aes(color=cate),size=4)+
    scale_color_manual(values = c('ca'=col.caRNA,'cyto'=col.cytoRNA))+
    geom_text_repel(data=pd.subfigC%>%filter(cate!='none'),aes(label=gene),point.padding = NA)+
    theme_bw()+theme(legend.position = 'none')+ coord_fixed()
  
  ggsave(filename = paste0(subfig_dir,'subfig7C_dat.pdf'),useDingbats=F,
         width = 2,height = 4,units = 'in',scale = 2)
  
  system(paste0('open ',subfig_dir,'subfig7C_dat.pdf'))
  
}


#  geom_vline(xintercept = .15)+  geom_hline(yintercept = .15)


saveRDS(pd.all.3,file='fig7B.rds')
  
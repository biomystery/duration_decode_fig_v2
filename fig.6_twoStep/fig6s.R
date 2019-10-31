source("../auxilary_functions.R")
subfig_dir <- "../figures/Fig.6/subfigs/"
fig_dir <- "../figures/Fig.6/"

# load data ---------------------------------------------------------------

rpkm.all<- read.csv(file='./data/rpkm.all.csv',stringsAsFactors = F,
                    row.names = 1)
# pd.tc 
pd.tc <- rpkm.all %>% group_by(Genotype,gene,Species)%>%
  mutate(max.rpkm = max(rpkm),
         frac.exp = rpkm/max(rpkm))
pd.tc$Time <- as.numeric(pd.tc$Time)



# score compare -----------------------------------------------------------

getScore <- function(f){
  fit <- readRDS(paste0(f,'allRes.Rds'))
  data.frame(gene=unlist(lapply(fit,function(x) x$plt$data$gene%>%unique)),
             minScore=unlist(lapply(fit,function(x) x$fitRes$value)))
             
}
scores.df.final <- getScore('./data/v1_')%>%rename('minScore'='Hill_3')%>%
  right_join(getScore('./data/v1_hill_1_')%>%rename('minScore'='Hill_1'))
    
fwrite(scores.df.final,'./data/final.score.cmp.csv')
rownames(scores.df.final)<- scores.df.final$gene
#scores.df.final <- read.csv("./data/final.score.csv",stringsAsFactors = F)
#rownames(scores.df.final)<- scores.df.final$gene
sum(scores.df.final$minScore<.15) #63


ggplot(scores.df.final%>%gather(key = 'model',value = 'score',2:3),aes(model,score))+
  geom_point(position = position_jitter(width = .1,height = 0))+geom_violin(fill=NA)+
  geom_boxplot(width=.25,fill=NA)+geom_hline(yintercept = .13,linetype=2,color='gold')+
  annotate(geom = 'text',x=1,y = .4,
           label=paste0("#gene_pass:\nHill=3 : ",sum(scores.df.final$Hill_3<0.13),
                        '\nHill=1 : ',sum(scores.df.final$Hill_1<0.13)),color='gold')

ggsave(filename = paste0(subfig_dir,'Fig6s_cmp_hill.pdf'),useDingbats=F,
       width = 6,height = 8,units = 'in'); system(paste0("open ",subfig_dir,'Fig6s_cmp_hill.pdf'))



# prepare table -----------------------------------------------------------

# prepare data (exp + best fit results)
all_res <- readRDS(file = 'data/v1_hill_1_allRes.Rds')


getPars <- function(g='Ccl5'){
  
  res <-   all_res[[g]]
  res.par <- 10^res$fitRes$par

  
  # base pars
  # with or wo shutdown
  pars<- list(k1=1,k2=1,k_1=log(2)/6,k_2=log(2)/480,
              Kd1=.5,Kd2=.5,n1=1,n2=1,kp=1,kt=1,
              kdeg=log(2)/30,nfkb_input=0.02
  )
    pars[c('k_1','k_2',"Kd1","Kd2",'kdeg')] <- res.par
  
  return(pars)
}

fig5.setting <- read_rds(path = "../fig.5_caRNA/data/pd.fig5.Rds")

pars <- do.call(rbind,lapply(fig5.setting$row_names,getPars))%>%as.data.frame
rownames(pars) <- fig5.setting$row_names

pars.final <- scores.df.final%>%dplyr::select(-Hill_3)%>%
  left_join(pars%>%rownames_to_column('gene'))%>%
  rename('Hill_1'='score')
fwrite(pars.final,file = "./data/model_v2_pars.full.csv");system("open ./data/model_v2_pars.full.csv")
fwrite(pars.final%>%
         dplyr::select(gene,score,k_1,k_2,Kd1,Kd2,kdeg),file = "./data/model_v2_pars.csv");system("open ./data/model_v2_pars.csv")

# B - example bestfit -----------------------------------------------------
plt_bestFit <- function(g="Slfn2"){
  # require 1) all_res 2) scores.df.final to tell which model 
  # got the plot 
  whichPlot <- all_res[[g]]
  whichPlot$plt$data$Species<- factor(whichPlot$plt$data$Species,levels = c("cytoRNA","caRNA"))
  whichPlot$plt + theme_bw()+ theme(legend.position = 'none',
                                    strip.text.x = element_blank(),
                                    plot.title = element_blank())+
    scale_color_brewer(palette = 'Set1') + 
    facet_grid(gene~Genotype+Stimuli)
}
p1 <- plt_bestFit(g="Ccl5")+theme(axis.text= element_blank(),title = element_blank())
p2 <- plt_bestFit(g="Cgn")+theme(axis.text = element_blank(),title = element_blank())
p3 <- plt_bestFit(g="Fpr1") +theme(axis.text = element_blank(),title = element_blank())

sapply(c('Ccl5','Cgn','Fpr1'), function(g) pars.final%>%filter(gene==g)%>%pull(score))

ggsave(filename = paste0(subfig_dir,'fig6B.pdf'),
  grid.arrange(grobs=list(p1,p2,p3)),
  width = 6,height = 4,units = 'in');system(paste0("open ",subfig_dir,'fig6B.pdf'))

if(F){
  pdf(file = paste0(fig_dir,'fig6s.pdf'),width = 6,height = 4)
  sapply(1:27, function(x){
    gs<- rownames(pd.scores)[((x-1)*3+1):(x*3)]
    p1 <- plt_bestFit(g=gs[1])+
      theme(axis.text= element_blank(),title = element_blank())
    p2 <- plt_bestFit(g=gs[2])+
      theme(axis.text = element_blank(),title = element_blank())
    p3 <- plt_bestFit(g=gs[3]) +
      theme(axis.text = element_blank(),title = element_blank())
    grid.arrange(grobs=list(p1,p2,p3))
  }
  )
  dev.off()
}


# C - heatmap of best Fit  ------------------------------------------------
# prepare heatmap (simulation) data matrix 

# 1. got the order of the genes from previous figure 
fig5.setting <- read_rds(path = "../fig.5_caRNA/data/pd.fig5.Rds")
fig5.setting$row_names

# 2. get the simulation matrix 
pd.fig6 <- list()
getSimMat <- function(g='Ccl5'){
  
  # get the best fit data for hm 
  
  res <-   all_res[[g]]
  getPDfromRes <- function(res){
    # duplicate the max timepoint 
    res.maxT <- res$pd%>% filter(Time==max(res$pd$Time))
    res.maxT$Time <- 8 
    
    res$pd <- rbind(res$pd,res.maxT)
    
    # time points 
    pd <- res$pd %>% unite(s_t,c("Species","Time"),remove=F) %>% 
      filter(s_t %in% (res$plt$data %>% unite(s_t,c("Species","Time")))$s_t)%>%
      mutate(type="Sim") %>% dplyr::select(-s_t)
    pd.exp <- res$plt$data%>% mutate(type="exp")
    return(pd%>% bind_rows(pd.exp[,colnames(pd)]))
  }
  getPDfromRes(res)%>% mutate(gene=g)
}

if(T){
  pd  <- lapply(fig5.setting$row_names,getSimMat)
  pd <- do.call(rbind,pd)
  
  # plot sim Hm 
  pd$gene <- factor(pd$gene,levels = fig5.setting$row_names)
  plotSimHm <- function(p_species="caRNA"){
    pd.ca.sim.mat <- pd%>% filter(type=="Sim"& Species==p_species) %>% ungroup()%>%
      dplyr::select(-c(Species,type))%>%
      unite(sample,c("Genotype","Stimuli","Time"))%>%
      spread(key = sample,value = frac.exp)
    pd.ca.sim.mat <- as.data.frame(pd.ca.sim.mat)
    rownames(pd.ca.sim.mat)<- pd.ca.sim.mat$gene;pd.ca.sim.mat$gene <- NULL 
    for(i in 1:4){
      setEPS()
      #wid <- ifelse(i==4,6,4);
      wid <- 2; lg <- ifelse(i==4 & p_species == "cytoRNA",T,F);
      p_colors <- ifelse(p_species=="caRNA","Blues","Reds")
      n_cols <- ifelse(p_species=="caRNA",6,7)
      ppd <- pd.ca.sim.mat[,((i-1)*n_cols+1):(i*n_cols)];rownames(ppd) <- rownames(pd.ca.sim.mat);
      postscript(paste0(subfig_dir,"subfig6c_sim",p_species,i,".eps"),onefile = F,
                 width = wid,height = 10)
      pheatmap(ppd,scale ="none",cluster_rows = F,cluster_cols = F,
               gaps_row = fig5.setting$row_seps-1,
               #gaps_col = seq(1,28,by = 6)-1,
               cellwidth = 12,fontsize = 10,
               show_rownames = lg,
               show_colnames = F,legend = F,
               color = colorRampPalette(brewer.pal(9,p_colors))(10))  
      dev.off()
    }
    pd.ca.sim.mat
  }
  
  pd.fig6$sim_ca_mat <- plotSimHm()
  pd.fig6$sim_cyto_mat <- plotSimHm(p_species = "cytoRNA")
  
  
  # legends 
  
  plotLegend<- function(cols,bks,fnames){
    require(RColorBrewer)
    #bks <- seq(round(min(rsums)),round(max(rsums))+1,length.out = 6)
    #cols<- colorRampPalette(c( "white", "blueviolet"))(5)
    #fnames<-'tmp.eps'
    setEPS()
    postscript(fnames,onefile = F,width = 0.1,height = .1*2*length(bks))
    par(mar=c(0.02, 0.04, 0.04,0.1))
    barplot(rep(1,length(bks)-1),width=diff(bks),space = 0,border = NA,
            col = colorRampPalette(brewer.pal(9,cols))(10),axes = F,horiz = T,
            xaxs = "i", yaxs = "i",xlim = c(0,1),
            ylim=c(0,sum(diff(bks)))
            #xlim=c(0,1)
    )
    axis(4,labels = F,tck = 0.2,at=cumsum(diff(bks)))
    box()
    dev.off()
  }
  plotLegend(cols="Reds",
             bks=seq(0,1,by = 0.1),
             fnames = paste0(subfig_dir,"subfig6c_simcytoRNA_lg.eps"))
  plotLegend(cols="Blues",
             bks=seq(0,1,by = 0.1),
             fnames = paste0(subfig_dir,"subfig6c_simcaRNA_lg.eps"))
  
  # scoreMatix 
  v1 <- readRDS(file = "../fig.4_modelv1/data/fig4.setting.rds")$hm_data[,1:2]
  pd.scores<- data.frame(  v2=(pars.final%>%column_to_rownames('gene'))[fig5.setting$row_names,'score'])
  pd.scores.v1 <- data.frame(v1[fig5.setting$row_names,"v1"])
  rownames(pd.scores.v1) <- rownames(pd.scores) <- fig5.setting$row_names
  
  setEPS()
  postscript(paste0(subfig_dir,"subfig6c_sim_score_v2.eps"),
             onefile = F,width = 1,height = 10)
  pheatmap(pd.scores,cluster_rows = F,cluster_cols = F,
           breaks =  c(0,.13,0.5),
           gaps_row = fig5.setting$row_seps-1,
           cellwidth = 12,
           show_rownames = F,
           show_colnames = F,legend = F,
           color = c("yellow","black"))  
  dev.off()
  
  pheatmap(pd.scores,cluster_rows = F,cluster_cols = F,
           breaks =  c(0,.13,0.5),
           gaps_row = fig5.setting$row_seps-1,
           cellwidth = 12,
           show_rownames = T,fontsize_row = 6,
           show_colnames = F,legend = F,
           color = c("yellow","black")) 
  
  setEPS()
  postscript(paste0(subfig_dir,"subfig6c_sim_score_v1.eps"),
             onefile = F,width = 1,height = 10)
  pd.v1 <- data.frame(as.numeric(pd.scores.v1[,1]))
  pheatmap(pd.v1,cluster_rows = F,cluster_cols = F,
           gaps_row = fig5.setting$row_seps-1,breaks = c(-0.1,0.5,1.2),
           cellwidth = 12,
           show_rownames = F,
           show_colnames = F,legend = F,
           color = c("black","yellow"))  
  dev.off()
  
  
  pd.fig6$sim_score$v2 <- pd.scores
  pd.fig6$sim_score$v1 <- pd.scores.v1
  
  saveRDS(pd.fig6,"pd.fig6.Rds")
  }





# init --------------------------------------------------------------------
# loading packages. 
source('../auxilary_functions.R')

maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 

v1.simData <- read.csv(file='./models/mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
tmp.simData <- read.csv(file = "./models/mRNA-Fit-avg-sp-v1e/bestFit.tc.csv",stringsAsFactors = F)
colnames(tmp.simData) <- sub('mRNA','normCnt.frac',colnames(tmp.simData)) 
v1.simData <- (rbind(v1.simData,tmp.simData))
# load input 
load(file = './models/mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata')

# load parmeters  
model.par <- read.csv(file='./models//mRNA-Fit-avg-sp-v1d/result.csv',stringsAsFactors = F,row.names = 1)
model.par <- read.csv(file='./table_model.csv',stringsAsFactors = F,row.names = 1)
gene.dic <- read.csv(file = "../data/mRNA.cluster_old.csv",stringsAsFactors = F,row.names = 2)
model.par$gene.2 <- gene.dic[model.par$ensembleID,"gene"]


# function 
runV1perturb<- function(eg.genes='Ccl7',checkNorm = F, s= "k_deg",pnew=log(2)/15){
  
  p <- plotEgFit(eg.genes = eg.genes)
  
  a<- p$data %>% group_by(gene,type,geno,sti) %>% dplyr::summarise(max= max(normCnt.frac)) %>%
    spread(key=sti,value=max) %>% mutate(LvT.sp = log2(lps/tnf))  %>% 
    ungroup %>% dplyr::select(gene,geno,LvT.sp,type)%>% 
    mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko'))) 
  a$type <- 'Sim.norm'
  
  b<- maxScale.mRNA %>% filter(gene %in% eg.genes) %>% group_by(gene,geno,sti) %>%
    dplyr::summarise(max=max(normCnt.frac)) %>% spread(key=sti,value=max) %>% mutate(LvT.sp=log2(lps/tnf)) %>%ungroup %>% dplyr::select(gene,geno,LvT.sp) %>%
    mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko')),type='Exp.') 
  
  pd.sp <- rbind(a,b)
  
  
  # load single gene parameters 
  pars <- model.par[eg.genes,5:10]
  
  # re-run the model 
  if(checkNorm){
    c <- fun.runSim.v1(pars = pars,g=eg.genes,input.nfkb = input.emsa.ifnar,
                       input.nfkb.mt = input.emsa.ifnarikba)
    c%>% group_by(gene,geno,sti) %>% dplyr::summarise(max= max(mRNA)) %>%
      spread(key=sti,value=max) %>% mutate(LvT.sp = log2(lps/tnf))  %>% 
      ungroup %>% dplyr::select(gene,geno,LvT.sp)%>% 
      mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko'))) 
    a 
  }
  
  # run all the perturbations 
  pars.p <- pars; pars.p[s] <- pnew
  c <- fun.runSim.v1(pars = pars.p,g=eg.genes,input.nfkb = input.emsa.ifnar,
                     input.nfkb.mt = input.emsa.ifnarikba)
  d<- c%>% group_by(gene,geno,sti) %>% dplyr::summarise(max= max(mRNA)) %>%
    spread(key=sti,value=max) %>% mutate(LvT.sp = log2(lps/tnf))  %>% 
    ungroup %>% dplyr::select(gene,geno,LvT.sp)%>% 
    mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko')),type=s) 
  pd.sp <- rbind(pd.sp,d)
}


# main loop ---------------------------------------------------------------
require(parallel)
pd.sp.perturb <- readRDS('~/Dropbox/Projects/DurationDecoding-code/notebooks/kdeg_15m.Rdata')
old.genes <- unique(as.character(pd.sp.perturb$gene))
genes <- model.par$gene.2
venn(list(old.genes,genes))
genes.new <- genes[!genes%in% old.genes]
res<- mclapply(genes,FUN = function(g) runV1perturb(eg.genes = g,
                                              pnew=log(2)/15),mc.cores = 7)
res.2 <- do.call(rbind,res)
saveRDS(file = 'kdeg_15m_new.Rdata',res.2)



# loop for tc data  -------------------------------------------------------
egs  <- c("Ccl5", "Gsap", "Rab15", "Mmp3", "Sod2", "Il1rl1") 
egs  <- c("Slc6a12", "Gsap", "Cd274", "Exoc3l4") 
egs  <- c("Bcl3", "Ccl7", "Tmem132e", "Slco3a1","Mmp9") 
runV1perturb_tc<- function(eg.genes='Ccl7',checkNorm = F, s= "k_deg",pnew=log(2)/15){
  
  p <- plotEgFit(eg.genes = eg.genes)
  a <- p$data %>% filter(type=="Sim.") %>% 
    mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko')))
   
  a$type <- 'Sim.norm'
  
  # load single gene parameters 
  pars <- model.par[eg.genes,5:10]
  
  # re-run the model 
  if(checkNorm){
    c <- fun.runSim.v1(pars = pars,g=eg.genes,input.nfkb = input.emsa.ifnar,
                       input.nfkb.mt = input.emsa.ifnarikba)
    c%>% group_by(gene,geno,sti) %>% dplyr::summarise(max= max(mRNA)) %>%
      spread(key=sti,value=max) %>% mutate(LvT.sp = log2(lps/tnf))  %>% 
      ungroup %>% dplyr::select(gene,geno,LvT.sp)%>% 
      mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko'))) 
    a 
  }
  
  # run all the perturbations 
  pars.p <- pars; pars.p[s] <- pnew
  c <- fun.runSim.v1(pars = pars.p,g=eg.genes,input.nfkb = input.emsa.ifnar,
                     input.nfkb.mt = input.emsa.ifnarikba)
  c <- (c%>% mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko')),type=s) )
  colnames(c) <- sub("mRNA","normCnt.frac",colnames(c))
  
  rbind(a,c)
}


res<- mclapply(egs,FUN = function(g) runV1perturb_tc(eg.genes = g,
                                                    pnew=log(2)/15),mc.cores = length(egs))
res.2 <- do.call(rbind,res)
saveRDS(file = 'kdeg_15m_tc_test.Rdata',res.2)

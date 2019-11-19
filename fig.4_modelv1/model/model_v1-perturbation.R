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
v1.simData$gene<- gdic[v1.simData$gene,"gene.2"]
maxScale.mRNA$gene <- gdic[maxScale.mRNA$gene,"gene.2"]

# load input 
ev <- new.env()
load(file = './models/mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata',ev)
input.emsa.ifnar <- ev$input.emsa.ifnar
input.emsa.ifnarikba <- ev$input.emsa.ifnarikba

# load parmeters  
model.par <- read_csv(file='./table_model.csv')%>%
  as.data.frame()%>%
  column_to_rownames("gene")


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
  pars <- model.par[eg.genes,2:7]
  
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

genes<- (model.par%>%
                   rownames_to_column("gene")%>%
  filter(nrmsd<0.13))$gene

i<-1
res<- lapply(genes,
             FUN = function(g) {
               print(i)
               runV1perturb(eg.genes = g,pnew=log(2)/15) 
               i <<- i+1
             })

res<- mclapply(genes,
               FUN = function(g) runV1perturb(eg.genes = g,
  pnew=log(2)/15),mc.cores = 6)
res <- do.call(rbind,res)
saveRDS(file = 'kdeg_15m_new.Rdata',res)



# loop for tc data  -------------------------------------------------------

runV1perturb_tc<- function(eg.genes='Ccl7',checkNorm = F, s= "k_deg",pnew=log(2)/15){
  
  p <- plotEgFit(eg.genes = eg.genes)
  a <- p$data %>% filter(type=="Sim.") %>% 
    mutate(geno=plyr::revalue(geno,c('wt'='ctrl','mt'='ko')))
   
  a$type <- 'Sim.norm'
  
  # load single gene parameters 
  pars <- model.par[eg.genes,2:7]
  
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


res<- mclapply(genes,
               FUN = function(g) runV1perturb_tc(eg.genes = g,
               pnew=log(2)/15),mc.cores = 6)
res.2 <- do.call(rbind,res)
saveRDS(file = './data/kdeg_15m_tc.Rdata',res.2)

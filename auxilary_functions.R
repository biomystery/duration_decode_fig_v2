
# librarys ----------------------------------------------------------------
require(pheatmap)
require(ggplot2)    
require(deSolve);require(gplots)
require(RColorBrewer)
require(reshape);
require(tidyverse); 
require(org.Mm.eg.db)
require(biomaRt)
require(pracma)
require(cowplot)
require(ggpubr)
require(data.table)
require(gridExtra)
require(tidyverse)


col.map.cap <- c(LPS="#f8766d",TNF="#00ba38",IL1="#619cff")
col.map <- c(lps="#f8766d",tnf="#00ba38",il1="#619cff")

col.sp <- rev(c(colorRampPalette(c( "#efee00", "azure"))(4),
              colorRampPalette(c( "azure", "#1080a3"))(4)))

col.cytoRNA <- colorRampPalette(c("mediumblue", "white", "firebrick1"))(20)
col.cytoRNA <- brewer.pal(n = 9,name = 'Set1')[1]
col.caRNA <- colorRampPalette(c("#06A199", "white", "#E0036C"))(20) # pink green
col.caRNA <- colorRampPalette(c("#B4CA03", "white", "#7F378B"))(20) # pink green
col.caRNA <- brewer.pal(n = 9,name = 'Set1')[2]
gdic <- read.csv('../data/mRNA.cluster_old.csv',stringsAsFactors = F)
rownames(gdic) <- gdic$gene

ccc_fun <- function(x,y){
  # calculate lin's concordance correlation coefficient
  #https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
  #x <- rand(100,1); y<- rand(100,1)
  2*cor(x,y) *sqrt(var(x))*sqrt(var(y))/(var(x)+var(y)+(mean(x)-mean(y))^2)
}

getEntrezID <- function(glist){
  symbol <- mapIds(org.Mm.eg.db,
                   keys=glist,
                   column="ENTREZID",
                   keytype="ENSEMBL",
                   multiVals="first")
  
  ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
  symbol[is.na(symbol)] <- getBM(attributes=c('ensembl_gene_id','entrezgene'),filters = 'ensembl_gene_id',
                                 values  =glist[is.na(symbol)], mart = ensembl)$entrezgene
  symbol
}

getSymbol <- function(glist){
  require(org.Mm.eg.db)
  require(biomaRt)
  ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
  symbol <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
                                 values  =glist, mart = ensembl)$mgi_symbol
  try(symbol[is.na(symbol)] <- mapIds(org.Mm.eg.db,
                       keys=glist[is.na(symbol)],
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first"),silent = T)
  
  
  names(symbol) <- glist
  symbol
}



# Preprocessing -----------------------------------------------------------
getSp<- function(metrics = rsums){
  deltaAuc<-sapply(c(1,4,7,10), function(i) (rsums[,i]/rsums[,i+1]) )# LPS/TNF
  deltaAuc.2 <-sapply(c(1,4,7,10), function(i) (rsums[,i]/rsums[,i+2]) )# LPS/IL1
  deltaAuc.3 <-sapply(c(1,4,7,10), function(i) (rsums[,i+1]/rsums[,i+2]) )# TNF/IL1
  
  deltaAuc.all <- data.frame(
    b1.LT.ctrl = deltaAuc[,1],
    b1.LI.ctrl = deltaAuc.2[,1],
    b1.TI.ctrl = deltaAuc.3[,1],
    b1.LT.ko = deltaAuc[,2],
    b1.LI.ko = deltaAuc.2[,2],
    b1.TI.ko = deltaAuc.3[,2],
    b2.LT.ctrl = deltaAuc[,3],
    b2.LI.ctrl = deltaAuc.2[,3],
    b2.TI.ctrl = deltaAuc.3[,3],
    b2.LT.ko = deltaAuc[,4],
    b2.LI.ko = deltaAuc.2[,4],
    b2.TI.ko = deltaAuc.3[,4]
  ) 
}

# overall -----------------------------------------------------------------
# from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# eg: multiplot(p1, p2, p3, p4, cols=2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# fig3 --------------------------------------------------------------------
getGdata <- function(gname='Ccl5',cnt =all.normCnt,
                     fit=all.fit.data,gd=gene.dic){
  g <- gname; names(g) <- rownames(gd)[which(gd==g)]
  g <- revert(g)
  list(normCntlog2 =cnt[g,],
       g = g,
       kdeg.fit = subset(fit,ensembleID == g))
}
revert <- function(g){
  tmp <- names(g)
  names(tmp)<-g
  tmp
}

plotActD.gene.log2.gg.key <- function(pd,sel= c('U','TNF','LPS')){
  pd.data <- data.frame(log2normCnt = as.numeric(pd$normCntlog2),
                        coldata)
  pd.data$Sti <- factor(pd.data$Sti,levels = c('U','TNF','LPS'))
  #print(pd.data$Sti)
  idx <- unlist(sapply(1:nrow(pd$kdeg.fit), 
                       function(i) return(pd$kdeg.fit$startID[i]:pd$kdeg.fit$endID[i])))
  pd.fit <- pd.data[idx,]
  list(pd.data=pd.data%>%filter(Sti%in%sel),
       pd.fit=pd.fit%>%filter(Sti%in%sel))
}

predict.rlm <- function (object, newdata = NULL, scale = NULL, ...){
  ## problems with using predict.lm are the scale and
  ## the QR decomp which has been done on down-weighted values.
  object$qr <- qr(sqrt(object$weights) * object$x)
  predict.lm(object, newdata = newdata, scale = object$s, ...)
}

plotSimEg<- function(hfs = hf,genotype='wt'){
  pd.1 <- get_fig3_pdata(hfs = hfs,is_ggplot = T,simdata = fun.runSim(hf=hfs,n=1,geno=genotype))
  pd.1$hf <-factor(pd.1$hf,levels = log10(hfs))
  
  p.lines<- ggplot(data=pd.1,aes(x=time,y=value,colour=stimuli)) +geom_line() + facet_grid(hf ~.)
  p.lines <- p.lines + xlab('Time (hr)')+ theme_bw() + theme(strip.text.y=element_blank(),
                                                             axis.text.y=element_blank(),
                                                             axis.title.y=element_blank(),
                                                             legend.position="none") 
  p.lines <- p.lines  +  scale_colour_manual(values = col.map.cap)
}
plotHMsim <- function(geno='mt',forFig=F) {
  # @geno = 'ctrl' or 'mt'
  pwd <- ifelse(forFig,'./data/','../fig.3_hf_hypothesis/data/')
  if(geno=='mt')
    load(paste0(pwd,'half-life-sim-mt.Rdata'))
  else
    load(paste0(pwd,'half-life-sim.Rdata'))
  
  ## ggplot method 
  pd.gg.2 <- get_fig3_pdata(is_ggplot = T,simdata = sim.data)
  p<-ggplot(pd.gg.2,aes(x=time,y=hf,fill=value,group=stimuli)) + geom_tile() + 
    scale_fill_gradientn(colours = brewer.pal(9,name = 'Blues'))
  p<- p + facet_wrap(~stimuli)
  #p + stat_contour(aes(z=value, colour =..level..), binwidth=.1,alpha=.6)
  p<-p + stat_contour(aes(z=value),colour ="azure",linetype=2) + xlab('Time (hr)')  + theme_bw() #half-life(log10 mins)
  p <- p + theme(axis.text.y=element_blank(),
                 axis.title.y=element_blank(),
                 strip.background = element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")
  p<- p + geom_hline(yintercept = c(log10(30),2),linetype=2,color='grey60')
  
  pwd <- "../"
  tf.pd <- readRDS(file = paste0(pwd,'fig.2_knockout/data/subfig2A.rds'))
  if(geno=='mt')   tf.pd.new <-rbind(transTFdata(geno = 'mutant',
                                                 tf.df = tf.pd$data, 
                                                 tps = unique(p$data$time)),
                                     transTFdata(sti="TNF",geno = 'mutant',
                                                 tf.df = tf.pd$data, 
                                                 tps = unique(p$data$time)))
  else 
    tf.pd.new <-     rbind(transTFdata(tf.df = tf.pd$data, 
                                       tps = unique(p$data$time)),
                           transTFdata(sti="TNF",
                                       tf.df = tf.pd$data, 
                                       tps = unique(p$data$time)))
  
  
  tf.p <-ggplot(tf.pd.new,aes(x=Time,y=genotype,fill=value)) + 
    geom_tile()+ scale_fill_gradientn(colours = brewer.pal(9,name = 'Blues')) + facet_grid(~Stimuli)
  tf.p <- tf.p  + theme_bw() 
  tf.p <- tf.p + theme(axis.text=element_blank(),
                       axis.title=element_blank(),
                       strip.text.x = element_blank(),
                       legend.position="none") 
  #
  return(list(p=p,tf.p=tf.p))
}
plotExp <- function(dtype="caRNA",savetofile=T,scale='none',show_fig=F){
  # doScale= 'none','all','geno','wt'
  ev <- new.env()
  if(dtype =="mRNA"){
    load(file = '~/Dropbox/Projects/DurationDecoding-code/Fig_code/fig.4_modelv1/models/mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata',ev)   
    normCnt <- ev$genes.normCount
    tps <- ev$genes.timepoints
  }else if(dtype == 'caRNA'){
    load(file='~/Dropbox/Projects/DurationDecoding-code/Fig_code/fig.5_caRNA/caRNAFit-avg-sp-v6b/caRNAFit-avg-sp-pre-v6.Rdata',ev)
    normCnt <- ev$caRNA.frac.avg
    tps <- ev$caRNA.time
  }
  n.gene <- nrow(normCnt); n.tp <- ncol(normCnt)
  normCnt <- switch (scale,
                     all = normCnt/apply(normCnt,1,max),
                     wt =  normCnt/apply(normCnt[,1:(ncol(normCnt)/2)],1,max),
                     geno = cbind(normCnt[,1:(ncol(normCnt)/2)]/apply(normCnt[,1:(ncol(normCnt)/2)],1,max),
                                  normCnt[,(ncol(normCnt)/2+1):ncol(normCnt)]/apply(normCnt[,(ncol(normCnt)/2+1):ncol(normCnt)],1,max)),
                     none = normCnt
  ) 
  maxScale.df <- data.frame(normCnt.frac = as.vector(t(normCnt)),
                            time = rep(tps,n.gene)/60,
                            gene = rep(rownames(normCnt),each=n.tp),
                            sti = rep(rep(c("lps","tnf",'lps',"tnf"),each=n.tp/4),n.gene),
                            geno = rep(rep(c("wt","wt",'mt',"mt"),each=n.tp/4),n.gene),
                            cond = rep(rep(c("wt.lps","wt.tnf",'mt.lps',"mt.tnf"),each=n.tp/4),n.gene))
  maxScale.df$gene <- as.character(maxScale.df$gene)
  maxScale.df$geno <- relevel(maxScale.df$geno,'wt')
  require(ggplot2)
  if(savetofile)    pdf(file=paste0(dtype,'.data.scale.',scale,'.pdf'),width = 12)
  tmp <- c(seq(1,n.gene,by=24),n.gene)
  for(i in 1:(length(tmp)-1)){
    pd <- maxScale.df[((tmp[i]-1)*n.tp+1):((tmp[i+1]-1)*n.tp),]
    p<- ggplot(pd,aes(x=time,y=normCnt.frac,colour=sti)) + geom_point(aes(shape=geno))+ geom_line(aes(linetype=geno)) + facet_wrap(~ gene, ncol = 8)
    if (show_fig) print(p)
  }
  if(savetofile)  {dev.off();system(paste0("open ",dtype,'.data.scale.',scale,'.pdf'))}
  maxScale.df
}

runModel.v1<- function (nfkb_input,pars,times=manytp){
  k0 <- pars$kb;kt<-pars$kt;k_deg<-pars$k_deg;
  kd <- pars$Kd;tau <- pars$tau;N <- pars$n
  
  mRNAini <- c(mRNA= (kt*nfkb_input$val[1]^N/(kd^N+nfkb_input$val[1]^N)+k0)/k_deg) # init,ss assumption
  
  ## subroutine: The model
  odeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      nfkb <-ifelse(Time<= tau,nfkb_input$val[1],pchipfun(nfkb_input$time,nfkb_input$val)(Time-tau))
      dmRNA    <- kt*nfkb^N/(kd^N+nfkb^N)+k0-k_deg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, odeModel, pars)
  out<-as.data.frame(out)
  out$mRNA[nrow(out)] <- out$mRNA[nrow(out)-1]
  out
}


runModel<- function (nfkb_input,pars,times=manytp){
  ## for simulation
  kb <- pars$k0;kt<-pars$kt;k_deg<-pars$k_deg;kd <- pars$kd;n <- pars$n
  
  mRNAini <- c(mRNA= (kb+kt*(nfkb_input$val[1])^n/(kd^n+(nfkb_input$val[1])^n))/k_deg) # init,ss assumption
  
  ## subroutine: The model
  odeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      #nfkb <-approxfun(nfkb_input$time,nfkb_input$val)(Time)
      nfkb <-pchipfun(nfkb_input$time,nfkb_input$val)(Time)
      dmRNA    <- kb + kt*nfkb^n/(nfkb^n+kd^n) -k_deg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, odeModel, pars)
  out<-as.data.frame(out)
}



fun.runSim.v1 <- function(pars,input.nfkb,input.nfkb.mt,g=gene){
  times <- seq(0,max(input.nfkb$Time),by=0.1);
  
  out  <- sapply(c("LPS","TNF"), function(x)
    runModel.v1(nfkb_input = data.frame(time= input.nfkb$Time,val = input.nfkb[,x]),
                pars= pars,
                times = times),simplify = F)
  
  out.mt  <- sapply(c("LPS","TNF"), function(x)
    runModel.v1(nfkb_input = data.frame(time= input.nfkb.mt$Time,val = input.nfkb.mt[,x]),
                pars= pars,
                times = times),simplify = F)
  
  
  # save output 
  plot.tp_step <- 5 # 5mins step
  
  all.tps <- out$LPS$time
  
  plot.tp_idx <- all.tps %in% seq(all.tps[1],all.tps[length(all.tps)],by=plot.tp_step)
  n.tps <- sum(plot.tp_idx)
  data.ts <- rbind(out$LPS[plot.tp_idx,],
                   out$TNF[plot.tp_idx,],
                   out.mt$LPS[plot.tp_idx,],
                   out.mt$TNF[plot.tp_idx,])
  data.ts$sti <- rep(c("lps","tnf","lps","tnf"),each=n.tps)
  data.ts$geno <- rep(c("wt","wt",'mt',"mt"),each=n.tps)
  data.ts$gene <- rep(g,n.tps*4)
  rownames(data.ts) <- NULL
  data.ts
}  

runModel.v2<-function (caRNA_input,pars,times=manytp){
  
  kt<- pars$kt; kdeg<-pars$k_deg;tau <- pars$tau;
  
  mRNAini <- c(mRNA= caRNA_input$val[1]*kt/kdeg) # init,ss assumption
  #  mRNAini <- c(mRNA=0)
  
  ## subroutine: The model
  ddeModel <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)),{
      #caRNA <-ifelse(Time<= tau,caRNA_input$val[1],approxfun(caRNA_input$time,caRNA_input$val)(Time-tau))
      caRNA <-ifelse(Time<= tau,caRNA_input$val[1],pchipfun(caRNA_input$time,caRNA_input$val)(Time-tau))
      
      dmRNA    <- kt*caRNA - kdeg*mRNA
      return(list(c(dmRNA)))
    })
  }
  
  out   <- ode(mRNAini, times, ddeModel, pars)
  out<-as.data.frame(out)
}

fun.runSim.v2<- function(pars,input.caRNA =pd.ca.exp,caRNA.time=ev$caRNA.time,
                         g="Ccl1"){

  caRNA_tmp <- matrix(input.caRNA[g,],ncol = 4)
  times <- seq(0,max(caRNA.time),by=0.1)
  out<-  sapply(1:4, function(x)
    runModel.v2(caRNA_input = data.frame(time=caRNA.time, 
                                         val =unlist(caRNA_tmp[,x])),
                pars=pars,times),simplify = F)
  
  # save output 
  plot.tp_step <- 5 # 5mins step
  
  all.tps <- out[[1]]$time
  
  plot.tp_idx <- all.tps %in% seq(all.tps[1],all.tps[length(all.tps)],by=plot.tp_step)
  n.tps <- sum(plot.tp_idx)
  data.ts <- rbind(out[[1]][plot.tp_idx,],
                   out[[2]][plot.tp_idx,],
                   out[[3]][plot.tp_idx,],
                   out[[4]][plot.tp_idx,])
  data.ts$sti <- rep(c("lps","tnf","lps","tnf"),each=n.tps)
  data.ts$geno <- rep(c("wt","wt",'mt',"mt"),each=n.tps)
  data.ts$gene <- rep(g,n.tps*4)
  rownames(data.ts) <- NULL
  data.ts
}  

fun.runSim<- function(n=1,hf = c(60,60*2.5),kd=.5,geno='wt'){
  require(parallel)
  
  dataFolder <- '../data/'
  ifnar_input <- read.csv(file=paste0(dataFolder,'Final_EMSA_ifnar.csv'))
  ifnarikba_input <- read.csv(file=paste0(dataFolder,'Final_EMSA_ifnarikba.csv'))
  
  emsa_max <- max(ifnar_input[,2:4])
  ifnar_input[,2:4] <- ifnar_input[,2:4]/emsa_max
  ifnarikba_input[,2:4] <- ifnarikba_input[,2:4]/max(ifnarikba_input[,2:4])
  
  manytp <- seq(0,480,by=0.1)
  times <- ifnar_input$Time
  idx <- which(manytp%in% times)
  idx[7]<- idx[7]-1
  N <- length(hf)
  kb = 0.001;kt = 0.5;#kd=.5;  #n=1;
  pars <- data.frame(k0 = rep(kb,N), 
                     kt = rep(kt,N),
                     kd = rep(kd,N),
                     n = rep(n,N),
                     k_deg = log(2)/hf)
  
  # run LPS 
  if(geno=='wt'){
    nfkb <- data.frame(time=ifnar_input$Time,
                       val =ifnar_input$LPS)
  }  else{
    nfkb <- data.frame(time=ifnarikba_input$time,
                       val =ifnarikba_input$LPS)
  }
  
  
  runParallel <- ifelse(length(hf)>1,T,F)
  n.core <- ifelse(length(hf)>14,14,length(hf))
  
  if(runParallel){
    sim_result <- mclapply(1:length(hf),function(x) runModel(nfkb_input = nfkb,
                                                             pars = pars[x,],times = manytp)$mRNA,
                           mc.cores = n.core)
    
  }else{
    sim_result <- sapply(1:length(hf),function(x) runModel(nfkb_input = nfkb,
                                                           pars = pars[x,],times = manytp)$mRNA)
  }
  
  
  lps_sim_result <- matrix(unlist(sim_result),nrow = N,byrow = T)
  
  #lps_sim_result<- lps_sim_result[,idx]
  # run tnf 
  if(geno=='wt')
    nfkb <- data.frame(time=ifnar_input$Time,
                       val =ifnar_input$TNF)
  else
    nfkb <- data.frame(time=ifnarikba_input$time,
                       val =ifnarikba_input$TNF)
  
  
  if(runParallel){
    sim_result <- mclapply(1:length(hf),function(x) runModel(nfkb_input = nfkb,
                                                             pars = pars[x,],times = manytp)$mRNA,
                           mc.cores = n.core)
    
  }else{
    sim_result <- sapply(1:length(hf),function(x) runModel(nfkb_input = nfkb,
                                                           pars = pars[x,],times = manytp)$mRNA)
  }
  
  tnf_sim_result <- matrix(unlist(sim_result),nrow = N,byrow = T)
  filter.tp <- manytp%in% seq(0,max(manytp),by=5)
  sim.result <- t(rbind(lps_sim_result[,filter.tp],tnf_sim_result[,filter.tp]))
}

get_fig3_pdata <- function(hfs=10^seq(0,3,length.out = 50),
                           simdata=sim.data,
                           tps=seq(0,480,by=5),is_ggplot=F){
  sim.data <- t(simdata)  
  manytp <- seq(0,480,by=0.1)
  pd <- cbind(sim.data[1:length(hfs),],sim.data[(length(hfs)+1):nrow(sim.data),])
  pd <- apply(pd, 1, function(x) x/max(x,na.rm = T))
  
  
  if(is_ggplot){
    pd.gg <- data.frame(t(pd))
    pd.gg <- pd.gg[,-c(ncol(pd.gg)/2,ncol(pd.gg))]
    pd.gg$hf <- log10(hfs)
    pd.gg.2 <- melt(pd.gg,id.vars = 'hf')
    pd.gg.2 <- pd.gg.2[order(pd.gg.2$hf),]
    pd.gg.2$time <- rep(tps[-length(tps)],2*length(hfs))/60
    pd.gg.2$stimuli <- rep(rep(c('LPS','TNF'),each=length(tps)-1),length(hfs))
    return(pd.gg.2)
  }
  t(pd)
}

plotFeaturevsHf<- function(b="b1",a=.5){
  m <- ggplot(pd.peak%>%filter(batch==b),aes(hf,rpkm,colour=stimuli))
  m<-m + geom_point(aes(shape=genotype),alpha=a) + facet_grid(feature~genotype)  + 
    scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')
  m<- m + stat_smooth(se = F,aes(fill=stimuli))
  m
}
plotFeaturevsHf.2<- function(b="b1",a=.5){
  # caRNA 
  m <- ggplot(pd.peak,aes(hf,rpkm,colour=stimuli))
  m<-m + geom_point(alpha=a) + facet_grid(feature~genotype)  + 
    scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')
  m<- m + stat_smooth(se = F,aes(fill=stimuli))
  m
}


transTFdata <- function(sti ="LPS",geno ="ctrl.",tf.df = tf.pd$data,
                        tps = unique(p$data$time)){
  tf.pd.sub <- subset(tf.df,genotype == geno & Stimuli==sti)
  data.frame(
    Time = tps,
    #value=with(tf.pd.sub,approxfun(Time,value)(tps)),
    value=with(tf.pd.sub,pchipfun(Time,value)(tps)),
    genotype= rep(geno,length(tps)),
    Stimuli = rep(sti,length(tps))
  )
}

# fig3S -------------------------------------------------------------------

plotAllEgFit <- function(i){
  genes.iter <- genes[tmp.idx[i]:(tmp.idx[i+1]-1)]
  
  p <- plotEgFit(genes.iter)  #pg <- ggplot_build(p)
  pd_lab <- cbind(p$data %>% group_by(gene) %>% dplyr::summarise(r2=unique(r2)),
                  sti='LPS',type="Sim.")
  pd_lab$gene <- as.character(pd_lab$gene)
  pd_lab$hf <- round(log(2)/fit.r2[pd_lab$gene,'k_deg']);pd_lab$gene <- factor(pd_lab$gene,levels = pd_lab$gene)
  p <- p + geom_text(data= pd_lab ,
                     aes(x = 0,y=1,label=paste0('R2=',signif(r2,2),'\n hf=',hf)),colour='red',
                     hjust=0,vjust=1) 
  
  # sim sp. calculation  & exp sp
  sp.mRNA.iter <- sp.mRNA[genes.iter,c(1,4)]
  sp.mRNA.iter<- sp.mRNA.iter %>% mutate(gene=rownames(sp.mRNA.iter)) %>% 
    gather(key=genotype,value=LvT.sp,-gene)
  sp.mRNA.iter$genotype<- sub('b1.LT.','',sp.mRNA.iter$genotype) 
  sp.mRNA.iter$genotype<- sub('b1.LT.','',sp.mRNA.iter$genotype) 
  sp.mRNA.iter$gene <- factor(sp.mRNA.iter$gene,levels =genes.iter )
  
  sp.mRNA.sim<- p$data %>% group_by(gene,type,geno,sti) %>% dplyr::summarise(mean= mean(normCnt.frac)) %>%
    spread(key=sti,value=mean) %>% mutate(LvT.sp = log2(LPS/TNF))  %>% 
    ungroup %>% select(gene,geno,LvT.sp,type)%>% 
    mutate(geno=revalue(geno,c('wt'='ctrl','mt'='ko'))) 
  names(sp.mRNA.sim) <- sub('geno','genotype',names(sp.mRNA.sim))
  sp.mRNA.iter <- rbind(sp.mRNA.sim, data.frame(sp.mRNA.iter,type='Exp.'))
  p.sp <-     ggplot(sp.mRNA.iter,aes(genotype,LvT.sp))+ 
    geom_bar(stat = "identity") + 
    facet_grid(gene~type,scales = 'free_y')
  
  
  # combine plots 
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 10)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  print(p,vp = vplayout(1,4:10))
  print(p.sp,vp = vplayout(1,1:3))
}

plotEgFit <- function(eg.genes,isFinal=F,use_facet=T,simData=v1.simData,
                      expData=maxScale.mRNA,wide=F){
  
  pd <- rbind(cbind(subset(simData,gene %in% eg.genes),type='Sim.'),
              cbind(subset(expData,gene %in% eg.genes)[,colnames(simData)],type='Exp.'))
  #pd$sti <- toupper(pd$sti);
  pd$geno <- factor(pd$geno,levels = c('wt','mt'))
  pd$gene <- factor(pd$gene,levels=eg.genes);pd$type <- relevel(pd$type,'Exp.')
  p <- ggplot(subset(pd,type== 'Sim.'),aes(x=time,y=normCnt.frac,colour=sti)) + 
    geom_line(aes(linetype=geno))
  p <- p+ geom_point(data = subset(pd,type=='Exp.'),aes(shape=geno)) + 
    geom_line(data = subset(pd,type=='Exp.'),aes(linetype=geno)) + 
    scale_colour_manual(values = col.map)
  if(use_facet){
    if(wide)
      p <- p + facet_grid( type~gene ) 
    else
      p <- p + facet_grid( gene ~type ) 
  } 
    
  
  
  if(isFinal)   return(p + theme_bw() + theme(axis.title.x=element_blank(),
                                       axis.title.y = element_blank(),legend.position='none',
                                       strip.text.x=element_blank(),strip.text.y=element_blank(),
                                       axis.text.x=element_blank(),axis.text.y=element_blank(),
                                       axis.ticks.length = unit(-0.1,'cm')))
  else p + ylab('Expression Level')+xlab('Time (hr)')
}

# fig2s -------------------------------------------------------------------
getAvgData <- function(do.mergebatch=F){
  pd.tc_stat <- sapply(1:8, function(x){
    pd.here  <- NULL 
    for(i in which(k_ord$cluster==x)){
      p <- plotExpSingleGene(i)  
      pd.here <- rbind(pd.here,p$data)
    }
    rownames(pd.here)<- NULL 
    if(do.mergebatch) 
      pd.here %>% group_by(stimuli,time,genotype) %>% summarise(mean_val=mean(value),std_value=sd(value)) %>% mutate(cluster=x)
    else 
      pd.here %>% group_by(stimuli,time,genotype,batch) %>% summarise(mean_val=mean(value),std_value=sd(value)) %>% mutate(cluster=x)
  })
  pd.tc_stat <- do.call(rbind,apply(pd.tc_stat,2,as.data.frame))
  pd.tc_stat
}

  
  getPeakMat <-function(par.windowsize = 100000,debug=F){
    rela.peak.mat <- matrix(0,nrow = nrow(rela.genes.anno),
                            ncol = par.windowsize*2+1);
    rela.gene.region <- data.frame(Chr= rela.genes.anno$Chr,
                                   Start = rela.genes.anno$TSS - par.windowsize,
                                   End = rela.genes.anno$TSS+ par.windowsize)
    
    rownames(rela.gene.region)<-rownames(rela.peak.mat) <- rela.genes 
    for(i in 1:nrow(peak.bed.raw)){
      #i <- 1  
      peak.i <- peak.bed.raw[i,]
      
      # first match chr 
      is.match.chr <- rela.gene.region$Chr == peak.i$Chr #length of 563
      peak.i.d_ps_s <- peak.i$Start - rela.gene.region$Start #length of 563
      peak.i.d_pe_s <- peak.i$End - rela.gene.region$Start #length of 563
      peak.i.d_ps_e <- peak.i$Start - rela.gene.region$End #length of 563
      idx.criteria.1 <- which(peak.i.d_ps_s<=0 & peak.i.d_pe_s>=0 & is.match.chr)
      idx.criteria.2 <- which(peak.i.d_ps_s>=0 & peak.i.d_ps_e<=0 & is.match.chr)
      
      if(length(idx.criteria.1)>0){
        # found indexes/genes for criteria 1
        for(g_id in idx.criteria.1){
          end.id <- min(peak.i$length+peak.i.d_ps_s[g_id],
                        ncol(rela.peak.mat)) 
          overlap.id <- 1:end.id
          rela.peak.mat[g_id,overlap.id] <- rela.peak.mat[g_id,overlap.id]+1
        }
      }else if(length(idx.criteria.2)>0) {
        # found indexes/genes for criteria 2
        for(g_id in idx.criteria.2){
          end.id <- min(peak.i.d_ps_s[g_id]+peak.i$length,ncol(rela.peak.mat))
          overlap.id <- (peak.i.d_ps_s[g_id]+1):end.id
          if(debug) print(data.frame(ps_s=peak.i.d_ps_s[g_id],
                                     ps_e=peak.i.d_ps_e[g_id],
                                     g_id = g_id,
                                     peak.i,
                                     rela.gene.region[g_id,]))
          
          rela.peak.mat[g_id,overlap.id] <- rela.peak.mat[g_id,overlap.id]+1
        }}}
    rela.peak.mat  
  }
 

# fig2 --------------------------------------------------------------------

fun.plotExp<- function(rpkm_all,selc,scale=F,wide=F){
  pd<- as.data.frame(rpkm_all[selc,c(1:14,22:35)])
  pd$gene <- rownames(pd)
  pd <- (melt(pd,id.vars = c("gene")))
  pd <- pd[order(pd$gene),]
  pd$Time <- rep(c(0,.5,1,2,3,5,8),length(selc)*2)
  pd$genotype <- rep(rep(c('ctrl.','mt'),each=14),length(selc))
  pd$stimuli <- rep(rep(c('LPS','TNF'),each=7),length(selc)*2)
  
  
  p<- ggplot(data = pd,aes(x=Time,y=value,colour=stimuli)) + geom_line(aes(linetype=genotype)) 
  if(wide) 
    p <- p+ facet_grid(.~gene,scales = 'free_y')
  else
    p <- p+ facet_wrap(~gene,ncol = 1,scales = 'free_y')
  p <- p + geom_point(aes(shape=genotype))
  p <- p + theme_bw() + scale_x_continuous(breaks = unique(pd$Time)) +  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = .5))
  p<- p + xlab('Time (hr)') + ylab('RPKM')+  scale_colour_manual(values = col.map.cap)
  if(scale){
    pd.new<- p$data %>% group_by(gene,genotype) %>% mutate(value=value/max(value))
    p%+% pd.new + ylab('Expression Level')
  }else p
}

rowMax <- function(x) apply(x, 1, max)
rowMin <- function(x) apply(x, 1, min)

plotKmeansk<- function(df,seed=30,...){
  kscans <-c(2:20,25,30,40,50)#,100,nrow(pd)-1)
    kmeans.withinss <- sapply(kscans, function(x) {
    set.seed(seed)
    kmeans(df,centers = x,...,
           iter.max = 2000,nstart = 25)$tot.withinss})
  k <- 9
  par(mfrow=c(1,2))
  plot(kscans,kmeans.withinss,type='b',col=ifelse(kscans==k,2,1),pch=16,
       xlim=c(0,50),xlab='k clusters')
  abline(v=c(k),col=2)
  plot(kscans[-1],diff(kmeans.withinss)/diff(kscans),type='b',col=ifelse(kscans[-1]==k,2,1),pch=16)
  abline(v=c(k),col=2)
}



loadExpData <- function(){
  dataFolder<-"~/Dropbox/Projects/DurationDecoding/data/RNAseq/processed_data/"
  load(file=paste0(dataFolder,'lrpkmD_clean_data_all.Rdata'))
  
  exp_data <- 2^lrpkm_all[1:28]; 
  exp_data <- exp_data[,c(1:7,15:21,8:14,22:28)]
  exp_data <- exp_data/rowMax(exp_data[,1:14])
  exp_data
}

genSimData <- function(){
  ### load simulation data 
  sim_data<- cbind(
    lps_sim_result,
    tnf_sim_result,
    lps_mt_sim_result,
    tnf_mt_sim_result
  )
  #dim(sim_data)
  colnames(sim_data)
  sim_data <- sim_data/rowMax(sim_data[,1:14])
  sim_data
}


plotLegend <- function(cols,bks,fnames){
  require(RColorBrewer)
  #bks <- seq(round(min(rsums)),round(max(rsums))+1,length.out = 6)
  #cols<- colorRampPalette(c( "white", "blueviolet"))(5)
  #fnames<-'tmp.eps'
  setEPS()
  postscript(fnames,onefile = F,width = 0.1,height = .1*2*length(bks))
  par(mar=c(0.02, 0.04, 0.04,0.1))
  barplot(rep(1,length(bks)-1),width=diff(bks),space = 0,border = NA,
          col = cols,axes = F,horiz = T,
          xaxs = "i", yaxs = "i",xlim = c(0,1),
          ylim=c(0,sum(diff(bks)))
          #xlim=c(0,1)
          )
  axis(4,labels = F,tck = 0.2,at=cumsum(diff(bks)))
  box()
  dev.off()
}

pheatmap.2 <- function(pd,bks,cols,cwid = 8,rn=F,...){
  pheatmap(pd, breaks = bks,
           cluster_cols = F,cluster_rows = F,cellwidth = cwid,
           legend_breaks = bks,
           gaps_row = seps-1,
           color = cols,
           scale='none',
           show_rownames = rn,show_colnames=F,...)
}


# fig1 --------------------------------------------------------------------


getSplist <- function(sti1='L',sti2='T',geno='ctrl',sp.th=2^.5,sp.mat = sp.mat.1){
  eval(parse(text = paste0("df=sp.mat$b1.",sti1,sti2,".",geno)))
  sp.mat.glist <- list(a = df >=sp.th,
                       b = df <= -sp.th,
                       c = (df<sp.th & df > -sp.th))
  # sp.mat.glist <- list(a = df >=sp.th,
  #                      b = df <=1/sp.th,
  #                      c = (df<sp.th & df >1/sp.th))
  # 
  names(sp.mat.glist) <- c(paste0(sti1,sti2,'.',sti1),
                           paste0(sti1,sti2,'.',sti2),
                           paste0(sti1,sti2,'.nosp'))
  sp.mat.glist
}

getVennData <- function(sp.th.val=2^.5,gtype='ctrl',sp.mat.1=sp.mat,simplify=F){
   if(simplify)
     pd <- getSplist(sp.th=sp.th.val,geno=gtype,sp.mat = sp.mat.1)
   else
     pd <-c(getSplist(sp.th=sp.th.val,geno=gtype,sp.mat = sp.mat.1),
            getSplist(sti1 = "L",sti2="I",sp.th=sp.th.val,geno=gtype),
            getSplist(sti1 = "T",sti2="I",sp.th=sp.th.val,geno=gtype))
               
  pd<-sapply(names(pd), function(x) pd[[x]])
  require(limma)
  tmp <- vennCounts(pd)
  pd.venn <- tmp[tmp[,'Counts']>=1,]
  pd.venn.counts <- pd.venn[,"Counts"]
  pd.venn <- pd.venn[,-ncol(pd.venn)]
  
  pd <- data.frame(category=sapply(1:nrow(pd.venn), function(x) paste(names(which(pd.venn[x,]==1)),collapse = '-')),
                   Counts= pd.venn.counts)
  pd  
}

donutPlot <- function(dat = pd.venn.ctrl,inside.v="LT.nosp-LI.nosp-TI.nosp",trans=T,group.v=NULL){
  if(trans){
    dat$category <- relevel(dat$category,ref=inside.v)
    x.mid <- dat$Counts[dat$category==inside.v]/sum(dat$Counts)
    tmp.row <- subset(dat,dat$category==inside.v)
    dat <- subset(dat,dat$category!=inside.v)
    
    dat$fraction = dat$Counts/ sum(dat$Counts)
    dat  <- dat[order(dat$fraction),]
    dat$ymax <- cumsum(dat$fraction)
    dat$ymin <- c(0,head(dat$ymax,n=-1))
    dat$xmax <- rep(1,nrow(dat))
    dat$xmin <- rep(x.mid,nrow(dat))
    
    tmp.row$fraction = 1; tmp.row$ymax=1;tmp.row$ymin=0;tmp.row$xmax=x.mid;tmp.row$xmin=0;
    dat <- rbind(dat,tmp.row)
  }
  require(RColorBrewer)
  p<-ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin)) +
    geom_rect() +
    coord_polar(theta="y") + 
    xlim(c(0, 1)) + 
    theme_bw() +
    theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank()) +
    #scale_fill_manual(values =c('grey',colorRampPalette(brewer.pal(9,"Set1"))(nlevels(dat$category)-1)))
    scale_fill_manual(values =c('grey',colorRampPalette(col.map[1:2])(nlevels(dat$category)-1)))
  if(!trans) eval(parse(text=paste0('p<-last_plot()+ facet_wrap(~ ',group.v,')')))
  p
} 

plot.spPie <- function(type='peak',gtype.1='ctrl',fig_folder) # or avg
{
  sp.mat.2 <- read.csv(file=paste0('./data/mRNA.sp.',type,'.csv'),row.names = 1,stringsAsFactors = F)
  pd.venn.ctrl <- getVennData(simplify = T,sp.th.val = .5,sp.mat = sp.mat.2,gtype = gtype.1)
  pd.venn.ctrl.th <- getVennData(sp.th.val = 1,simplify = T,sp.mat=sp.mat.2,gtype = gtype.1)
  
  
  # threshold = 1.4
  setEPS()
  postscript(paste0(fig_folder,"subfig1c_pie_th1_",type,"_gtype_",gtype.1,".eps"),
             onefile = F,width = 3,height = 3)
  pie(x = pd.venn.ctrl$Counts,labels = "",col = c(gray(.8),col.map['TNF'],col.map['LPS'])) #mRNA.sp.peak.csv
  dev.off()
  # threshold = 2
  setEPS()
  postscript(paste0(fig_folder,"subfig1c_pie_th2_",type,"_gtype_",gtype.1,".eps"),
             onefile = F,width = 3,height = 3)
  pie(x = pd.venn.ctrl.th$Counts,labels = "",col = c(gray(.8),col.map['TNF'],col.map['LPS'])) #mRNA.sp.peak.csv
  par(mfrow=c(1,1))
  dev.off()
  
  print(pd.venn.ctrl)
  print(pd.venn.ctrl.th)
  
}
# fig1s -------------------------------------------------------------------

plotExpSingleGene <- function(gene_idx=1,maxFrac=T){
  require(ggplot2)
  g<- rownames(pd)[gene_idx]
  pd.g <- data.frame(value = as.numeric(pd[gene_idx,]),
                     gene = g,
                     time = rep(c(0,.5,1,2,3,5,8),12),
                     stimuli  = rep(rep(c('lps','tnf','il1'),each=7),4),
                     genotype = rep(rep(c('ctrl.','mt'),each=21),2),
                     batch = rep(c('b1','b2'),each=42))
  pd.g$stimuli <- factor(pd.g$stimuli,levels=c('lps','tnf','il1'))
  ggplot(pd.g,aes(time,value,colour=stimuli)) + 
    geom_line(aes(linetype=genotype)) + 
    geom_point(aes(shape=genotype)) + 
    facet_wrap(~batch) + ggtitle(g) 
  if(maxFrac) 
    last_plot() + ylab('maxFrac')
  else
    last_plot() + ylab('rpkm')
  
}

plotExpSingleGene.2 <- function(gene_idx=1,maxFrac=T){
  # for caRNA 
  require(ggplot2)
  g<- rownames(pd)[gene_idx]
  pd.g <- data.frame(value = as.numeric(pd[gene_idx,]),
                     gene = g,
                     time = rep(c(0,.5,1,2,4,8),4),
                     stimuli  = rep(rep(c('lps','tnf'),each=6),2),
                     genotype = rep(rep(c('ctrl.','mt'),each=12),2) )
  pd.g$stimuli <- factor(pd.g$stimuli,levels=c('lps','tnf'))
  ggplot(pd.g,aes(time,value,colour=stimuli)) + 
    geom_line(aes(linetype=genotype)) + 
    geom_point(aes(shape=genotype)) + 
    facet_wrap(~batch) + ggtitle(g) 
  if(maxFrac) 
    last_plot() + ylab('maxFrac')
  else
    last_plot() + ylab('rpkm')
  
}


# fig4 --------------------------------------------------------------------
cmp_adjR2.2 <- function(v1.simData,rep.num=5,weigh_pt=T){ # weigth on both peak and peak location
  
  # norm residue 
  v1.simData.exp <- v1.simData%>%filter(time %in% maxScale.mRNA$time)
  maxScale.mRNA <- maxScale.mRNA%>% filter(gene %in% v1.simData$gene)
  maxScale.mRNA <- maxScale.mRNA[,colnames(v1.simData.exp)]
  maxScale.mRNA$sti <- as.character(maxScale.mRNA$sti)
  maxScale.mRNA$geno <- as.character(maxScale.mRNA$geno)
  norm.res <- maxScale.mRNA
  norm.res$res <- v1.simData.exp$normCnt.frac - maxScale.mRNA$normCnt.frac
  
  # repeated / weighted peak residuals 
  
  peak.sim <-  v1.simData%>% group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%dplyr::select(1:5) %>% 
    dplyr::slice(rep(1:n(),each=rep.num))
  
  peak.exp <-   maxScale.mRNA%>% group_by(paste(sti,geno,gene,sep="-"))%>% 
    filter(normCnt.frac==max(normCnt.frac))%>% #max in sti,gene,geno
    ungroup()%>%dplyr::select(1:5) %>% 
    dplyr::slice(rep(1:n(),each=rep.num))
  
  all.equal(with(peak.sim,paste(sti,geno,gene,sep="-")),
            with(peak.exp,paste(sti,geno,gene,sep="-")))
  
  # move peak time as data.points
  if(weigh_pt){ 
    peak.sim<- (peak.sim %>% rbind(peak.sim%>%mutate(normCnt.frac=time)))
    peak.exp<- (peak.exp %>% rbind(peak.exp%>%mutate(normCnt.frac=time)))}
  peak.exp$res <- peak.sim$normCnt.frac - peak.exp$normoCnt.frac;
  all.res <- rbind(peak.exp,norm.res)
  
  
  pd<- all.res%>% group_by(gene) %>% 
    summarise(rmsd = sqrt(sum(res^2)/length(res)),
              mean_data= mean(normCnt.frac),
              range_data=max(normCnt.frac)-min(normCnt.frac),
              r2.old = 1-var(res)/var(normCnt.frac),
              r2.new = 1- sum(res^2)/sum((normCnt.frac - mean(normCnt.frac))^2))%>%
    mutate(nrmsd_m = rmsd/mean_data,
           nrmsd_r = rmsd/range_data)
  pd %>% gather(r2.type,r2,2:3)
  print(signif(quantile(pd$r2.new-pd$r2.old,c(0,.25,.5,.75,1)),2))
  p<-ggplot(pd,aes(rmsd,r2.old)) + geom_point()+ 
    geom_point(aes(rmsd,r2.new),colour=2,shape=2)+
    ylab("R2") + ggtitle("R2 new (red) vs. R2 old (black)")
  list(p=p,pd=pd)
}
plotSp<- function(genotype='wt'){
  plot((pd.sp %>% filter(type=='Exp.'& geno==genotype))$sp[ord],1:n.gene,
       pch=16,xlab='',ylab='',xlim=range(pd.sp$sp),
       yaxt='n')
  sapply(0:n.gene,function(x) abline(h=x+.5,cex=.75,lty=1,col='grey'))
  points((pd.sp %>% filter(type=='Sim.'& geno==genotype))$sp[ord],1:n.gene,
         pch=1,col=2)
  abline(v=0)
}

plotSp.2<- function(genotype='ctrl',pd=pd.sp.perturb){
  plot((pd %>% filter(type=='Exp.'& geno==genotype))$LvT.sp[ord],1:n.gene,
       pch=16,xlab='',ylab='',xlim=range(pd$LvT.sp),
       yaxt='n')
  sapply(0:n.gene,function(x) abline(h=x+.5,cex=.75,lty=1,col='grey'))
  points((pd %>% filter(type=='Sim.norm'& geno==genotype))$LvT.sp[ord],1:n.gene,
         pch=1,col=2)
  points((pd %>% filter(type=='k_deg'& geno==genotype))$LvT.sp[ord],1:n.gene,
         pch=1,col=3)
  abline(v=0)
}

plotSp.3<- function(genotype='ctrl',pd=pd.sp.perturb,pd.down = pd.sp.perturb.2){
  plot((pd %>% filter(type=='Exp.'& geno==genotype))$LvT.sp[ord],1:n.gene,
       pch=16,xlab='',ylab='',xlim=range(pd$LvT.sp),
       yaxt='n')
  sapply(0:n.gene,function(x) abline(h=x+.5,cex=.75,lty=1,col='grey'))
  points((pd %>% filter(type=='Sim.norm'& geno==genotype))$LvT.sp[ord],1:n.gene,
         pch=1,col=1)
  points((pd %>% filter(type=='k_deg'& geno==genotype))$LvT.sp[ord],1:n.gene,
         pch=1,col=2)
  points((pd.down %>% filter(type=='k_deg'& geno==genotype))$LvT.sp[ord],1:n.gene,
         pch=1,col=3)
  
  abline(v=0)
}

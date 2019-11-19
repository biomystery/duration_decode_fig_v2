## load data 
#setwd("./mRNA-Fit-avg-sp-v1c")
load(file='mRNA-Fit-avg-sp-v1c-pre.Rdata')#,list=c("genes.kdeg","genes.normCount",
#        "input.emsa.ifnarikba",
#        "input.emsa.ifnar", "genes.timepoints"
source(file = "mRNA-Fit-avg-sp-v1c-funs.R") #rowMax, runModel.v1, runOptim.v1





# main 
# round 2: fit control + knockout -----------------------------------
N=2 
x <- g <- "Bcl3"

lb <- c(-3,-3,-3,log10(genes.kdeg.lb[x]),-5);
ub<- c(3,3,3,log10(genes.kdeg.ub[x]),log10(120));
initP <- c(-1,-1,-1,log10(genes.kdeg[x]),-1)
initP <- g.fit$par; initP[5] <- log10(45)
g.fit <-  runOptim.v1(initP = initP,lb = lb,ub = ub,
                      g.name = x,input.nfkb = input.emsa.ifnar,
                      input.nfkb.mt = input.emsa.ifnarikba)

# plot  -------------------------------------------------------------------
require(ggplot2)
pd.exp <- data.frame(g.fit$expdata,
                     sti=rep(c('lps','tnf','lps','tnf'),each=7),
                     genotype=rep(c('ctrl','mt'),each=14))
head(g.fit$bestFit.mt$LPS)
d1 <- data.frame(val=g.fit$bestFit.mt$LPS$mRNA/g.fit$simScale.mt,
           time=g.fit$bestFit.mt$LPS$time,
           genotype='mt',
           sti='lps')
d2 <- data.frame(val=g.fit$bestFit.mt$TNF$mRNA/g.fit$simScale.mt,
                 time=g.fit$bestFit.mt$TNF$time,
                 genotype='mt',
                 sti='tnf')
d3 <- data.frame(val=g.fit$bestFit.ctrl$TNF$mRNA/g.fit$simScale.wt,
                 time=g.fit$bestFit.ctrl$TNF$time,
                 genotype='ctrl',
                 sti='tnf')
d4 <- data.frame(val=g.fit$bestFit.ctrl$LPS$mRNA/g.fit$simScale.wt,
                 time=g.fit$bestFit.ctrl$LPS$time,
                 genotype='ctrl',
                 sti='lps')
pd.sim <- rbind(d1,d2,d3,d4); pd.exp$genotype <- relevel(pd.exp$genotype,'ctrl')
pd.sim$genotype <- relevel(pd.sim$genotype,'ctrl')
ggplot(pd.exp,aes(time,val,colour=sti)) + geom_point(aes(shape=genotype))
last_plot() + geom_line(data = pd.sim,aes(linetype=genotype))

rm(list=ls())
require(RColorBrewer)
require(pheatmap)

load(file = './caFit.Rdata')

calR2 <- function(x){
  tmp <- caFit[[x]]
  1 - var(tmp$residue)/var(tmp$expdata$val)
  
}

caFit.R2 <- sapply(1:length(caFit), calR2)


mRNAfit.R2 <- read.csv(file='../fig.4_modelfit/data/v1-ctrl-goodness-of-fit.csv',
         header = T,row.names = 1)


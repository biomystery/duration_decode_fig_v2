rm(list=ls())
require(RColorBrewer)
require(pheatmap)

source(file = './caRNAFit-avg-sp-v2/caRNAFit-avg-sp-funs-v2.R')
load(file = './caRNAFit-avg-sp-v2/caRNAFit-avg-sp-pre-v2.Rdata')
load(file = './caRNAFit-avg-sp-v2/caRNAFit-avg-sp-run-v2.Rdata')

nm <- 1:length(caFit)
for(i in 1:max(nm)) nm[i] <- caFit[[i]]$gene
names(caFit) <- nm

caFit1.R2 <- sapply(1:length(caFit), calR2)

files <- paste0("./caRNAFit-avg-sp-v",2:5,"/caRNAFit-avg-sp-run-v",2:5,'.Rdata')

for (i in 1:length(files)){
  load(file = files[i])
  eval(expr = parse(text = paste0("caFit",i,".R2 <- sapply(1:length(caFit), calR2)")))
}



# v2 - fit kt; v3 - delay + kt , v4 - delay + kt + kdeg , v5 - kt + kdeg
models <- c('kt','tau+kt','tau+kt+kdeg','kt+kdeg')
cols <- c('#999999','#E69F00',"dodgerblue","blueviolet","brown1","brown4")
plot(density(caFit3.R2),col=cols[3],main='R2 distribution',xlab = "R2",lwd=2)
lines(density(caFit2.R2),col=cols[2])
#lines(density(caFit4.R2),col=cols[3])
lines(density(caFit1.R2),col=cols[1])

sum(caFit3.R2>=0.75)



library(easyGgplot2)
pd <- data.frame( r2 = c(caFit1.R2,caFit2.R2,caFit3.R2),
                  models = rep(c('kt',"tau_kt","tau_kt_kdeg"),each=79))

plot<-ggplot2.histogram(data=pd, xName='r2', groupName='models',
                        addMeanLine=TRUE, showLegend=T,
                        groupColors=c('#999999','#E69F00',"dodgerblue"), alpha=0.5,
                        backgroundColor="white")

plot<-ggplot2.customize(plot, xtitle="R2 ", ytitle="Count",
                        showLegend=T, 
                        mainTitle="R2 by models")  + stat_bin()
plot
?ggplot2.histogram

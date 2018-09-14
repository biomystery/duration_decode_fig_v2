
# plot caRNA seq (max scale)----------------------------------------------------------
rm(list=ls());

cpm.avg.gene <- read.csv(file='../data/caRNA-cpm-max-scale.csv',header = T,
                         row.names = 1,stringsAsFactors = F) 

# plot data 

pd <-cpm.avg.gene;
fun.scale <- function(pd){
  idx.I <- (substr(colnames(pd),1,1)== "I")
  rmax.I <- apply(pd[,idx.I],1,max)
  idx.D <- (substr(colnames(pd),1,1)== "D")
  rmax.D <- apply(pd[,idx.D],1,max)
  return(cbind(pd[,idx.I]/rmax.I,pd[,idx.D]/rmax.D))
}

pd <- fun.scale(cpm.avg.gene) # scale to each genotpye

write.csv(file='../data/caRNA-cpm-max-geno-scale.csv',pd)


# caRNA zscale sp ---------------------------------------------------------

rsums <- sapply(0:3, function(i)
  rowSums(pd[,(i*6+1):((i+1)*6)]))

sp.tl <-sapply(c(1,3), function(i) (rsums[,i+1]/rsums[,i]) )# TNF/LPS 
colnames(sp.tl) <- c("wt","mt")
write.csv(file='caRNA-max-avg-sp.csv',sp.tl)

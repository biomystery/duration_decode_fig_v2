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
  barplot(rep(1,length(bks)),space = 0,border = NA,
          col = cols,axes = F,horiz = T,
          xaxs = "i", yaxs = "i",
          ylim=c(0,length(bks)-1)
          #xlim=c(0,1)
          )
  axis(4,labels = F,tck = 0.2)
  box()
  dev.off()
}

pheatmap.2 <- function(pd,bks,cols,cwid = 8,...){
  pheatmap(pd, breaks = bks,
           cluster_cols = F,cluster_rows = F,cellheight = 1,cellwidth = cwid,
           legend_breaks = bks,
           gaps_row = seps-1,
           color = cols,
           scale='none',
           show_rownames = F,show_colnames=F,...)
}

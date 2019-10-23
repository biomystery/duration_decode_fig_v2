
# merge_scores ------------------------------------------------------------
fits <- c('./1st/','./2nd/','./3rd/')
getScore <- function(f){
  fit <- readRDS(paste0(f,'allRes.Rds'))
  unlist(lapply(fit,function(x) x$fitRes$value))
}

scores <- lapply(fits,getScore)
scores.df <- do.call(cbind,scores)
colnames(scores.df)<-c('V1','V2','V3')

write.csv(file='fit_scores_cmp.csv',scores.df)



# statistics ------------------------------------------------------------------
table(scores.df[,1]< scores.df[,3])
which(scores.df[,1]> scores.df[,3])
require(dplyr)
require(tidyr)
scores.df <- as.data.frame((scores.df))
scores.df$gene=rownames(scores.df)
scores.df.2<- scores.df%>%gather(model,score,1:2)
scores.df.final <- scores.df.2%>%group_by(gene)%>%
  summarise(minScore = min(score),
         minModel = model[which.min(score)])

write.csv(file='./data/final.score.csv',scores.df.final)

# heatmap 
require(ggplot2)
rownames(scores.df.final)<- scores.df.final$gene
ggplot(scores.df.final,aes(minScore))+geom_histogram(color="white")


# parameters --------------------------------------------------------------
getPars <- function(f){
  fit <- read.csv(paste0(f,'pars.csv'),stringsAsFactors = F)
  fit
}

pars.all <- lapply(fits,getPars)
names(pars.all)<- c('V1','V2','V3')
final.score <- read.csv(file='./data/final.score.csv',row.names = 1,stringsAsFactors = F)

final.pars <- apply(final.score, 1, function(s) 
{
  s=as.data.frame(t(s),stringsAsFactors = F)
  s%>%left_join(pars.all[[s$minModel]])%>%select(-one_of("k3"))
})
  
final.pars<-do.call(rbind,final.pars)

fwrite(final.pars,'/Users/frank/Dropbox/Projects/DurationDecoding-code/supplemental_tables/table_model_v2.csv')

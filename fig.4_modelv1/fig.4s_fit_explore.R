source("../auxilary_functions.R")
# selected genes  -------------------------------------------------------
ev <- new.env()
load('./clustering_cateIII_clean.Rdata',ev)
glist.v1 <- with(ev, unique(pd$gene[pd$v1]))
glist.v1.left <- with(ev,unique(pd$gene[!pd$v1]))
res <- read.csv(file='./v1f_v1de_cmp.csv',stringsAsFactors = F)

selc.genes <- (res%>%filter(v2<v1))$gene


# load best fit data emsa
v1.simData <- read.csv(file='./mRNA-Fit-avg-sp-v1d/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v1.simData) <- sub("mRNA","normCnt.frac",names(v1.simData)) 
maxScale.mRNA <- plotExp(dtype = 'mRNA',scale = "geno",savetofile = F)
maxScale.mRNA$cond <- NULL 
selc.genes <- selc.genes[selc.genes%in% v1.simData$gene]
p <- plotEgFit(eg.genes = selc.genes,wide = F)



# load caRNA fit data 
v2.simData <- read.csv(file='../fig.5_caRNA/caRNAFit-avg-sp-v6a/bestFit.tc.csv',
                       stringsAsFactors = F)
names(v2.simData) <- sub("mRNA","normCnt.frac",names(v2.simData)) 
v2.simData$sti <- toupper(v2.simData$sti)
head(p$data);head(v2.simData)
p.2 <- p %+% rbind(p$data,data.frame(v2.simData%>%filter(gene %in%selc.genes),type='caSim.'))


# load the caRNA data
ev <- new.env()
load("../fig.5_caRNA/caRNAFit-avg-sp-v6a/caRNAFit-avg-sp-pre-v6.Rdata",ev)
pd.ca.exp <- ev$caRNA.frac.avg[selc.genes,]
pd.ca.exp <- pd.ca.exp%>% mutate(gene=rownames(pd.ca.exp))%>%
  gather(key = sample,value = normCnt.frac,1:24)%>%
  mutate(time=sapply(sample,function(x) as.numeric(unlist(strsplit(x,split = "_"))[3])),
         sti = sapply(sample,function(x) unlist(strsplit(x,split = "_"))[2]),
         geno= sapply(sample,function(x) unlist(strsplit(x,split = "_"))[1]),
         type="caExp.")
geno.dic <- c(I="wt",D="mt");sti.dic <- c(L='LPS',T='TNF')
pd.ca.exp$sti <-sti.dic[pd.ca.exp$sti]
pd.ca.exp$geno <-geno.dic[pd.ca.exp$geno]
pd.ca.exp$gene <- factor(pd.ca.exp$gene,levels = selc.genes)

p.2 + geom_point(data = pd.ca.exp,
                 aes(time,normCnt.frac,shape=geno,colour=sti)) +
  geom_line(data = pd.ca.exp,
            aes(linetype=geno))+ theme_bw() + theme(axis.title = element_blank(),axis.text = element_blank(),
                                                    legend.position = "none")


#ggsave(p.2,filename = paste0(fig_folder,'v2_improve.png'),height = 4,width = 9)

# r2 plot 
pd.r2$gene <- droplevels(pd.r2$gene);
pd.r2$gene <- factor(pd.r2$gene,levels = levels(p.2$data$gene))
p.r2 <- ggplot(pd.r2,aes(type,r2)) + geom_bar(stat = 'identity') + facet_grid(.~gene)
p.r2 + theme_bw() + theme(axis.title = element_blank(),
                          legend.position = "none")+ scale_y_continuous(position = 'right')

# overlay with half-life 
ev <- new.env()
load("./mRNA-Fit-avg-sp-v1d/mRNA-Fit-avg-sp-v1c-pre.Rdata",ev)
pd.hfs <- round(log(2)/ev$genes.kdeg)
pd.hfs[selc.genes]

# overlay with changing half-life 
ev <- new.env()
load("../fig.5_caRNA/mRNA-Fit-avg-sp-v2a/mRNA-Fit-avg-sp-v2-pre.Rdata",ev)
pd.hfs <- cbind(ev$genes.kdeg.lps,ev$genes.kdeg.tnf)
pd.hfs <- data.frame(pd.hfs[selc.genes,],gene=selc.genes)
pd.hfs <- pd.hfs %>% gather(sample,hfs,1:10) %>% mutate(hfs=round(log(2)/hfs))

pd.hfs$time  <- sapply(pd.hfs$sample,function(x) 
  as.numeric(unlist(strsplit(x,split="_"))[2]))
pd.hfs$sti  <- sapply(pd.hfs$sample,function(x) unlist(strsplit(x,split="_"))[1])
pd.hfs$time[is.na(pd.hfs$time)] <- 0 
pd.hfs$sti[pd.hfs$sti=="U"] <- "LPS"; pd.hfs$sti[pd.hfs$sti=="U.1"] <- "TNF"
pd.hfs$hfs[is.na(pd.hfs$hfs)] <- min(pd.hfs$hfs,na.rm = T)
pd.hfs$gene <- factor(pd.hfs$gene,levels = selc.genes)
p.hf <- ggplot(pd.hfs,aes(time,hfs,colour=sti,group=sti)) + 
  geom_line()

p.hf + theme_bw() + facet_grid(.~gene)+
  theme(axis.title = element_blank(),
        legend.position = "none",axis.text.x = element_blank())+ 
  scale_y_continuous(position = 'right') + 
  scale_colour_manual(values = col.map)

## load v2a sim
v1a.simData <- read.csv(file='../fig.5_caRNA/mRNA-Fit-avg-sp-v2a/bestFit.tc.csv',
                        stringsAsFactors = F)
names(v1a.simData) <- sub("mRNA","normCnt.frac",names(v1a.simData)) 
v1a.simData$sti <- toupper(v1a.simData$sti)

p.3 <- p.2 %+% rbind(p.2$data,data.frame(v1a.simData%>%filter(gene %in%selc.genes),type='V1aSim.'))

v2a.simData <- read.csv(file='../fig.4_modelv1/caRNA-Fit-avg-sp-v7/bestFit.tc.csv',
                        stringsAsFactors = F)
names(v2a.simData) <- sub("mRNA","normCnt.frac",names(v2a.simData)) 
v2a.simData$sti <- toupper(v2a.simData$sti)
p.4 <- p.3 %+% rbind(p.3$data,data.frame(v2a.simData%>%filter(gene %in%selc.genes),type='V2aSim.'))


## load 
# nfkbie  -----------------------------------------------------------------
tnf3hr.hf.data <- data.frame( time=c(0,.5,1,3,6),
                              val=c(11.27,11.1,10.92,9.62,8.15))                          
tnf3hr.hf.data <- data.frame( time=c(0,.5,1,3,6),
                              val=c(6.37,7.93,7.78,6.84,4.72))                          
tnf3hr.hf.data <- data.frame( time=c(0,.5,1,3,6),
                              val=c(10.73,11.06,10.78,8.67,7.25))                          

tnf3hr.hf.fit <- lm(val~time,tnf3hr.hf.data)
par(mfrow=c(2,2))
plot(tnf3hr.hf.fit.2)
par(mfrow=c(1,1))
tnf3hr.hf.fit.2 <- lm(val~time,tnf3hr.hf.data[1:3,])
plot(tnf3hr.hf.data)
abline(tnf3hr.hf.fit)
abline(tnf3hr.hf.fit.2,col=2)
summary(tnf3hr.hf.fit); summary(tnf3hr.hf.fit.2)
car::Anova(tnf3hr.hf.fit,tnf3hr.hf.fit.2)

tnf3hr.hf.data$val <- round(2^tnf3hr.hf.data$val -1) 
glm.fit<- glm(val ~ time, data=tnf3hr.hf.data,  family="poisson"(link="log"))
tnf3hr.hf.data$val <- log(tnf3hr.hf.data$val)
glm.fit.2<- lm(val ~ time, data=tnf3hr.hf.data)

plot(tnf3hr.hf.data)

lines(0:6,exp(predict(glm.fit,newdata = data.frame(time=0:6))))
#http://stats.stackexchange.com/questions/160253/compare-two-lm-where-one-is-calculated-on-a-subset-of-the-data
tnf3hr.hf.data.2 <- rbind(data.frame(tnf3hr.hf.data,group="A"),
                          data.frame(tnf3hr.hf.data[1:3,],group="B"))

lmtest <- lm(val ~ time*group,data = tnf3hr.hf.data.2)
summary(lmtest)
log(2)/(-tnf3hr.hf.fit$coefficients[2])
summary(tnf3hr.hf.fit)
summary(tnf3hr.hf.fit.2)
log(2)/(-tnf3hr.hf.fit.2$coefficients[2])*60

anova(c(tnf3hr.hf.fit,tnf3hr.hf.fit.2))


# inconsistency found between mannual computing and the browser 


####metabolite correlation per row with sites

sites<-read.delim(file="~/R24/finalsitelvl/merged_sites_1597changes.txt",header = T,stringsAsFactors = F,row.names = 1)
sites[is.na(sites)]<-0
head(sites)

ac.sites<-sites[,1:10]
suc.sites<-sites[,11:20]

??pca
plotPCA(ac.sites)
head(sites)
?prcomp
ac.pca.s <- prcomp(ac.sites,
                 center = TRUE,
                 scale. = TRUE) 
ac.pca <- prcomp(ac.sites,
                   center = TRUE,
                   scale. = FALSE) 
plot(ac.pca)
summary(ac.pca)
summary(ac.pca.s)
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)



acar<-read.delim(file="~/R24/metabolites/filtered.ac.FC.txt",header = T,stringsAsFactors = F)

acar[is.na(acar)]<-0
head(acar)
?cor
covar.ac.acar<-diag(cov(t(acar),t(ac.sites)))

test1<-cor.test(t(acar)[,1],t(ac.sites)[,1])
test1$parameter
?cor.test
#### function to compute the correlation between each row in the matrix with each metabolite

metabolite.cor=function(metabolites=acar,sites=ac.sites){
  head(acar)
  head(sites)
  metabolites[is.na(metabolites)]<-0
  n.compare<-nrow(metabolites)
  n.sites<-nrow(sites)
  sites[is.na(sites)]<-0
  outnames<-rownames(metabolites)
  outnames
  #i=1
  #j=3
  rownames(sites)[j]
  
  outlist.pos<-list()
  outlist.neg<-list()
  corval<-list()
  for(i in 1:n.compare){
    print(outnames[i])
    for(j in 1:n.sites){
      #print(j)
      ### keep all p values for q correction 
      ### store only if pvalue<0.05
      ### 
      if(sum(as.numeric(sites[j,]))!=0 & sum(as.numeric(metabolites[i,]))!=0){
        outlist.pos[[outnames[i]]][[rownames(sites)[j]]]<-cor.test(as.numeric(sites[j,]),as.numeric(metabolites[i,]),method = "pearson")$p.value
        #outlist.neg[[outnames[i]]][[rownames(sites)[j]]]<-cor.test(as.numeric(sites[j,]),-1*as.numeric(metabolites[i,]),method = "pearson")$p.value
        
        #if(tmp<=0.05){
        #  outlist[[outnames[i]]][[rownames(sites)[j]]]<-tmp
        #}
      }

      #outlist[[outnames[i]]][[rownames(sites)[j]]]<-cor.test(as.numeric(sites[j,]),as.numeric(metabolites[i,]),method = "spearman")$p.value
    }
    #why<-sapply(1:nrow(sites),function(j) cor.test(as.numeric(sites[j,]),as.numeric(metabolites[i,])))
    #warnings()
  }
  
  hist(unlist(outlist))
  qvalue(unlist(outlist))
  hist(qvalue(unlist(outlist)),xlim=c(0,0.1))
  tmp.q.pos<-qvalue(unlist(outlist.pos))$qvalue
  tmp.q.neg<-qvalue(unlist(outlist.neg))$qvalue
  
  tmp.q.pos[tmp.q.pos<=0.01]
  tmp.q.neg[tmp.q.neg<=0.01]
  
  ac.sites["Q9R0H0_643",]
  metabolites["C2",]
  ac.sites["Q61425_192",]
  metabolites["C5",]
  ac.sites["P24270_449",]
  metabolites["C4.Ci4",]
  ac.sites["P12710_46",]
  ac.sites["P31786_51",]
  
  
  
  metabolites["C4.Ci4",]
  
  
  setwd("C:/Users/jmeyer/Documents/R24/metabolites/correlations")
  write.table(tmp.q[tmp.q<=0.05],file="AC.pos.q0pt05.txt",sep="\t")
  ?qvalue
  ### return some structure that gives the significant correlations
  ### q-value correct these
  metabolites[1,]
  min(na.omit(outlist[[1]]))
  library(qvalue)
  qvalue(na.omit(outlist[[1]]))$qvalue[qvalue(na.omit(outlist[[1]]))$qvalue<=0.01]
  na.omit(outlist[[1]])[na.omit(outlist[[1]])<=0.01]
  
  
  ac.sites["Q9R0H0_643",]
  ac.sites["P62806_9_13_17",]
  acar[1,]
  
  min(outlist[[1]])
  
  }


sapply(1:nrow(df1), function(i) cor(df1[i,], df2[i,]))

covar.ac.acar
head(acar)



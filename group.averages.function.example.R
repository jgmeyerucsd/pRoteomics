
#take group averages, plot scaled heatmap, and return list of averages and scaled group averages 
###### first part gets table in correct format - individual measures must be in the column labels

head(t(oa))
head(oa)
oat<-oa[,2:8]
gns<-c("con","fr","gl","hfd","hfd_fr","hfd_gl")
allgroups<-c("con2w","fr2w","gl2w","hfd2w","hfd_fr2w","hfd_gl2w","con10w","fr10w","gl10w","hfd10w","hfd_fr10w","hfd_gl10w")
gns2<-unlist(lapply(FUN=rep,gns,times=8))
gns10<-unlist(lapply(FUN=rep,gns,times=6))
group.names<-c(gns2,gns10)
row.names(oat)<-oa[,1]
row.names(oat)<-paste(group.names,row.names(oat),sep="")
head(t(oat),8)




group.aves=function(table=t(oat),groups=allgroups){
  if(length(groups)==1){
    groups<-unique(substr(names(table),start=2,stop=nchar(names(table))-2))
    }
  averaged<-data.frame(rep(0,times=nrow(table)))
  rownames(averaged)<-rownames(table)
  i=1
  for(x in groups){
    tmp.grp<-grep(colnames(table),pattern=x)
    if(i==1){
      averaged[,1]<-rowMeans(table[,tmp.grp])
      names(averaged)[i]<-groups[i]
    }
    if(i>=2){
      averaged<-cbind(averaged,rowMeans(table[,tmp.grp]))
      names(averaged)[i]<-groups[i]
    }
    i=i+1
  }
  head(averaged)
  scaled<-t(scale(t(averaged)))
  head(scaled)
  sd(scaled[1,])
  library(gplots)
  library(RColorBrewer)
  mycolors<-colorRampPalette(c("blue","white","red"))(n=50)
  nrow(scaled)
  ?heatmap.2
  heatmap.2(scaled,scale="none",col = mycolors,trace="none",Colv = F,Rowv = F,dendrogram = "none",main=paste(nrow(scaled),"rows, scaled"),margins = c(10,10),tracecol = "black",key.title = "")
  
  return(list(ave=averaged,scaled.ave=scaled))
}




group.aves(table=excludecontam)->testave
names(testave)

testave[["scaled.ave"]][1:100,]

oa.scaled<-group.aves(table=t(oat),groups=allgroups)
oa.scaled

### also for the other metabolites
head(t(aa))
head(aa)
aat<-aa[,1:ncol(aa)]

gns<-c("con","fr","gl","hfd","hfd_fr","hfd_gl")
allgroups<-c("con2w","fr2w","gl2w","hfd2w","hfd_fr2w","hfd_gl2w","con10w","fr10w","gl10w","hfd10w","hfd_fr10w","hfd_gl10w")
gns2<-unlist(lapply(FUN=rep,gns,times=8))
gns10<-unlist(lapply(FUN=rep,gns,times=6))
group.names<-c(gns2,gns10)
row.names(aat)
row.names(aat)<-paste(group.names,row.names(aat),sep="")
head(t(aat),8)

### get group averages
aa.scaled<-group.aves(table=t(aat),groups=allgroups)

head(aa.scaled)

write.table(oa.scaled[])

head(ac)

act<-ac[,1:ncol(ac)]
row.names(act)
row.names(aat)
row.names(act)<-paste(group.names,row.names(act),sep="")
head(act,8)
colnames(act)<-paste("ac",colnames(act),sep="")

ac.scaled
ac.scaled<-group.aves(table=t(act),groups=allgroups)


acoa<-read.delim(file="~/R24/metabolite_corr/acylcoA.txt",header = T, row.names = 1,stringsAsFactors = F)
head(acoa)


row.names(acoa)<-paste(group.names,row.names(acoa),sep = "")
colnames(acoa)
colnames(acoa)<-paste("coA.",colnames(acoa),sep="")

acoa.scaled<-group.aves(table=t(acoa),groups=allgroups)

allmetab<-rbind(
  acoa.scaled[[2]],
  oa.scaled[[2]],
  aa.scaled[[2]],
  ac.scaled[[2]]
)


write.table(sep="\t",x=allmetab,file="all.metabolites.zscaled.txt",quote = F, row.names = T,col.names = T)

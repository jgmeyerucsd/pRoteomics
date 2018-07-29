
#take group averages, plot scaled heatmap, and return list of averages and scaled group averages 
###### first part gets table in correct format - individual measures must be in the column labels
list.files("~/R24/metabolites/")

aaoa<-read.delim(file="~/R24/metabolites/acylcoA.txt",header = T, row.names = 1,stringsAsFactors = F)
aaoa<-read.delim(file="~/R24/metabolites/AC.txt",header = T, row.names = 1,stringsAsFactors = F)

#separate aa and oa

#aaoa<-read.delim(file="~/R24/metabolites/OA.txt",header = T, row.names = 1,stringsAsFactors = F)

head(aaoa)
aaoa=aaoa+0.0000001

#aaoat<-aaoa[,2:22]
gns<-c("con","fr","gl","hfd","hfd_fr","hfd_gl")
allgroups<-c("con2w","fr2w","gl2w","hfd2w","hfd_fr2w","hfd_gl2w","con10w","fr10w","gl10w","hfd10w","hfd_fr10w","hfd_gl10w")
gns2<-unlist(lapply(FUN=rep,gns,times=8))
gns10<-unlist(lapply(FUN=rep,gns,times=6))
group.names<-c(gns2,gns10)





my.t.test.pv <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.t.test.est <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$estimate)
}


row.def<-list(1:8,9:16,17:24,25:32,33:40,41:48,49:54,55:60,61:66,67:72,73:78,79:84)
head(aaoa)

multi.test3=function(input=aaoa,row.def=row.def){
  n.col<-ncol(input)
  head(input)
  pval<-c()
  fc<-c()
  all.pv<-data.frame()
  all.fc<-data.frame()
  #input[,i]
  for(i in 1:n.col){
    pval<-c()
    fc<-c()
    for(j in 2:6){
      pval[j-1]<-my.t.test.pv(input[row.def[[1]],i],input[row.def[[j]],i])
      fc[j-1]<-log(my.t.test.est(input[row.def[[1]],i],input[row.def[[j]],i])[2]/my.t.test.est(input[row.def[[1]],i],input[row.def[[j]],i])[1],base=2)
    }
    for(j in 8:12){
      pval[j-2]<-my.t.test.pv(input[row.def[[7]],i],input[row.def[[j]],i])
      fc[j-2]<-log(my.t.test.est(input[row.def[[7]],i],input[row.def[[j]],i])[2]/my.t.test.est(input[row.def[[7]],i],input[row.def[[j]],i])[1],base=2)
    }
    all.pv<-rbind(all.pv,pval)
    rownames(all.pv)[i]<-paste(colnames(input)[i])
    all.fc<-rbind(all.fc,fc)
    rownames(all.fc)[i]<-paste(colnames(input)[i])
  }
  colnames(all.fc)<-paste(1:10)
  colnames(all.pv)<-paste(1:10)
  all<-list(all.fc=all.fc,all.pv=all.pv)
  return(all)
}

all<-multi.test3(input=aaoa,row.def=row.def)
warnings()
library(qvalue)

all.fc<-all$all.fc
all.pv<-all$all.pv
all.qmat.aaoa<-matrix(qvalue(p=unlist(all$all.pv))$qvalues,byrow = FALSE,ncol=ncol(all.pv))
hist(qvalue(p=unlist(all$all.pv)))
colnames(all.qmat.aaoa)<-colnames(all$all.pv)
colnames(all.fc)<-paste(c(rep("2w",times=5),rep("10w",times=5)),gns[2:6])
row.names(all.qmat.aaoa)<-row.names(all$all.pv)
head(all.qmat.aaoa)

### those with 0.01
q01<-all.fc
q05<-all.fc

#all.fc[all.qmat.aaoa>0.01]<-NA

q05[all.qmat.aaoa>0]<-""
q05[all.qmat.aaoa<0.05]<-"*"
q05[all.qmat.aaoa<0.01]<-"**"

all.fc[all.qmat.aaoa>0.05]<-NA
#all.qmat.aaoa
#all.fc[abs(all.fc)<0.58]<-NA

### heatmap of acyl coA changes
?heatmap.2
library(ggplot2)
library(gplots)
my_palette<-colorRampPalette(c("blue","white","red"))(n=50)

heatmap.2(as.matrix(all.fc),trace="none",scale="none",Colv = F,Rowv=F,
          na.col="grey",dendrogram = "none",col=my_palette,tracecol = "black",
          key.title = "",
          margins = c(10,10), colsep = 5,
          cellnote=q05,  notecol="black")

head(all.fc)


fnames(testave)

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
heatmap.2(allmetab,scale="none",col = mycolors,trace="none",Colv = F,Rowv = F,dendrogram = "none",main=paste(nrow(scaled),"rows, scaled"),margins = c(10,10),tracecol = "black",key.title = "")


write.table(sep="\t",x=allmetab,file="all.metabolites.zscaled.txt",quote = F, row.names = T,col.names = T)



### now read in the raw signals for the acetyl/succinyl
getwd()
setwd("~/R24/metabolite_corr/")
list.files()

ack<-read.delim("acK_sitelvl_signals.txt",stringsAsFactors = F)

head(ack)

colnames(ack)
rownames(ack)<-ack[,1]

sitegrp<-c("con2w","fr2w","gl2w","hfd2w","hfd_fr2w","hfd_gl2w","con10w","fr10w","gl10w","hfd10w","hfd_fr10w","hfd_gl10w")
sgn<-unlist(lapply(FUN=rep,sitegrp,times=5))
colnames(ack)[2:60]<-sgn[1:59]
acK<-ack[,2:60]

ack.scaled<-group.aves(table=acK,groups=allgroups)

ack.scaled[[2]]->ack.z
3.63349245





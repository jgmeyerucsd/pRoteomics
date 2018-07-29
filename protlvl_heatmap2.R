prot<-read.delim(file="C:/users/jmeyer/documents/R24/20161021_R24_protlvl_candidates_all.xls")


head(prot)
library(reshape2)
prot.castFC <- dcast(data = prot, formula = Group ~ Comparison..group1.group2.,value.var=c("AVG.Log.Ratio"))
head(prot.castFC)
names(prot.castFC)
prot.castFDR <- dcast(data = prot, formula = Group ~ Comparison..group1.group2.,value.var=c("Qvalue"))
head(prot.castFDR)

order<-c(57,60,62,67,66,2,3,4,5,6)
#### invert the labels on the upside down rows
prot.FC.fix1<-prot.castFC[,c(57,60,62,67,66,2,3,4,5,6)]
head(prot.FC.fix1)
prot.FC.fix2<-prot.FC.fix1
prot.FC.fix2[,c(4,6:10)]<-prot.FC.fix2[,c(4,6:10)]*-1
head(prot.FC.fix2)
head(prot.castFDR[,order])
rownames(prot.FC.fix2)<-substr(prot.castFC[,1],start=4,stop=9)

#### make true/false matrix of significant values
prot.sig<-prot.castFDR[,order]<=0.05
#prot.sig<-prot.sig[,order]>=0.01
#### template matrix of values I want so I don't overwrite other datmat in case of errors
prot.onlysig<-prot.FC.fix2

#### replace non-significant values in the matrix with NA
prot.onlysig[prot.sig==FALSE]<-0
prot.onlysig[abs(prot.onlysig)<=0.58]<-0
prot.onlysig[is.na(prot.onlysig)]<-0
prot.dist<-dist(prot.onlysig)
prot.clust<-hclust(prot.dist)
#abs(-1)
head(prot.onlysig)
prot.onlysig
### put names back
#row.names(kac.FC.order)<-as.character(kac.castFC[,1])
nrow(prot.onlysig)
prot.filt<-prot.onlysig[rowSums(prot.onlysig==0)<10,]
nrow(prot.filt)
prot.filt.dist<-dist(prot.filt)
prot.filt.clust<-hclust(prot.filt.dist)
prot.filt[prot.filt==0]<-NA
nrow(prot.filt)
### plot heatmap
library(gplots)
heatmap.2(as.matrix(prot.filt),main="protein level, 420 proteins",scale="none", 
          key.par=list(),trace="none",  keysize = 2,
          key.title = "",tracecol="black",Rowv=as.dendrogram(prot.filt.clust), 
          Colv=F, na.col="grey",dendrogram="none",col=bluered)
dev.off()
write.table()
### try plotting only those with suk/ack sites
head(suk.castFC)
acsuk<-c(substr(suk[,1],start=1,stop=6),substr(kac[,1],start=1,stop=6))
mod.prot.names<-unique(acsuk)
match("Q3U186",sig.protnames)
match("Q3U186",mod.prot.names)

sig.protnames<-row.names(prot.onlysig)
match(sig.protnames,mod.prot.names)

mod.rows<-match(mod.prot.names,sig.protnames)


"Q3U186"
length(na.omit(match(sig.protnames,mod.prot.names)))
modmat<-as.matrix(prot.onlysig[mod.rows,])
prot.onlysig[mod.rows[1],]
modmat[is.na(modmat)]<-0
nrow(modmat)
head(modmat)
prot.onlysig["O08756",]  ### 
which(row.names(prot.onlysig)=="O08756") ### row 18

which(mod.prot.names=="Q64459") ### 6


### plot heatmap - only modified proteins
heatmap.2(modmat,scale="none", trace="none", Rowv=T, Colv=F, na.col="grey",dendrogram="none",col=bluered)
heatmap.2(filt,scale="none", trace="none", Rowv=T, Colv=F, na.col="grey",dendrogram="none",col=bluered)

nrow(modmat)
?dist
final<-as.matrix(prot.onlysig[mod.rows[c(1:10)],])
submat[complete.cases(submat),]
filt<-modmat[rowSums(modmat==0)!=10,]
nrow(filt)
filt.dist<-dist(filt)
filt.clust<-hclust(filt.dist)
filt[filt==0]<-NA
clara(filt,k=1)
daisy(filt)
hclust(daisy(as.matrix(prot.onlysig[mod.rows,])))


### plot heatmap -- any change
heatmap.2(filt,scale="none", trace="none", Rowv=as.dendrogram(filt.clust), Colv=F, na.col="grey",dendrogram="none",col=bluered,main="107 modified proteins")

prot.FC.fix2
?hist
?heatmap.2
dist(as.matrix(prot.FC.fix2[na.omit(match(mod.prot.names,row.names(prot.FC.fix2))),]))
hclust(dist(as.matrix(prot.FC.fix2[na.omit(match(mod.prot.names,row.names(prot.FC.fix2))),])))

heatmap.2(as.matrix(prot.FC.fix2[na.omit(match(mod.prot.names,row.names(prot.FC.fix2))),]),scale="none", trace="none", Rowv=T, Colv=F, na.col="grey",col=bluered)
heatmap(as.matrix(prot.FC.fix2[na.omit(match(mod.prot.names,row.names(prot.FC.fix2))),]),scale="none",na.col="grey",col=bluered)



#### write table of proteins with at least one change above 1.5 fold
write.table(file="Qpt05_fc1pt5.txt",sep="\t",filt)
write.table(file="protlvl_420changes_Q05_fc58.txt",sep="\t",prot.filt)
setwd("~/R24/protlvlreport/")
getwd()





#### prot level heatmap of top 2 MCODE clusters
c1<-read.csv(file="cluster1.csv",stringsAsFactors = F,header=T,row.names = 1)
head(c1)
c2<-read.csv(file="cluster2.csv",stringsAsFactors = F,header=T,row.names = 1)
head(c2)

c1[is.na(c1)]<-0
c1
c1.dist<-dist(c1)
c1.clust<-hclust(c1.dist)
c1[c1==0]<-NA
head(c1)
heatmap.2(as.matrix(c1),Rowv = F, Colv = F)
heatmap.2(as.matrix(c1),main="protein level, top 2 clusters",scale="none", 
          key.par=list(),trace="none",  keysize = 2,
          key.title = "",tracecol="black",Rowv=as.dendrogram(c1.clust), 
          Colv=F, na.col="grey",dendrogram="none",col=bluered)

### clust 2 heatmap
c2[is.na(c2)]<-0
c1
c2.dist<-dist(c2)
c2.clust<-hclust(c2.dist)
c2[c2==0]<-NA
head(c1)
heatmap.2(as.matrix(c2),Rowv = F, Colv = F)
heatmap.2(as.matrix(c2),main="protein level, top 2 clusters",scale="none", 
          key.par=list(),trace="none",  keysize = 2,
          key.title = "",tracecol="black",Rowv=as.dendrogram(c2.clust), 
          Colv=F, na.col="grey",dendrogram="none",col=bluered)

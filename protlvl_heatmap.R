
### new quant

prot<-read.delim(file="Z:/R24/resource/proteinlvl/protlvl_spectronaut_qpt01filt.txt")


head(prot)
library(reshape2)
prot.castFC <- dcast(data = prot, formula = Group ~ Comparison..group1.group2.,value.var=c("AVG.Log2.Ratio"))
head(prot.castFC)
names(prot.castFC)
prot.castFDR <- dcast(data = prot, formula = Group ~ Comparison..group1.group2.,value.var=c("Qvalue"))
head(prot.castFDR)
names(prot.castFC)
names(prot.castFDR)

order<-c(56,59,61,66,65,2,3,4,5,6)
order<-c(56,59,61,66,65,2,3,4,5,6,22)  ### including the 10wGluc/10wHFD

#order<-c(2,3,4,5,6)
prot.FC1<-prot.castFC[,order]
head(prot.FC1)

rownames(prot.FC1)<-substr(prot.castFC[,1],start=1,stop=6)
rownames(prot.FC1)
prot.FDR1<-prot.castFDR[,order]
rownames(prot.FDR1)<-prot.castFDR[,1]

#### invert the labels on the upside down rows
prot.FC2<-prot.FC1

prot.FC2[,c(4,6:10)]<-prot.FC2[,c(4,6:10)]*-1
head(prot.FC1)
head(prot.FC2)

#prot.sig<-prot.sig[,order]>=0.01
#### template matrix of values I want so I don't overwrite other datmat in case of errors
prot.onlysig<-prot.FC2
prot.FC2

#### replace non-significant values in the matrix with NA
prot.onlysig[is.na(prot.FC2)]<-0



### put names back
#row.names(kac.FC.order)<-as.character(kac.castFC[,1])
nrow(prot.onlysig)
#rownames(prot.filt)
prot.filt<-prot.onlysig[rowSums(prot.onlysig==0)<11,]
nrow(prot.filt)
prot.filt.dist<-dist(prot.filt)
prot.filt.clust<-hclust(prot.filt.dist)
prot.filt[prot.filt==0]<-NA
nrow(prot.filt)
head(prot.filt,10)
write.table(file="Z:/R24/resource/proteinlvl/254candidates.onlysig.HFDvGluc.txt",sep="\t",prot.filt)
#### start here #####

read.delim(file="Z:/R24/resource/proteinlvl/254candidates.onlysig.txt",header = T)->prot.filt



### plot heatmap
library(RColorBrewer)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=11)
library(gplots)
?heatmap.2

hm<-heatmap.2(as.matrix(prot.filt),main="protein level, 254 proteins",scale="none", 
          key.par=list(),trace="none",  keysize = 2,
          key.title = "",tracecol="black",Rowv=as.dendrogram(prot.filt.clust), 
          Colv=F, na.col="grey90",dendrogram="none",col=my_palette,colsep=c(0,5,10),sepcolor = 1,margins = c(5,15))
dev.off()
hm$rowInd

write.table(file="Z:/R24/resource/proteinlvl/254candidates.hmorder.onlysig.txt",sep="\t",prot.filt[hm$rowInd,])

uniprot<-rownames(prot.filt)


getwd()
go<-read.delim("Z:/R24/resource/proteinlvl/HFDonly_goMF.txt",stringsAsFactors = F)

head(go)
plotme<-cbind(go[,"Term"],go[,"Benjamini"])
?barplot
barplot(-log(go[1:14,"Benjamini"],base=10),horiz = F, ylim=c(0,10))

dev.off()


#### add gene names
source("http://bioconductor.org/biocLite.R")
biocLite("UniProt.ws")
biocLite("S4Vectors")
biocLite("IRanges")

library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]
genemap <- select(up,uniprot, "GENES")
genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]
geneid.filt <- genes[match(rownames(prot.filt), genemap[,"UNIPROTKB"])]


#### get only the cyp - p450s
cypset<-grep("Cyp",geneid)
rownames(prot.filt)<-geneid
geneid
rownames()
row.names(cypset)



hm.ccyp<-heatmap.2(as.matrix(prot.filt[cypset,]),main="protein level, 254 proteins",scale="none", 
              key.par=list(),trace="none",  keysize = 2,
              key.title = "",tracecol="black",Rowv=F, 
              Colv=F, na.col="white",dendrogram="none",col=my_palette,colsep=c(0,5,10),sepcolor = 1,margins = c(5,15))


cypset<-grep("Cyp",geneid)
geneid
translocon<-c("Hspa5","Sec63","Sec61a1","Sec61b","Ssr4","Bcap31")
match(translocon,rownames(prot.filt))

rownames(prot.filt)<-geneid
geneid
rownames()
row.names(cypset)



hm.tl<-heatmap.2(as.matrix(prot.filt[match(translocon,rownames(prot.filt)),]),main="protein level, 254 proteins",scale="none", 
                   key.par=list(),trace="none",  keysize = 2,
                   key.title = "",tracecol="black",Rowv=F, 
                   Colv=F, na.col="white",dendrogram="none",col=my_palette,colsep=c(0,5,10),sepcolor = 1,margins = c(5,15))



length(geneid)
length(uniprot)

pf.genes<-cbind(geneid.filt,prot.filt)
nchar(geneid)
head(pf.genes)
pf.genes.narm<-pf.genes[which(is.na(pf.genes[,1])==FALSE),]
nrow(pf.genes.narm)
head(pf.genes.narm)
pf<-cbind(pf.genes.narm[,2:7])
rownames(pf)<-pf.genes.narm[,1]

### generate barplot of proteins increased or decreased
nrow(prot.filt)
geneid
rownames(prot.filt)
na.omit(geneid.filt)

nrow(pf)




count.updown=function(datmat=pf){
  ncols<-ncol(datmat) 
  changes.list<-list()

  group<-c()
  x<-c()
  y.up<-c()
  y.down<-c()
  
  for(i in 1:ncols){
    x<-c(x,i)
    y.up<-c(y.up,length(na.omit(datmat[,i][datmat[,i]>0])))
    y.down<-c(y.down,length(na.omit(datmat[,i][datmat[,i]<0])))
  }
  group<-c(rep("up",times=ncols),rep("down",times=ncols))
  y<-c(y.up,-y.down)
  dat<-data.frame(group=group,
                  x=x,
                  y=y,
    stringsAsFactors=F)
  
  return(dat)
  
}
count.updown()->dat

library(ggplot2)
sp<-ggplot(dat, aes(x=x, y=y, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual("legend", values = c("up" = "red", "down" = "blue")) +
  coord_flip() +   scale_x_reverse() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.y=element_blank()) + 
  ylab("number of changes")  

sp+scale_y_continuous(limits=c(-220,60),breaks=c(-200,-100, 0, 50))


?theme
dat
?element_line

###### volcano plots
par(cex=1.5)
plot(x=tempset[,"log2FC"],
     y=tempset[,"log_oddsDE"],
     pch=20,
     main=paste(x),
     xlab="log2(Fold Change)", ylab="odds differential expression")
points(x=tempset.filt[,"log2FC"],y=tempset.filt[,"log_oddsDE"],pch=20,col="red")

source("http://bioconductor.org/biocLite.R")
biocLite("mygene")




library(mygene)

res <- queryMany(geneid, scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')

res$go.BP
res$go.MF
write.table(res$go.MF,file="MF.txt")
write.table(mfunlist,file="MF.txt")

mfunlist<-unlist(res$go.MF)
bpunlist<-unlist(res$go.BP)

pmatch(table=mfunlist,x="protein")
mfunlist[grepl(mfunlist,pattern="steroid")]
mfunlist[grepl(mfunlist,pattern="oxida")]
bpunlist[grepl(bpunlist,pattern="steroid")]
bpunlist[grepl(bpunlist,pattern="oxida")]

mitofun<-regexpr(pattern="mito",mfunlist)

mitofun[na.omit(mitofun==4)]

mfunlist[87]

barplot()

datmat[,x]

))))}}}}}


getwd()
write.table(file="~/R24/protlvl/MQlibQcompl_clustered.txt",sep="\t",prot.filt[hm$rowInd,])

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

getwd()



setwd("C:/BAT/")
list.files()->files
files


### read in acetyl mapDIA output
sites<-read.delim(file="BAT_analysis_output.txt")
head(sites)
library(reshape2)

### cast acetyl
castFC <- dcast(data = sites, formula = Protein ~ Label2,value.var=c("log2FC"))
castFDR <- dcast(data = sites, formula = Protein ~ Label2,value.var=c("FDR"))





### what is the order of columns I want? manually...
names(castFC)
names(castFC)

### new orders for the columns
#### order: RT sirt5 KO, CA sirt5 KO
order<-c(5,3, 2,4)


### order the succinyl sites as desired and name the rows
FDR.order<-castFDR[,order]
FC.order<-castFC[,order]
row.names(FDR.order)<-as.character(castFDR[,1])
row.names(FC.order)<-as.character(castFC[,1])

head(FDR.order)

### filter suk by FDR <0.01
significant<-castFDR[,order]<=0.01

onlysig<-FC.order

#### replace non-significant values in the matrix with NA
onlysig[significant==FALSE]<-NA
head(onlysig)
#### replace the values not greater than 1 from the matrix
onlysig[abs(onlysig)<1]<-NA
onlysig<-onlysig[rowSums(is.na(onlysig))<4,]


onlysig.4clust<-onlysig
onlysig.4clust[is.na(onlysig)]<-0
head(onlysig.4clust)
dist<-dist(onlysig.4clust)
clust<-hclust(dist)

rownames(onlysig)



library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=11)


heatmap.2(as.matrix(onlysig),scale="none", trace="none", tracecol = "black", Rowv=as.dendrogram(clust), Colv=T, na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)
nrow(onlysig)
s.index.unfilt<-grep(row.names(onlysig),pattern="K100")

succinyl.rows<-onlysig[s.index.unfilt,]
succinyl.raw<-row.names(succinyl.rows)
s.p<-substr(succinyl.raw,start=1, stop=6)
unique(s.p)[1]->x
site.df<-cbind(s.p, succinyl.raw)

ss<-c()
sites.list<-list()



for(x in unique(s.p)){
  tmp.position<-which(s.p==x)
  tmp.sites<-unique(unlist(strsplit(succinyl.raw[tmp.position],split=paste(x,"_",collapse="",sep=""))))
  tmp.sites<-tmp.sites[nchar(tmp.sites)>1]
  tmp.sites<-unlist(strsplit(tmp.sites,split="_"))
  tmp.sites<-unique(tmp.sites)
  ss<-c(ss,length(grep(tmp.sites,pattern="K100")))
  sites.list[[x]]<-tmp.sites
}

### how many total succinyl sites
sum(ss)
### how many total succinyl proteins
length(sites.list)
??venn
venn(list())
tmp.prot

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

onlysig[onlysig>-log(base=2, 2)]<-NA
onlysig<-onlysig[rowSums(is.na(onlysig))<4,]


head(onlysig)
nrow(onlysig)
onlysig.4clust<-onlysig
onlysig.4clust[is.na(onlysig)]<-0
head(onlysig.4clust)
dist<-dist(onlysig.4clust)
clust<-hclust(dist)

rownames(onlysig)



library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("white", "red"))(n=11)


heatmap.2(as.matrix(onlysig),scale="none", trace="none", tracecol = "black", Rowv=as.dendrogram(clust), Colv=F, na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)
nrow(onlysig)
s.index.unfilt<-grep(row.names(onlysig),pattern="K100")

onlysig[1,1:2]
onlysig[1,3:4]

succinyl.rows.5bko<-onlysig[s.index.unfilt,1:2]
succinyl.rows.ca<-onlysig[s.index.unfilt,3:4]

sr5<-succinyl.rows.5bko[rowSums(is.na(succinyl.rows.5bko))<2,]
src<-succinyl.rows.ca[rowSums(is.na(succinyl.rows.ca))<2,]


onlysig.4clust<-src
onlysig.4clust[is.na(src)]<-0
head(onlysig.4clust)
dist<-dist(onlysig.4clust)
clust<-hclust(dist)

heatmap.2(as.matrix(src),scale="none", trace="none", tracecol = "black", Rowv=as.dendrogram(clust), Colv=F,  na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)

onlysig.4clust<-sr5
onlysig.4clust[is.na(sr5)]<-0
head(onlysig.4clust)
dist<-dist(onlysig.4clust)
clust<-hclust(dist)

heatmap.2(as.matrix(sr5),scale="none", trace="none", tracecol = "black", Rowv=as.dendrogram(clust), Colv=F,  na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)

nrow(sr5)
nrow(src)
succinyl.rows.5bko

succinyl.raw.s5<-row.names(sr5)
succinyl.raw.ca<-row.names(src)

s.p.s5<-substr(succinyl.raw.s5, start=1, stop=6)
s.p.ca<-substr(succinyl.raw.ca, start=1, stop=6)

#unique(s.p)[1]->x


ss.s5<-c()
ss.ca<-c()

sites.list.s5<-list()
sites.list.ca<-list()

s.p<-s.p.s5
succinyl.raw<-succinyl.raw.s5


for(x in unique(s.p)){
  tmp.position<-which(s.p==x)
  tmp.sites<-unique(unlist(strsplit(succinyl.raw[tmp.position],split=paste(x,"_",collapse="",sep=""))))
  tmp.sites<-tmp.sites[nchar(tmp.sites)>1]
  tmp.sites<-unlist(strsplit(tmp.sites,split="_"))
  tmp.sites<-unique(tmp.sites)
  ss.s5<-c(ss.s5,length(grep(tmp.sites,pattern="K100")))
  sites.list.s5[[x]]<-tmp.sites
}


s.p<-s.p.ca
succinyl.raw<-succinyl.raw.ca

for(x in unique(s.p)){
  tmp.position<-which(s.p==x)
  tmp.sites<-unique(unlist(strsplit(succinyl.raw[tmp.position],split=paste(x,"_",collapse="",sep=""))))
  tmp.sites<-tmp.sites[nchar(tmp.sites)>1]
  tmp.sites<-unlist(strsplit(tmp.sites,split="_"))
  tmp.sites<-unique(tmp.sites)
  ss.ca<-c(ss.ca,length(grep(tmp.sites,pattern="K100")))
  sites.list.ca[[x]]<-tmp.sites
}

### how many total succinyl sites
sum(ss.s5)
sum(ss.ca)

### how many total succinyl proteins
length(sites.list.s5)
length(sites.list.ca)

sites.list.s5.down<-sites.list.s5
sites.list.ca.down<-sites.list.ca

sites.list.s5[1:10]



sites.list.s5[1:10]

### enumerate the values of the lists and append to their list entry name

list<-sites.list.s5
list<-sites.list.ca

expanded<-c()
names(list)
for(x in names(list)){
  for(y in list[[x]]){
    expanded<-c(expanded,paste(x, y, sep="_",collapse = ""))
  }
}

s5.expanded<-expanded
ca.expanded<-expanded
length(s5.expanded)
length(ca.expanded)

s5.expanded<-s5.expanded[grep(s5.expanded, pattern="K100")]
ca.expanded<-ca.expanded[grep(ca.expanded, pattern="K100")]


?setdiff

s5.ca.intersect<-intersect(ca.expanded, s5.expanded)
s5.ca.intersect
s5.ca.intersect[grep( s5.ca.intersect,pattern= "K100")]

unique(substr(s5.ca.intersect,start=1, stop=6))
write.table(s5.ca.intersect, "increasedsites.sirt5bko.ca.overlap.txt", quote=F, row.names = F, col.names = F)
?venn
library(VennDiagram)
??VennDiagram
venn(list("cold-acclimated"=ca.expanded, "Sirt5"=s5.expanded))
dev.off()
draw.pairwise.venn(60,657,17, ext.text=F)

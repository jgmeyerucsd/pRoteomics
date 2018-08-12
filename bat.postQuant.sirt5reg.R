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
names(castFC)[order]

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
head(onlysig)

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



### remove rows without at least one value
naomit<-onlysig[rowSums(is.na(onlysig))<4,]

#### split into acetyl and succinyl
s.index<-grep(row.names(naomit),pattern="K100")
a.index<-grep(row.names(naomit),pattern="K42")



both.index<-intersect(s.index,a.index)
?intersect
?rm
a.index!=both.index
only.a.index<-setdiff(a.index, both.index)

s.set<-naomit[s.index,]
a.set<-naomit[only.a.index,]
nrow(s.set)

unique(substr(rownames(s.set),start=1,stop=6))

head(naomit)
s.naomit4clust<-s.set
s.naomit4clust[is.na(s.naomit4clust)]<-0
s.nadist<-dist(s.naomit4clust)
s.naclust<-hclust(s.nadist)
head(s.set)
dev.off()
hm1<-heatmap.2(as.matrix(s.set),scale="none", trace="none", tracecol = "black", 
               Rowv=as.dendrogram(s.naclust), Colv=T, 
               na.col="grey",dendrogram="none",
               col=my_palette,cexCol=0.8,labRow = F,key.title = "",
               cexRow = 0.8)


write.table(s.set[hm1$rowInd,],quote=F,sep="\t",file="succinylonly_heatmapOrder.txt")

s.5bko.set<-s.set[,1:2]
s.5bko.set.filtered<-s.5bko.set
s.5bko.set.filtered[s.5bko.set.filtered<0]<-NA
s.5bko.set.filtered<-s.5bko.set.filtered[rowSums(is.na(s.5bko.set.filtered))<2,]



nrow(s.5bko.set.filtered)
s.5bko.set.filtered->s5
s5[is.na(s5)]<-0
s5dist<-dist(s5)
s5clust<-hclust(s5dist)



hm2<-heatmap.2(as.matrix(s.5bko.set.filtered),scale="none", trace="none", tracecol = "black", 
               Rowv=as.dendrogram(s5clust), Colv=T, 
               na.col="grey",dendrogram="none",
               col=my_palette,cexCol=0.8,labRow = F,key.title = "",
               cexRow = 0.8)

unique(substr(rownames(s5),start=1,stop=6))

hm1$colInd



ordered.naomit<-naomit[hm1$rowInd,]
heatmap(as.matrix(ordered.naomit),Rowv=F,Colv = F,scale=F,)
hm1<-heatmap.2(as.matrix(ordered.naomit),scale="none", trace="none", tracecol = "black", Rowv=F, Colv=F, na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)

uniprot<-substr(rownames(ordered.naomit),start=1,stop=6)

#### get gene names per row
uniprot<-substr(rownames(s.set),start=1, stop=6)

library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]
genemap <- select(up,unique(uniprot) , "GENES")



genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]



geneids<-toupper(geneid)






write.table(file="BAT_succinyl4hm.txt",sep="\t",cbind(geneids,s.set),quote=F)

#### overlap with gene sets
list.files()

bOx<-read.delim("customTCA.txt",skip = 1,head=T,stringsAsFactors = F)
bOx<-read.delim("betaOx.txt",skip = 1,head=T,stringsAsFactors = F)
bOx<-read.delim("ETC.txt",skip = 1,head=T,stringsAsFactors = F)
bOx<-read.delim("AAmetab.txt",skip = 1,head=T,stringsAsFactors = F)

toupper(bOx)


tmp.index<-is.na(match(geneids,toupper(unlist(bOx))))
tmpset<-s.set[tmp.index==FALSE,]
rownames(tmpset)<-paste(geneids[tmp.index==FALSE],substr(rownames(tmpset),start=8,stop=nchar(rownames(tmpset))))
tca<-tmpset
tca<-tca[order(row.names(tca)),]


library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=31)
#?heatmap.2
m.hm<-heatmap.2(as.matrix(tca),scale="none", trace="none", tracecol="black", colsep=c(5,10,15), Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette,cexRow = 0.4,
                margins = c(15,15))



fbox.aves<-group.aves(tmpset[,31:59], groups=NA)
nrow(box.aves[[2]])

tiff(filename = "customFAbetaOx.tif",
     width = 8, height = 11, units = "in", pointsize = 12,
     compression = c("none"),
     bg = "white", res = 600, family = "", restoreConsole = TRUE,
     type = c("windows"))

my_palette <- colorRampPalette(c("blue","white", "red"))(n = 30)
heatmap.2(as.matrix(box.aves$scaled.ave),Rowv = F,Colv = F,
          na.color = "grey",
          dendrogram = "none",
          trace = "none",
          col=my_palette,
          tracecol = "black",
          cexRow=0.4,
          margins = c(15,15))


dev.off()































nrow(s.set)

#### pull out only sets with succinyl, only sets with acetyl, only sites with both

succk<-naomit[grep("K100",row.names(naomit)),]
ack<-naomit[grepl("K42",row.names(naomit)),]

both<-naomit[grepl("K42",row.names(naomit))&grepl("K100",row.names(naomit)),]
nrow(naomit)

s4c<-succk
s4c[is.na(s4c)]<-0
s4cdist<-dist(s4c)
s4cclust<-hclust(s4cdist)
heatmap.2(as.matrix(succk),scale="none", trace="none", tracecol = "black", Rowv=as.dendrogram(s4cclust), Colv=T, na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)
nrow(succk)
#### acetyl
a4c<-ack
a4c[is.na(a4c)]<-0
a4cdist<-dist(a4c)
a4cclust<-hclust(a4cdist)
nrow(ack)
heatmap.2(as.matrix(ack),scale="none", trace="none", tracecol = "black", Rowv=as.dendrogram(a4cclust), Colv=T, na.col="grey",dendrogram="none",col=my_palette,cexCol=0.8,labRow = F,key.title = "",cexRow = 0.8)

nrow(both)


dev.off()

heatmap.2(as.matrix(suk.onlysig),scale="none", trace="none", Rowv=as.dendrogram(suk.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)

### remove rows with all NA in either set and plot heatmaps
?hist
hist(10-rowSums(is.na(kac.naomit)),breaks=seq(from=0,to=10,by=1),lwd=2,ylab="number of sites",xlab="number of diets",main="diets causing acetyl site changes")

kac.naomit<-kac.onlysig[rowSums(is.na(onlysig))<10,]
naomit<-onlysig[rowSums(is.na(onlysig))==0,]


par(cex=2)
shist<-hist(10-rowSums(is.na(suk.naomit)),lwd=2,breaks=seq(from=0,to=10,by=1),ylab="number of sites",xlab="number of diets",main="diets causing succinyl site changes")
shist
suk.naomit<-suk.onlysig[rowSums(is.na(suk.onlysig))<10,]
suk.naomit<-suk.onlysig[rowSums(is.na(suk.onlysig))<10,]



heatmap.2(as.matrix(kac.onlysig),scale="none", trace="none", Rowv=as.dendrogram(kac.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)
heatmap.2(as.matrix(suk.onlysig),scale="none", trace="none", Rowv=as.dendrogram(suk.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)





###make merged heatmap
?merge
merged<-merge(kac.onlysig,suk.onlysig,all=T,by=0)
#testmerge<-merge(kac.onlysig[1:5,],suk.onlysig[1:5,],all=T,by=0)
#testmerge
rownames(testmerge)
head(merged)
rownames(merged)<-merged[,1]
merged.all<-merged[,2:21]
### plot of the number of changes per site
hist(20-rowSums(is.na(merged.all)),breaks=0:20,main="frequency of sites changing per diet",xlab="# of diets causing change")
merged<-merged[rowSums(is.na(merged))<20,]

head(mergedNA)
mergedNA<-merged[,2:21]
merged<-merged[,2:21]
rownames(mergedNA)
merged[is.na(merged)]<-0
head(merged)
merged.dist<-dist(merged)
merged.clust<-hclust(merged.dist)
merged.naomit<-mergedNA[rowSums(is.na(mergedNA))<20,]
names(merged.naomit)
heatmap.2(as.matrix(mergedNA),scale="none", trace="none", Rowv=as.dendrogram(merged.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette)
nrow(mergedNA)
### plot all three 
dev.off()
par(mfcol=c(3,1))
as.matrix(kac.onlysig,colnames=paste(1:10))[1:5,]
colnames(kac.onlysig)<-1:10
colnames(suk.onlysig)<-1:10
colnames(merged.naomit)<-c(1:10,1:10)
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 100)

dev.off()
par(cex=2)
heatmap.2(as.matrix(merged.naomit),scale="none", trace="none", colsep=c(5,10,15), labRow = F,key.title = "",tracecol="black",Rowv=as.dendrogram(merged.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette)
nrow(merged.naomit)
top10pct<-merged.naomit[rowSums(is.na(merged.naomit))<=11,]
top5pct<-merged.naomit[rowSums(is.na(merged.naomit))<=9,]
nrow(top5pct)
merged.naomit-> merged.diff
ac.count<-10-rowSums(is.na(merged.diff[,1:10]))
suc.count<-10-rowSums(is.na(merged.diff[,11:20]))

difference<-ac.count-suc.count
difference
diff.hist<-hist(difference,breaks=20)
diff.hist
min(difference)
max(difference)
merged.diff[,1:10]
most.succ<-merged.naomit[difference<=-4,]
most.ac<-merged.naomit[difference>=8,]
nrow(most.succ)
nrow(most.ac)
heatmap.2(as.matrix(most.succ),na.color="grey",scale="none",RowV=FALSE,ColV=FALSE,dendrogram="none",col=my_palette)
heatmap.2(as.matrix(most.succ),na.color="grey",scale="none",RowV=F,ColV=F,col=my_palette)
heatmap.2(as.matrix(most.succ),scale="none", trace="none", colsep=c(5,10,15), labRow = F,key.title = "",tracecol="black",Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette)
heatmap.2(as.matrix(most.ac),scale="none", trace="none", colsep=c(5,10,15), labRow = F,key.title = "",tracecol="black",Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette)


write.table(most.ac,file="74most_acetyl_diff.txt", col.names = T,sep="\t")
write.table(most.succ,file="69most_succinyl_diff.txt", col.names = T,sep="\t")

top10pct
nrow(top10pct)
top10pct.4clust<-top10pct
top10pct.4clust[is.na(top10pct)]<-0
top10dist<-dist(top10pct.4clust)
top10clust<-hclust(top10dist)

heatmap.2(as.matrix(top10pct),scale="none", trace="none", labRow = F,colsep=c(5,10,15),tracecol="black",key.title = "",Rowv=as.dendrogram(top10clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette)

#### top 5pct clustered heatmap
nrow(top10pct)
top5pct.4clust<-top5pct
top5pct.4clust[is.na(top5pct)]<-0
top5dist<-dist(top5pct.4clust)
top5clust<-hclust(top5dist)

heatmap.2(as.matrix(top5pct),scale="none", trace="none", labRow = F,colsep=c(5,10,15),tracecol="black",key.title = "",Rowv=as.dendrogram(top5clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette)


nrow(top10pct)
?par
par(lwd=4,cex.lab=2,cex.axis=2,mfcol=c(1,3),mai=c(1,1,1,0.2))
library(gage)
??gage
histobj<-hist(20-rowSums(is.na(merged.naomit)),lwd=3,breaks=0:20,main="frequency of sites changing per diet",xlab="number of diets")
hist(10-rowSums(is.na(kac.naomit)),breaks=seq(from=0,to=10,by=1),lwd=3,ylab="number of sites",xlab="number of diets",main="diets causing acetyl site changes")
shist<-hist(10-rowSums(is.na(suk.naomit)),lwd=3,breaks=seq(from=0,to=10,by=1),ylab="number of sites",xlab="number of diets",main="diets causing succinyl site changes")
histobj<-hist(20-rowSums(is.na(merged.naomit)),lwd=3,breaks=0:20,main="# of diets causing ",xlab="number of diets",ylab="number of sites")
dev.off()
histobj
sum(histobj$counts)
write.table(merged.naomit,file="merged_sites_changes_naomit2.txt",sep="\t")
write.table(top10pct,file="merged_sites_changes_top10pct2.txt",sep="\t")

nrow(merged.naomit)
##### correlation matrix between acetyl and succinyl matricies


full.ack<-merged[,1:10]
head(full.ack)
full.suk<-merged[,11:20]
cor(full.ack,full.suk)
cor(full.ack)
### cor btw suk/ack
heatmap.2(as.matrix(cor(full.ack,full.suk)),scale="none", trace="none", Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)
### cor within ack
heatmap.2(as.matrix(cor(full.ack)),scale="none", trace="none", Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)
### cor within ack
heatmap.2(as.matrix(cor(full.suk)),scale="none", trace="none", Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)
### cor within ack
heatmap.2(as.matrix(cor(merged)),scale="none", trace="none", Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=bluered)



kac.rowvals<-as.character(kac.castFC[,1])
suk.rowvals<-as.character(suk.castFC[,1])

kac.protchar<-substr(kac.rowvals,start=1,stop=6)
suk.protchar<-substr(suk.rowvals,start=1,stop=6)

#both<-cbind(suk.castFC)



#suk.onlysig[suk.significant==FALSE]<-0
merge(kac.FC.order,suk.FC.order,by=0,all=TRUE)->mergeFC
head(kac.FC.order)
head(suk.FC.order)

head(mergeFC)

merge(kac.FDR.order,suk.FDR.order,by=0,all=TRUE)->mergeFDR

mergeFDR[1:10,1:3]
mergeFC[1:10,1:3]


mergedFC.uniprot<-substr(mergeFC[,1],start=1,stop=6)
mergedFDR.uniprot<-substr(mergeFDR[,1],start=1,stop=6)

cbind(mergedFC.uniprot,mergeFC)->mergeFC
head(mergeFC)
cbind(mergedFDR.uniprot,mergeFDR)->mergeFDR



fa.FC.sitelines<-mergeFC[which(is.na(match(mergeFC[,1],unlist(fa_ox)))==FALSE),3:22]
fa.FDR.sitelines<-mergeFDR[which(is.na(match(mergeFC[,1],unlist(fa_ox)))==FALSE),3:22]

unlist(fa_ox)

fa.fdr.TF<-fa.FDR.sitelines<=0.01
fa.FC.sitelines[fa.fdr.TF==FALSE]<-NA

rownames(fa.FC.sitelines)<-rownames(TAGsites)

rownames(fa.FC.sitelines)<-mergeFC[which(is.na(match(mergeFC[,1],unlist(fa_ox)))==FALSE),2]



heatmap.2(as.matrix(fa.FC.sitelines),scale="none", trace="none", Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=bluered,margins = c(15,10))


heatmap.2(as.matrix(kac.castFC[,order]),scale="row", trace="none", Rowv=TRUE, Colv=FALSE, dendrogram="none")

head(suk.onlysig)




protchar<-substr(rowvals,start=1,stop=6)
mapdia.castFC1<-cbind(protchar,mapdia.castFC)
mapdia.castFDR1<-cbind(protchar,mapdia.castFDR)


head(mapdia.castFC1)
head(mapdia.castFDR1)

mapdia.castFC1[which(mapdia.castFC1[,1]=="Q91Y97"),]
mapdia.castFDR1[which(mapdia.castFC1[,1]=="Q91Y97"),]


which(mapdia.castFC1)

control10w<-levels(mapdia[,"Label2"])[c(2,4,6,8,10)]

cont10w.rows<-mapdia[which(mapdia[,"Label2"]==control10w),]

which(as.character(cont10w.rows[,2])==as.character(targetIDs))
proteinbyrow<-as.character(cont10w.rows[,2])



which(proteinbyrow==targetIDs[1])
rowvec<-c()
for(x in targetIDs){
  rowvec<-c(rowvec,which(protchar==x))
}
onlyTAG<-mapdia.castFC1[rowvec,]
onlyTAGfdr<-mapdia.castFDR[rowvec,]

pt01fdr<-onlyTAGfdr[,2:17]<=0.01
pt01fdr<-onlyTAGfdr[,2:17]>=0.01

pt01fdr[,1]
onlyTAG[,1]
onlyTAG[,2]

for(i in 1:16){
  onlyTAG[,i+2][pt01fdr[,i]==TRUE]<-NA
}

pt01fdr
onlyTAG[is.na(onlyTAG)]
a<-onlyTAG[,order]
a[!is.na(a)]<-"*"

onlyTAGfdr
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=100)
graphics.off() 
tiff(filename = "test13_TAGbiosynth.tiff", width = 2000, antialias ="default", height =10000, units = "px", pointsize = 32)
onlyTAG[,2]
rownames(onlyTAG)<-onlyTAG[,2]

3+c(2,4,6,8,10)
order<-c(12,6,10,14,8)

heatmap.2(as.matrix(onlyTAG[,order]),symm=TRUE, trace="none",margins = c(20, 10), Rowv=FALSE, Colv=FALSE, na.rm=FALSE, na.color="grey",dendrogram="none",cexCol=1.0, lwid=c(4,4), lhei=c(.1,4), cexRow=.8,col=my_palette)
graphics.off()
heatmap.2(as.matrix(onlyTAG[,3:18]),scale="none", na.rm=FALSE, na.color="grey",
          trace="none",margins = c(15, 10), Rowv=FALSE, Colv=FALSE, dendrogram="none",cexCol=1.8, lwid=c(4,4), lhei=c(.1,4), cexRow=.8)



ncol(onlyTAGfdr)
ncol(onlyTAG)
library(gplots)
pdf()

?heatmap.2

rowsplit<-strsplit(as.character(cont10w.rows[,1]),split="_")
which()
unlist(rowsplit)


columns<-colnames(mapdia)
write.csv(suk.castFC, "succinyl.FC.dcast.csv", row.names=FALSE)
write.csv(kac.castFC, "acetyl.FC.dcast.csv", row.names=FALSE)


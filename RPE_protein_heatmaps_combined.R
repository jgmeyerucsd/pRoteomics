
### new quant
getwd()

setwd("~/3d_eye/manuscript/")

list.files()

#rpe<-read.delim(file="intracellular/RPEjune2016_Candidates.tsv", sep = "\t")
rpe<-read.delim(file= "intracellular/unfiltered_intracellular_RPE_Candidates.txt", sep = "\t")
pr<-read.delim(file="PR_Candidates.tsv")
sec<-read.delim(file="Secretome_Candidates.tsv")


#head(rpe)
#colnames(rpe)
#rpe[1,""]

head(prot)
library(reshape2)
##############################################################
### rpe set
colnames(rpe)
rpe.castFC <- dcast(data = rpe, formula = Group ~ Comparison..group1.group2.,value.var=c("AVG.Log2.Ratio"))
rpe.castFDR <- dcast(data = rpe, formula = Group ~ Comparison..group1.group2.,value.var=c("Qvalue"))
head(rpe.castFC)
rpe.castFC[1,]
row.names(rpe.castFC)
### reorganize to have uniprot as row name
rpe.castFC.1<-data.frame(rpe.castFC[,2:ncol(rpe.castFC)])
head(rpe.castFC.1)
row.names(rpe.castFC.1)
row.names(rpe.castFC.1)<-substr(rpe.castFC[,1],start=1,stop=6)
rpe.castFDR.1<-data.frame(rpe.castFDR[,2:ncol(rpe.castFC)])
row.names(rpe.castFDR.1)<-substr(rpe.castFDR[,1],start=1,stop=6)
colnames(rpe.castFDR.1)<-unlist(lapply(FUN=paste, colnames(rpe.castFDR.1), "intra", sep=".", collapse = ""))
colnames(rpe.castFC.1)<-unlist(lapply(FUN=paste, colnames(rpe.castFC.1), "intra", sep=".", collapse = ""))

head(rpe.castFC.1)
rpe.castFC[1,]
#rpe.castFC.1<-data.frame(rpe.castFC[,7])*-1  #### flip the signs because the ratio is backwards
#rownames(rpe.castFC.1)<-substr(rpe.castFC[,1],start=1,stop=6)
#rpe.castFDR.1<-data.frame(rpe.castFDR[,7])
#rownames(rpe.castFDR.1)<-substr(rpe.castFDR[,1],start=1,stop=6)


##############################################################
### pr set
pr.castFC <- dcast(data = pr, formula = Group ~ Comparison..group1.group2.,value.var=c("AVG.Log2.Ratio"))
pr.castFDR <- dcast(data = pr, formula = Group ~ Comparison..group1.group2.,value.var=c("Qvalue"))
head(pr.castFC)
typeof(pr.castFC)

### reorganize to have uniprot as row name
pr.castFC.1<-data.frame(pr.castFC[,2])
row.names(pr.castFC.1)<-substr(pr.castFC[,1],start=1,stop=6)
pr.castFDR.1<-data.frame(pr.castFDR[,2])
row.names(pr.castFDR.1)<-substr(pr.castFDR[,1],start=1,stop=6)

head(pr.castFC.1)
head(pr.castFDR.1)
typeof(pr.castFC.1)


###########################################################
### secr set
sec.castFC <- dcast(data = sec, formula = Group ~ Comparison..group1.group2.,value.var=c("AVG.Log2.Ratio"))
sec.castFDR <- dcast(data = sec, formula = Group ~ Comparison..group1.group2.,value.var=c("Qvalue"))
head(sec.castFC)
ncol(sec.castFC)
typeof(sec.castFC)
sec.castFC[1,]


sec.castFC.1<-data.frame(sec.castFC[,c(2,4)])*-1  #### flip the signs because the ratio is backwards
sec.castFC.1<-cbind(sec.castFC[,3],sec.castFC.1)
head(sec.castFC.1)
sec.castFC.1<-sec.castFC.1[,c(1,3,2)]
rownames(sec.castFC.1)<-substr(sec.castFC[,1],start=1,stop=6)
colnames(sec.castFC.1)<-c("1wPQ.sec","3wPQ.sec","3wPQ/1wPQ.sec")
sec.castFDR.1<-data.frame(sec.castFDR[,c(3,4,2)])
rownames(sec.castFDR.1)<-substr(sec.castFDR[,1],start=1,stop=6)
colnames(sec.castFDR.1)<-c("1wPQ.sec","3wPQ.sec","3wPQ/1wPQ.sec")


#############################3
### commands to check the lists
head(sec.castFC.1)

head(rpe.castFC.1)
head(rpe.castFDR.1)
head(pr.castFC.1)



###### how many rows each?

nrow(rpe.castFC.1)  ### 3257
nrow(pr.castFC.1)  ### 3429
nrow(sec.castFC.1) # 1319

### how many overlap?
?union
setdiff(row.names(rpe.castFC.1),row.names(pr.castFC.1))
intersect(row.names(rpe.castFC.1),row.names(pr.castFC.1))

length(intersect(row.names(rpe.castFC.1),row.names(pr.castFC.1)))


3429-2530
3257-2530
length(intersect(row.names(rpe.castFC.1),row.names(sec.castFC.1)))
1319-1054
length(intersect(row.names(pr.castFC.1),row.names(sec.castFC.1)))
1319-936
Reduce(intersect, list(row.names(rpe.castFC.1),
                                 row.names(pr.castFC.1),
                                          row.names(sec.castFC.1))
)

length(unique(c(row.names(rpe.castFC.1),
       row.names(pr.castFC.1),
       row.names(sec.castFC.1))))

1054-906
2530-906
579+1624+869+148+906+30+235

###### merge the dataframes ----   fold changes
#   ?merge
merged1<-merge(rpe.castFC.1,sec.castFC.1,all = TRUE, by=0)
head(merged1)
merged2<-cbind(merged1[,2:ncol(merged1)])
row.names(merged2)<-merged1[,1]

head(merged2)
merged3<-merge(merged2,pr.castFC.1,all = TRUE, by=0)
head(merged3)
merged.FC<-cbind(merged3[,2:ncol(merged3)])
rownames(merged.FC)<-merged3[,1]
head(merged.FC)
colnames(merged.FC)<-c("PQ1.i","PQ3.i","PQ3/PQ1.i","PQ1.s","PQ3.s","PQ3/PQ1.s","PR")


###### merge the FDRs
merged1<-merge(rpe.castFDR.1,sec.castFDR.1,all = TRUE, by=0)
head(merged1)
head(pr.castFDR.1)
merged2<-cbind(merged1[,2:ncol(merged1)])
row.names(merged2)<-merged1[,1]

head(merged2)
merged3<-merge(merged2,pr.castFDR.1,all = TRUE, by=0)

merged.FDR<-cbind(merged3[,2:ncol(merged3)])
rownames(merged.FDR)<-merged3[,1]
colnames(merged.FDR)<-c("PQ1.i","PQ3.i","PQ3/PQ1.i","PQ1.s","PQ3.s","PQ3/PQ1.s","PR")

head(merged.FC)
head(merged.FDR)
head(merged)


## ##########################       #################This is out of place, requires genemap.all
#geneids.unfiltered <- genes[match(row.names(merged.FC), genemap.all[,"UNIPROTKB"])]
#row.names(merged.FC)[12]
#geneids.unfiltered[12]
#row.names(merged.FC)<-geneids.unfiltered
#write.table(file="merged3sets_FCs_unfiltered.txt",sep="\t",merged.FC)

#### replace non-significant values in the matrix with NA
head(rpe.castFC.1)
head(rpe.castFDR.1)
prot.onlysig<-merged.FC
nrow(merged.FC)
prot.onlysig[is.na(merged.FC)]<-0
prot.onlysig[merged.FDR>=0.01]<-0

##prot.onlysig[abs(merged.FC)<log(1.2,base=2)]<-0   ## filter for a fold change
#head(prot.onlysig)
#log(1.2,base=2)

### put names back
#row.names(kac.FC.order)<-as.character(kac.castFC[,1])
nrow(prot.onlysig)
ncol(prot.onlysig)
#rownames(prot.filt)
prot.filt<-prot.onlysig[rowSums(prot.onlysig==0)<7,]
nrow(prot.filt)
prot.filt.dist<-dist(prot.filt)
prot.filt.clust<-hclust(prot.filt.dist)
prot.filt[prot.filt==0]<-NA
nrow(prot.filt)
head(prot.filt,10)



write.table(file="20180409_newSecr_mergedFC_q01_noFC.txt",sep="\t",prot.filt)
#### start here #####

read.delim(file="Z:/R24/resource/proteinlvl/254candidates.onlysig.txt",header = T)->prot.filt


#### gene mapping, only need to run once per set of input
library(RSQLite)
library(UniProt.ws)

up <- UniProt.ws(taxId=9606)
taxId(up) <- 9606
split_sp <- function(x) unlist(strsplit(x, " "))[1]

uniprot<-row.names(merged.FC)
genemap.all<-select(up, uniprot, "GENES")

filt.uniprot<-row.names(prot.filt)
#genemap.filt <- select(up,filt.uniprot, "GENES")
genes <- unlist(lapply(genemap.all[,2], split_sp))
geneid <- genes[match(filt.uniprot, genemap.all[,"UNIPROTKB"])]

filt.uniprot[1]
row.names()

#row.names(prot.filt)[1]
#geneid[1]
row.names(prot.filt)<-geneid


### plot heatmap
library(RColorBrewer)
?colorRampPalette
my_palette <- colorRampPalette(c("red", "white", "forestgreen"))(n=17)
library(gplots)
?heatmap.2
install.packages("viridis")

breaks=seq(-2, 2, by=0.05) #41 values
#now add outliers
breaks=append(breaks, 10)
breaks=append(breaks, -10, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="green",mid="white",high="red")
#now plot heatmap - I'm using heatmap.2, but am assuming this will work with Obi's heatmap.3


hm<-heatmap.2(as.matrix(prot.filt),main = "1749 proteins",scale="none", 
              key.par=list(),trace="none",  keysize = 2,
              key.title = "",tracecol="black",Rowv=as.dendrogram(prot.filt.clust), 
              Colv=F, na.col="grey",dendrogram="none",col=mycol,
              colsep=c(0,3,6),
              sepcolor = 1,
              margins = c(5,15),
              breaks=breaks,
              cexRow = 0.2)

dev.off()
max(na.omit(unlist(prot.filt)))


###################################33
####################33
#separate subset of the intracellular #### 1212 proteins at this cutoff

filt.intra<-prot.filt[,1:3]
filt.extra<-prot.filt[,4:6]
filt.pr<-prot.filt[,7]


filt.intra[is.na(filt.intra)]<-0
filt.intra<-filt.intra[rowSums(filt.intra==0)<3,]
nrow(filt.intra)
filt.intra.dist<-dist(filt.intra)
filt.intra.clust<-hclust(filt.intra.dist)
filt.intra[filt.intra==0]<-NA
nrow(filt.intra)
head(prot.filt,10)


###########################3
#my_palette <- colorRampPalette(c("red", "white", "blue"))(n=17)
hm<-heatmap.2(as.matrix(filt.intra),
              main="311 proteins",
              scale="none", 
              key.par=list(),
              trace="none",
              keysize = 2,
              key.title = "",
              tracecol="black",
              Rowv=as.dendrogram(filt.intra.clust), 
              Colv=F, 
              na.col="white",
              dendrogram="none",
              col=mycol,
              colsep=c(0,1,2,3),
              breaks=breaks,
              
              sepcolor = 1,margins = c(5,15))



write.table(file="intra_FC1p5_q05_2015data.txt",sep="\t",filt.intra[hm$rowInd,])


###################################33
####################33
#separate subsets of the secreted ###712 changes secreted


#filt.extra<-prot.filt[,5:7]

filt.extra[is.na(filt.extra)]<-0
filt.extra<-filt.extra[rowSums(filt.extra==0)<3,]
nrow(filt.extra)
filt.extra.dist<-dist(filt.extra)
filt.extra.clust<-hclust(filt.extra.dist)
filt.extra[filt.extra==0]<-NA
nrow(filt.extra)
head(prot.filt,10)

###########################3
hm<-heatmap.2(as.matrix(filt.extra),main="311 proteins",scale="none", 
              key.par=list(),trace="none",  keysize = 2,
              key.title = "",tracecol="black",Rowv=as.dendrogram(filt.extra.clust), 
              Colv=F, na.col="white",dendrogram="none",col=my_palette,colsep=c(0,1,2,3),sepcolor = 1,margins = c(5,15))



write.table(file="intra_FC1p5_q05_2015data.txt",sep="\t",filt.intra[hm$rowInd,])

###################################33
####################33
#separate subsets of the PR


filt.pr<-data.frame(prot.filt[,7])
row.names(filt.pr)<-row.names(prot.filt)
filt.pr[is.na(filt.pr)]<-0
filt.pr<-filt.pr[(filt.pr==0)==FALSE]
nrow(filt.pr)
filt.pr.dist<-dist(filt.pr)
filt.pr.clust<-hclust(filt.pr.dist)
filt.pr[filt.pr==0]<-NA
nrow(filt.extra)
head(prot.filt,10)

###########################3
hm<-heatmap.2(as.matrix(filt.pr),main="311 proteins",scale="none", 
              key.par=list(),trace="none",  keysize = 2,
              key.title = "",tracecol="black",Rowv=as.dendrogram(filt.extra.clust), 
              Colv=F, na.col="white",dendrogram="none",col=my_palette,colsep=c(0,1,2,3),sepcolor = 1,margins = c(5,15))



write.table(file="intra_FC1p5_q05_2015data.txt",sep="\t",filt.intra[hm$rowInd,])




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
up <- UniProt.ws(taxId=9606)
taxId(up) <- 9606
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

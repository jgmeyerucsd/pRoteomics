setwd("C:/users/jmeyer/documents/R24/dietstudy/")



read.delim(file="ack_analysis_output.txt")->k
read.delim(file="suck_analysis_output.txt")->s




### need sites and separately the significant ones



library(reshape2)
### cast acetyl
k.castFC <- dcast(data = k, formula = Protein ~ Label2,value.var=c("log2FC"))
k.castFDR <- dcast(data = k, formula = Protein ~ Label2,value.var=c("FDR"))
s.castFC <- dcast(data = s, formula = Protein ~ Label2,value.var=c("log2FC"))
s.castFDR <- dcast(data = s, formula = Protein ~ Label2,value.var=c("FDR"))
head(k.castFC)
head(s.castFC)

dev.off()
plot(princomp(as.matrix(k.fc.order)))
colnames(k.fc.order)<-unlist(strsplit(colnames(k.fc.order),split="/"))[seq(from=1,to=19,by=2)]
plot.new()
plot(xlim=c(-1,1),ylim=c(-1,1))
head(m.fc)
m.fc[is.na(m.fc)]<-0
pca3d(t(as.matrix(k.fc.order)),show.labels = TRUE,fancy =F,ellipse.ci = 0.95,components = 1:3)
par()
pca2d(t(as.matrix(k.fc.order)),show.labels = TRUE,fancy =F,ellipse.ci = 0.95,radius=2)
pca2d(t(as.matrix(k.fc.order)),show.labels = TRUE,fancy =F,ellipse.ci = 0.95,radius=2)
pca2d(t(as.matrix(k.fc.order)),radius=2)


colnames(s.fc.order)<-unlist(strsplit(colnames(s.fc.order),split="/"))[seq(from=1,to=19,by=2)]

pca2d(t(as.matrix(s.fc.order)),radius=2)
dev.off()
head(s.fc.order)
sall.hm<-heatmap.2(as.matrix(s.fc.order), scale="none",trace=c("none"),
          margins = c(10,10),
          tracecol="black",
          col=bluered(5),
          dendrogram = "column",
          key.title = "",
          labRow = F)
aall.hm<-heatmap.2(as.matrix(k.fc.order), scale="none",trace=c("none"),
          margins = c(10,10),
          tracecol="black",
          col=bluered(5),
          dendrogram = "column",
          key.title = "",
          labRow = F)

#### kegg enrichment
library(clusterProfiler)
search_kegg_organism('musculus', by='scientific_name') #mmu
?enrichKEGG
head(k.fc.order)
head(s.fc.order)
#need entrez geneid
aall<-k.fc.order[aall.hm$rowInd,]
aall<-k.fc.order[aall.hm$rowInd,aall.hm$colInd]

head(aall)

a.uniprot<-substr(row.names(aall),start=1,stop=6)


#### additional pathways
library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]
?select
keytypes(up)
genemap <- select(up,unique(a.uniprot) , "ENTREZ_GENE")
genes <- unlist(lapply(genemap[,2], split_sp))
a.entrez <- genes[match(a.uniprot, genemap[,"UNIPROTKB"])]
a.entrez
write.table(aall, file="a.all.clust.txt",sep="\t")
#### check the subset
heatmap.2(as.matrix(aall[1:92,]), scale="none",trace=c("none"),
                   margins = c(10,10),
          Rowv = NA,
          Colv = NA,
                   tracecol="black",
                   col=bluered(5),
                   key.title = "",
                   labRow = F)

### first 92 are down cluster
a.entrez[1:92]
a.uniprot[1:92]
a.kk <- enrichKEGG(gene         = a.entrez[1:92],
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(a.kk)
mkk2 <- gseMKEGG(geneList = na.omit(a.uniprot[1:92]))
test.gl
sort

  typeof(geneList)
  names(geneList)
kk2 <- gseKEGG(geneList     = na.omit(a.uniprot[1:92]),
               organism     = 'mmu',
               keyType = "uniprot",
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
library(DOSE)
data(geneList)
de <- names(geneList)[geneList > 1]

heatmap.2(as.matrix(k.fc.order), scale="none",trace=c("none"), margins = c(10,10), col=bluered(11))
?heatmap.2
dev.off()
pca3d(t(as.matrix(m.fc)),show.labels = TRUE,fancy =F,ellipse.ci = 0.95,components = 1:3)

comp(k.castFC)
#### order the columns so they are the order 2w vs control, 10w vs control, 10w vs 2w
order2w<-c(14,15,16,17,18)
order10w<-c(3,5,7,9,11)
order2v10<-c(2,4,6,8,10,13)
#order<-c(order2w,order10w,order2v10) #all comparisons
order<-c(order2w,order10w) # only 2w or 10w
order<-order2v10 ###only 2v10


k.fc.order<-k.castFC[,order]
k.fdr.order<-k.castFDR[,order]
s.fc.order<-s.castFC[,order]
s.fdr.order<-s.castFDR[,order]
?pca

head(k.fc.order)
head(k.fdr.order)

head(s.fc.order)
head(s.fdr.order)

# make the row names match
row.names(k.fc.order)<-gsub(k.castFC[,1],pattern=".K\\+42.011",replacement = "")
row.names(k.fdr.order)<-gsub(k.castFC[,1],pattern=".K\\+42.011",replacement = "")

row.names(s.fc.order)<-gsub(s.castFC[,1],pattern=".K\\+100.016",replacement = "")
row.names(s.fdr.order)<-gsub(s.castFC[,1],pattern=".K\\+100.016",replacement = "")

#row.names(sl)<-gsub(sl[,1],pattern=".K\\+42.0105",replacement = "")

#### set nonsign. to NA, count rows with sig
a.onlysig<-k.fc.order
a.onlysig[k.fdr.order>=0.01]<-NA
a.onlysig[abs(a.onlysig)<1]<-NA
a.filtered<-a.onlysig[rowSums(is.na(a.onlysig))!=ncol(a.onlysig),]
nrow(a.onlysig)
nrow(a.filtered) ## 653 if using >2 fold, FDR<0.01
unique(substr(row.names(a.filtered),start=1,stop=6))  ### 251 proteins of those have at least one significant change in acylation site

#### same for succinyl
s.onlysig<-s.fc.order
s.onlysig[s.fdr.order>=0.01]<-NA
s.onlysig[abs(s.onlysig)<1]<-NA
s.filtered<-s.onlysig[rowSums(is.na(s.onlysig))!=ncol(s.onlysig),]
nrow(s.onlysig)
nrow(s.filtered) ## 653 if using >2 fold, FDR<0.01
unique(substr(row.names(s.filtered),start=1,stop=6))  ### 251 proteins of those have at least one significant change in acylation site




uniprot<-substr(row.names(m.filtered),start=1,stop=6)
m<-m.filtered
m[is.na(m)]<-0
m.dist<-dist(m,method = "manhattan")
m.dist<-dist(m,method = "maximum")

m.clust<-hclust(m.dist)
?dist
m[m==0]<-NA
library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white","red"))(n = 11)

m.hm<-heatmap.2(as.matrix(s.fc.order),scale="none", trace="none", 
                tracecol="black", colsep=c(6), 
                sepcolor = "black",
                na.col="grey",
                Rowv=T,
                dendrogram="none",col=my_palette, margins = c(10,10))

#### let them order themselves
m.hm<-heatmap.2(as.matrix(m[,16:20]),scale="none", trace="none", 
                col=my_palette, margins = c(15,15), density.info="none")

# including time dependent
m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", tracecol="black", colsep=c(5,10,16,21,26), Rowv=as.dendrogram(m.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette,margins = c(10,10))
## including only 2v10
m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", tracecol="black", colsep=c(6), Rowv=as.dendrogram(m.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette,margins = c(10,10))

m.hm$rowInd
t<-cbind(uniprot,m)
t[m.hm$rowInd,]
write.table(t[m.hm$rowInd,],sep="\t",file="PIQED_hm2w10w_sKfix.txt")
m<-read.delim(file="mergedsites_PIQED_stringent.txt")
head(m)

uniprot<-substr(m[,1],start=1,stop=6)

#### additional pathways
library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]
genemap <- select(up,unique(uniprot) , "GENES")


genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]



geneids<-toupper(geneid)

m[354,]
geneids[geneids=="PC"]
geneids[geneids=="SDHD"]



##### combine the pathways into one ---
list.files()
bOx<-read.delim("customTCA.txt",skip = 1,head=T,stringsAsFactors = F)

unlist(bOx)
tmp.index<-is.na(match(geneids,unlist(bOx)))
tmpset<-m[tmp.index==FALSE,]
rownames(tmpset)<-paste(geneids[tmp.index==FALSE],substr(rownames(tmpset),start=8,stop=nchar(rownames(tmpset))))
tca<-tmpset
tca<-tca[order(row.names(tca)),]


library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=11)
#?heatmap.2
m.hm<-heatmap.2(as.matrix(tca[,c(6:10,16:20)]),scale="none", trace="none", tracecol="black", colsep=c(5,10,15), Rowv=F, Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette,cexRow = 0.4)




fbox.aves<-group.aves(tmpset[,31:59], groups=NA)
nrow(box.aves[[2]])






?tiff
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
geneids[218:220]


kac.FC<-cbind(geneid,kac.castFC)


?select
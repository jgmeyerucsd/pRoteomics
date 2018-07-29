

setwd("C:/users/jmeyer/documents/R24/dietstudy/")


read.delim(file="site_level_aK.txt", row.names = 1)->k.sl
read.delim(file="site_level_sK.txt", row.names = 1)->s.sl
read.delim(file="protlvl/proteinlevels.txt", row.names = 1)->p


head(k.sl)
head(p)

k.sl.c<-k.sl[,1:(ncol(k.sl)-2)]
sort(colnames(p))
p<-p[,order(colnames(p))]
colnames(p)
head(p)

p.cor<-cor(p)
dev.off()
par(mfcol=c(1,3))
p.cor[lower.tri(p.cor,diag=T)]<-NA
min(na.omit(as.numeric(p.cor))) ### protein lvl min=0.757
mean(na.omit(as.numeric(p.cor))) ### protein lvl mean=0.96
corcol<-colorRampPalette(c("white","red"))

heatmap.2(p.cor,scale="none", trace="none",
          breaks=seq(from=0.69,to=1,by=0.01),
          key.title = "",
          tracecol = "black",
          Rowv = NA, 
          Colv = F, dendrogram = "none", col = corcol,sepcol="black",
          colsep=seq(from=5, to=55,by=5))



a.cor.mat<-cor(k.sl.c)
a.cor.mat[lower.tri(a.cor.mat,diag=T)]<-NA
nrow(cor.mat)
min(na.omit(as.numeric(a.cor.mat))) ### protein lvl min=0.69
mean(na.omit(as.numeric(a.cor.mat))) ### protein lvl mean=0.91
mode(na.omit(as.numeric(a.cor.mat))) ### acetyl lvl mean=0.91


heatmap.2(a.cor.mat,scale="none", trace="none",
          breaks=seq(from=0.69,to=1,by=0.01),
          key.title = "",
          tracecol = "black",
          Rowv = NA, 
          Colv = F, dendrogram = "none", col = corcol,sepcol="black",
            colsep=seq(from=5, to=55,by=5))
#### heatmap of the 
as.matrix(log(k.sl.c))->log.kac
is.nan
length(log.kac[is.na(log.kac)])
log.kac[log.kac==-Inf]<-NA
rowSums(log.kac,is)
?heatmap.2

k.d<-dist(t(log.kac))
k.c<-hclust(k.d)

heatmap.2(log.kac,scale="none", trace="none",
          key.title = "", 
          tracecol = "black",
          dendrogram = "none", 
          sepcol="black",
          distfun = dist)
#succinyl cor hm
head(k.sl)

head(s.sl)
ncol(s.sl)-2
s.sl.c<-s.sl[,1:(ncol(s.sl)-2)]
s.sl[1,1]
head(s.sl.c[,1])
corcol<-colorRampPalette(c("white","red"))
s.cor.mat<-cor(s.sl.c)
s.cor.mat[lower.tri(s.cor.mat,diag=T)]<-NA
nrow(cor.mat)
min(na.omit(as.numeric(s.cor.mat))) ### protein lvl min=0.728
mean(na.omit(as.numeric(s.cor.mat))) ### succinyl mean lvl min=0.88


?cor
?heatmap.2
heatmap.2(s.cor.mat,scale="none", trace="none",Rowv = NA, 
          Colv = F, 
          breaks=seq(from=0.69,to=1,by=0.01),
          key.title="",
          tracecol = "black",
          dendrogram = "none", 
          col = corcol,
          sepcol="black",
          colsep=sepn)
sepn<-c(3,7,11,15,20,25,30,34,38,43,48,52)

### need sites and separately the significant ones


##### Protein level correlation map






















library(reshape2)
### cast acetyl
k.castFC <- dcast(data = k, formula = Protein ~ Label2,value.var=c("log2FC"))
k.castFDR <- dcast(data = k, formula = Protein ~ Label2,value.var=c("FDR"))
s.castFC <- dcast(data = s, formula = Protein ~ Label2,value.var=c("log2FC"))
s.castFDR <- dcast(data = s, formula = Protein ~ Label2,value.var=c("FDR"))
head(k.castFC)
head(s.castFC)

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



head(k.castFC)[order2w]
head(k.castFC)[order10w]

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

m.fc<-merge(k.fc.order,s.fc.order,all=T,by=0)
m.fdr<-merge(k.fdr.order,s.fdr.order,all=T,by=0)

head(m.fc)
head(m.fdr)

row.names(m.fc)<-m.fc[,1]
row.names(m.fdr)<-m.fdr[,1]
m.fc<-m.fc[,c(2:ncol(m.fc))]
m.fdr<-m.fdr[,c(2:ncol(m.fdr))]


nrow(m.fc)
unique(substr(row.names(m.fc),start=1,stop=6))  ### 402 unique proteins with at least one mod quantified



m.onlysig<-m.fc
m.onlysig[m.fdr>=0.01]<-NA
m.onlysig[abs(m.onlysig)<1]<-NA


m.onlysig[is.na(m.onlysig)]<-0   #### from NA to zero
m.onlysig[m.onlysig==0]<-NA  #### from zero to NA

#colnames(m.onlysig)<-c("Fruc","Gluc","HFD","HFD+Fruc","HFD+Gluc","HFD+gluc/HFD+Fruc")

### remove the rows that don't have all significant 
m.filtered<-m.onlysig[rowSums(is.na(m.onlysig))!=ncol(m.onlysig),]


nrow(m.onlysig)
nrow(m.filtered) ## 653 if using >2 fold, FDR<0.01
unique(substr(row.names(m.filtered),start=1,stop=6))  ### 251 proteins of those have at least one significant change in acylation site

uniprot<-substr(row.names(m.filtered),start=1,stop=6)
m<-m.filtered
m[is.na(m)]<-0
m.dist<-dist(m)
m.clust<-hclust(m.dist)


m[m==0]<-NA

library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=100)
my_palette <- colorRampPalette(c("blue", "white","red"))(n=100)
my_palette <- colorRampPalette(c("red", "orange","yellow"),bias=10)(n=11)
?RColorBrewer
colors = c(seq(-4,-1.01,length=1),seq(-1,1,length=10),seq(1.01,4,length=1))

my_palette <- colorRampPalette(c("blue", "white","red"))(n = 11)

?colorRampPalette
#names(m)
#?heatmap.2
## excluding time dependent
m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", 
                tracecol="black", colsep=c(5,10,15), 
                sepcolor = "black",
                Rowv=as.dendrogram(m.clust), Colv=FALSE, 
                na.col="white", 
                dendrogram="none",col=my_palette)
# including time dependent
m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", tracecol="black", colsep=c(5,10,16,21,26), Rowv=as.dendrogram(m.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette,margins = c(10,10))
## including only 2v10
m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", tracecol="black", colsep=c(6), Rowv=as.dendrogram(m.clust), Colv=FALSE, na.col="grey",dendrogram="none",col=my_palette,margins = c(10,10))

m.hm$rowInd
t<-cbind(uniprot,m)
t[m.hm$rowInd,]
write.table(t[m.hm$rowInd,],sep="\t",file="PIQED_hm2v10_anyFC.txt")
m<-read.delim(file="mergedsites_PIQED_stringent.txt")
head(m)

uniprot<-substr(m[,1],start=1,stop=6)
m<-m[,2:ncol(m)]
head(m)




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

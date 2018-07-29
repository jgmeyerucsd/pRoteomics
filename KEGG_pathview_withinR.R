#### kegg enrichment
library(clusterProfiler)
search_kegg_organism('musculus', by='scientific_name') #mmu
search_kegg_organism('Escherichia', by='scientific_name') #mmu

?enrichKEGG
head(k.fc.order)
head(s.fc.order)
#need entrez geneid
aall<-k.fc.order[aall.hm$rowInd,]
aall<-k.fc.order[aall.hm$rowInd,aall.hm$colInd]

head(aall)

a.uniprot<-substr(row.names(aall),start=1,stop=6)

#### get list of sites regulated by given name
split_gene=function(x)unlist(strsplit(x, "_"))[1]
split_num=function(x)unlist(strsplit(x, "_"))[1]
bg.list<-substr(row.names(m),start=1,stop=6)

get.regulated.by.column=function(input.column="YfiQ",input.df=cbind(seq.windows,m)){
  split_gene=function(x)unlist(strsplit(x, "_"))[1]
  split_num=function(x)unlist(strsplit(x, "_"))[1]
  motifs<-input.df[,1]
  tmp.df<-input.df[,grep(input.column,colnames(input.df))]
  head(tmp.df)
  tmp.df2<-motifs[rowSums(is.na(tmp.df))<2]
  #head(tmp.df2)

  #uniprot.vec<-substr(row.names(tmp.df2),start=1,stop=6)
  #return(list(uniprot.vec,tmp.df2))
  return(tmp.df2)
}

head(m)
### yfiq

get.regulated.by.column(input.column="YfiQ",input.df=cbind(seq.windows,m))->yfiq
head(yjab)
head()
na.omit(yfiQ)
write.table(yfiq[[1]], file="yfiq_genes.txt",sep="\t",quote=F, row.names = F)
write.table(na.omit(yfiQ), file="yfiq_motifs.txt",sep="\t",quote=F, row.names = F)

get.regulated.by.column(input.column="YjaB",input.df=cbind(seq.windows,m))->yjab
head(yjab)
write.table(na.omit(yjab), file="yjab_motifs.txt",sep="\t",quote=F, row.names = F)

get.regulated.by.column(input.column="YiaC",input.df=cbind(seq.windows,m))->yiac
yiac
write.table(na.omit(yiac), file="yiac_motifs.txt",sep="\t",quote=F, row.names = F)

?enrichKEGG
yfiq.kegg <- enrichKEGG(gene         = yfiq[[1]],
                        organism     = 'eco',
                        keyType="uniprot")
head(yfiq.kegg)
-log(yfiq.kegg$qvalue,base=10)->tmp.x

divide=function(x)as.numeric(unlist(strsplit(x,"/")))[1]/as.numeric(unlist(strsplit(x,"/")))[2]
tmp.bg<-lapply(yfiq.kegg$BgRatio, divide)
tmp.gr<-lapply(yfiq.kegg$GeneRatio, divide)
tmp.y<-unlist(tmp.gr)/unlist(tmp.bg)
par(cex=1)
plot(tmp.x,tmp.y,pch=20,xlim=c(0,30),main="YfiQ, KEGG pathways", xlab="-log10(FDR)",ylab="fold enrichment")
text(tmp.x,tmp.y, yfiq.kegg$Description,pos=4)

# yjab
get.regulated.by.column(input.column="YjaB",input.df=m)->yjab
head(yjab)
yjab[[1]]


yjab.kegg <- enrichKEGG(gene         = yjab[[1]],
                   organism     = 'eco',
                   keyType="uniprot")
head(yjab.kegg)
-log(yjab.kegg$qvalue,base=10)->tmp.x

as.numeric(yfiq.kegg$BgRatio)
x<-yfiq.kegg$BgRatio[1]
divide=function(x)as.numeric(unlist(strsplit(x,"/")))[1]/as.numeric(unlist(strsplit(x,"/")))[2]
tmp.bg<-lapply(yjab.kegg$BgRatio, divide)
tmp.gr<-lapply(yjab.kegg$GeneRatio, divide)
tmp.y<-unlist(tmp.gr)/unlist(tmp.bg)

plot(tmp.x,tmp.y,pch=20,xlim=c(0,30),main="YfiQ, KEGG pathways", xlab="-log10(FDR)",ylab="fold enrichment")
text(tmp.x,tmp.y, yjab.kegg$Description,pos=4)


# yiac
get.regulated.by.column(input.column="YiaC",input.df=m)->yiac
head(yjab)
yiac[[1]]
write.table(yiac[[1]], file="yiac_genes.txt",sep="\t",quote=F, row.names = F)
write.table(bg.list, file="bg_genes.txt",sep="\t",quote=F, row.names = F)


yiac.kegg <- enrichKEGG(gene         = yiac[[1]],
                        organism     = 'eco',
                        keyType="uniprot",pvalueCutoff = 0.2)

yiac.kegg
?enrichGO

source("https://bioconductor.org/biocLite.R")
biocLite("org.EcK12.eg.db")
a
yiac.go <- enrichGO(gene         = yiac[[1]],
                    'org.EcK12.eg.db',
                        keytype="uniprot")

-log(yiac.kegg$qvalue,base=10)->tmp.x

divide=function(x)as.numeric(unlist(strsplit(x,"/")))[1]/as.numeric(unlist(strsplit(x,"/")))[2]
tmp.bg<-lapply(yiac.kegg$BgRatio, divide)
tmp.gr<-lapply(yiac.kegg$GeneRatio, divide)
tmp.y<-unlist(tmp.gr)/unlist(tmp.bg)
plot(tmp.x,tmp.y,pch=20,xlim=c(0,30),main="YiaC, KEGG pathways", xlab="-log10(FDR)",ylab="fold enrichment")
text(tmp.x,tmp.y, yiac.kegg$Description,pos=4)

plot(tmp.x,tmp.y)

head(yiac.kegg)

fa.kk$Description
mkk2 <- gseMKEGG(geneList = na.omit(a.uniprot[1:92]))
test.gl
sort

typeof(geneList)
names(geneList)
?gseKEGG
?gseKEGG
yiac.gse <- gseKEGG(geneList     = yiac[[1]],
               organism     = 'eco',
               keyType = "uniprot",
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

library(DOSE)
deg = names(geneList)[abs(geneList) > 1]
?enrichDO
do = enrichDO(yiac[[1]])
dotplot(do, showCategory=20)



library(DOSE)
data(geneList)

de <- names(geneList)[geneList > 1]
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


########## bacteria dataset, gene back to uniprot

library(UniProt.ws)
split_sp <- function(x)unlist(strsplit(x, " "))[1]

up <- UniProt.ws(taxId=83333)
taxId(up) <- 83333
### genes
?select
keytypes(up)
genemap <- select(up,unique(uniprot), "GENES")
genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]
geneids <- toupper(geneid)

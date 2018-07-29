



library("Biobase")
library("annotate")
library("Category")
library("hgu95av2")
library("genefilter")

source("https://bioconductor.org/biocLite.R")
biocLite("GSEABase")
biocLite("UniProt.ws")



setwd("C:/users/jmeyer/documents/GitHub/pRoteomics/")
read.csv(file="acK_1280sites_only.csv")->kac

### need sites and separately the significant ones
sl<-read.delim(file="protein_level.txt")
kac<-read.delim(file="analysis_output.txt")

head(sl)
head(kac)




library(reshape2)
### cast acetyl
kac.castFC <- dcast(data = kac, formula = Protein ~ Label2,value.var=c("log2FC"))
kac.castFDR <- dcast(data = kac, formula = Protein ~ Label2,value.var=c("FDR"))
head(kac.castFC)
order<-c(2,3,4,5,6,7)
kac.fc.order<-kac.castFC[,order]
kac.fdr.order<-kac.castFDR[,order]


# make the row names match
row.names(kac.fc.order)<-gsub(kac.castFC[,1],pattern=".K\\+42.0105",replacement = "")
row.names(kac.fdr.order)<-gsub(kac.castFC[,1],pattern=".K\\+42.0105",replacement = "")
row.names(sl)<-gsub(sl[,1],pattern=".K\\+42.0105",replacement = "")

row.names(kac.fc.order)
row.names(sl)
unique(substr(row.names(sl),start=1,stop=6))


kac.significant<-kac.fdr.order<=0.01

kac.onlysig<-kac.fc.order


head(kac.fc.order)

kac.onlysig[kac.fdr.order>=0.01]<-NA
head(kac.onlysig)
head(kac.fdr.order)
kac.onlysig[abs(kac.onlysig)<0.58]<-NA


kac.onlysig[is.na(kac.onlysig)]<-0   #### from NA to zero
kac.onlysig[kac.onlysig==0]<-NA  #### from zero to NA

head(kac.onlysig)
colnames(kac.onlysig)<-c("Fruc","Gluc","HFD","HFD+Fruc","HFD+Gluc","HFD+gluc/HFD+Fruc")

### remove the rows that don't have all significant 
kac.filtered<-kac.onlysig[rowSums(is.na(kac.onlysig))!=6,]

head(kac.onlysig)
head(kac.filtered)

nrow(kac.onlysig)
row
nrow(kac.filtered)
unique(substr(row.names(kac.filtered),start=1,stop=6))

uniprot<-substr(row.names(kac.filtered),start=1,stop=6)

slf<-sl[match(row.names(kac.filtered),row.names(sl)),]

head(slf)

head(kac.filtered)

ncol(slf)
slf<-slf[,2:60]
head(slf)
head(slf[,31:32])

slf.ave<-group.aves(slf[,31:59],groups=NA)


slf.aves<-slf.ave[[2]]
slf.aves



uniprot<-substr(row.names(slf),start=1,stop=6)

setwd("Z:/R24/finalsitelvl/")
getwd()
pathway.files<-list.files()

pathway.files


#### additional pathways

biocLite("RSQLite")
library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]
genemap <- select(up,unique(uniprot) , "GENES")


genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]



geneids<-toupper(geneid)


geneids[geneids=="EHHADH"]
geneids[geneids=="CAT"]



##### combine the pathways into one ---
list.files()
bOx<-read.delim("customTCA.txt",skip = 1,head=T,stringsAsFactors = F)
bOx<-read.delim("C:/users/jmeyer/documents/github/proteomics/customMitoRos.txt",skip = 1,head=T,stringsAsFactors = F)
bOx<-read.delim("customATPsynth.txt",skip = 1,head=T,stringsAsFactors = F)

unlist(bOx)
tmp.index<-is.na(match(geneids,unlist(bOx)))

library(ggplot2)
library(gplots)
library(RColorBrewer)
tmpset
colnames(tmpset)

hsp90b1<-"P08113" ### HSP90b1
geneids
rownames(slf)

write.table(geneids,file="geneids.txt",sep="\t")


tmpset<-slf[tmp.index==FALSE,]
rownames(tmpset)<-paste(geneids[tmp.index==FALSE],substr(rownames(tmpset),start=8,stop=nchar(rownames(tmpset))))


box.aves<-group.aves(tmpset[,31:59], groups=NA)
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

bOx




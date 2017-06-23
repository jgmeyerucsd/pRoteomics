



library("Biobase")
library("annotate")
library("Category")
library("hgu95av2")
library("genefilter")

source("https://bioconductor.org/biocLite.R")
biocLite("GSEABase")




setwd("C:/users/jmeyer/documents/R24/sirt35overlap/")
read.csv(file="acK_1280sites_only.csv")->kac
kac<-read.delim(file="analysis_output.txt")
head(kac)
library(reshape2)
### cast acetyl
kac.castFC <- dcast(data = kac, formula = Protein ~ Label2,value.var=c("log2FC"))
kac.castFDR <- dcast(data = kac, formula = Protein ~ Label2,value.var=c("FDR"))
head(kac.castFC)
order<-c(2,3,4,5,6,7)
kac.fc.order<-kac.castFC[,order]
kac.fdr.order<-kac.castFDR[,order]

row.names(kac.fc.order)<-gsub(kac.castFC[,1],pattern=".K\\+42.0105",replacement = "")
row.names(kac.fdr.order)<-gsub(kac.castFC[,1],pattern=".K\\+42.0105",replacement = "")


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





setwd("Z:/Jesse/R24/finalsitelvl/")
pathway.files<-list.files(pattern="GSEA")




#### additional pathways


library(UniProt.ws)


geneids<-toupper(geneid)





##### combine the pathways into one ---

bOx<-read.delim(pathway.files[6],skip = 1,head=T,stringsAsFactors = F)
unlist(bOx)
tmp.index<-is.na(match(geneids,unlist(bOx)))

library(ggplot2)
library(gplots)
library(RColorBrewer)


tmpset<-kac.onlysig[tmp.index==FALSE,]
rownames(tmpset)<-paste(geneids[tmp.index==FALSE],substr(rownames(tmpset),start=8,stop=nchar(rownames(tmpset))))
?tiff
tiff(filename = "Rplot.tif",
     width = 8, height = 11, units = "in", pointsize = 12,
     compression = c("none"),
     bg = "white", res = 600, family = "", restoreConsole = TRUE,
     type = c("windows"))

my_palette <- colorRampPalette(c("blue","white", "red"))(n = 30)
heatmap.2(as.matrix(tmpset),Rowv = F,Colv = F,
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




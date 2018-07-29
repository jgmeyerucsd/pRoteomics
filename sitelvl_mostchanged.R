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
#### at least n changes
n=6
m.filtered<-m.onlysig[rowSums(is.na(m.onlysig))<=(ncol(m.onlysig)-n),]

hist(20-rowSums(is.na(m.filtered)))


nrow(m.onlysig)
nrow(m.filtered) ## 623 if using >2 fold, FDR<0.01
unique(substr(row.names(m.filtered),start=1,stop=6))  ### 251 proteins of those have at least one significant change in acylation site

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

m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", 
                tracecol="black", colsep=c(5,10,15), 
                sepcolor = "black",
                Rowv=as.dendrogram(m.clust), Colv=FALSE, 
                na.col="white",
                row.names<-"", 
                dendrogram="none",col=my_palette)
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
#write.table(t[m.hm$rowInd,],sep="\t",file="PIQED_hm2w10w_top5pct.txt")
write.table(t[m.hm$rowInd,],sep="\t",file="PIQED_hm2w10w_sKfix_manh.txt")


m<-read.delim(file="PIQED_hm2w10w_sKfix_manh.txt")
head(m)

uniprot<-substr(m[,1],start=1,stop=6)
#uniprot<-substr(row.names(t[m.hm$rowInd,]), start=1, stop=6)
m<-m[,2:ncol(m)]
head(m)




#### additional pathways
library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]


targs<-geneids

genemap <- select(up,unique(uniprot) , "GENES")


genes <- unlist(lapply(genemap[,2], split_sp))
geneids <- genes[match(uniprot, genemap[,"UNIPROTKB"])]


unique(geneid)
tmp.index<-is.na(match(geneids,targs))
tmpset<-m[tmp.index==FALSE,]

rownames(tmpset)

rownames(tmpset)<-paste(geneids[tmp.index==FALSE],substr(rownames(tmpset),start=8,stop=nchar(rownames(tmpset))))
mostchange<-tmpset
mostchange<-mostchange[order(row.names(mostchange)),]

top5pct.hm<-heatmap.2(as.matrix(mostchange[,c(6:10,16:20)]),scale="none", trace="none", tracecol="black", colsep=c(5,10,15), Rowv=F, Colv=FALSE, na.col="grey90",dendrogram="none",col=my_palette,cexRow = 0.4)









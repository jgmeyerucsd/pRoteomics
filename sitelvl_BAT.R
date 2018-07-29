setwd("C:/users/jmeyer/documents/R24/dietstudy/")



read.delim(file="Z:/BAT/BAT_acKsuccK_sitelvl_analysis_output.txt")->k
#read.delim(file="suck_analysis_output.txt")->s

head(k)


### need sites and separately the significant ones



library(reshape2)
### cast acetyl
k.castFC <- dcast(data = k, formula = Uniprot_site ~ Label2,value.var=c("log2FC"))
k.castFDR <- dcast(data = k, formula = Uniprot_site ~ Label2,value.var=c("FDR"))
#s.castFC <- dcast(data = s, formula = Protein ~ Label2,value.var=c("log2FC"))
#s.castFDR <- dcast(data = s, formula = Protein ~ Label2,value.var=c("FDR"))
head(k.castFC)
head(s.castFC)

#### order the columns so they are the order 2w vs control, 10w vs control, 10w vs 2w
#order2w<-c(14,15,16,17,18)
#order10w<-c(3,5,7,9,11)
#order2v10<-c(2,4,6,8,10,13)
#order<-c(order2w,order10w,order2v10) #all comparisons
#order<-c(order2w,order10w) # only 2w or 10w
#order<-order2v10 ###only 2v10

#### which sites are only regulated by sirt5 in cold acclimated 

uniprot<-substr(k.castFC[,1],start=1,stop=6)

head(k.castFC)


### find the sirt5 regulated sites in cold acclimated
ca.fc.index<-which(abs(k.castFC[,3])>1 & k.castFDR[,3]<0.01)
rt.fc.index<-which(abs(k.castFC[,5])>1 & k.castFDR[,5]<0.01)

intersect(ca.fc.index, rt.fc.index)
setdiff(ca.fc.index, rt.fc.index)-> only.ca
setdiff(rt.fc.index, ca.fc.index)-> only.rt

### find which proteins in these 

only.ca.prot<-uniprot[only.ca]
only.rt.prot<-uniprot[only.rt]

all.ca.prot<-uniprot[ca.fc.index]
all.rt.prot<-uniprot[rt.fc.index]

ca.uniq.prot<-setdiff(only.ca.prot,all.rt.prot)
rt.uniq.prot<-setdiff(only.rt.prot,all.ca.prot)

ca.rt.unique<-rbind(cbind(ca.uniq.prot,rep("#add8e6,black",times=length(ca.uniq.prot))),
  cbind(rt.uniq.prot,rep("red,black",times=length(rt.uniq.prot))))




write.table(file="Z:/BAT/sitelvl/only.ca.rt.sirt5reg.txt",ca.rt.unique, row.names = F,
            col.names=F, quote=F)

library(UniProt.ws)
up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
split_sp <- function(x)unlist(strsplit(x, " "))[1]

### find the genes only in cold ac
genemap <- select(up,ca.uniq.prot, "GENES")
genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]
geneid.filt <- genes[match(rownames(prot.filt), genemap[,"UNIPROTKB"])]

unlist(genemap[[2]])


x1<-c(1,2,3,4)
x2<-c(2,3,4,5)

setdiff(x1,x2)
setdiff(x2,x1)


?intersect
k.castFC[ca.fc.index,]


ca.fdr.index<-which(k.castFDR[,3]<0.01)




proteinchanges<-read.delim(file="Z:/BAT/20161222_BAT_protlvl_q10pct_candidates.xls")
head(proteinchanges)
pc<-proteinchanges[proteinchanges[,"Absolute.AVG.Log.Ratio"]>0.58,]
pc[,1]=="CAwt / CA5bko" | pc[,1]=="RT5bko / RTwt" 
pcs<-as.character(pc[pc[,1]=="CAwt / CA5bko" | pc[,1]=="RT5bko / RTwt" ,"UniProtIds"])
pcs


overlapping<-intersect(pcs,c(ca.uniq.prot,rt.uniq.prot))

match(overlapping,as.character(pc[,"UniProtIds"]))

proteinchanges[130,]


proteinchanges


temp<-cbind(rnorm(10),rnorm(10))
row.names(temp)<-c("Q04837","P0C0L4","P0C0L5","O75379","Q13068","A2MYD1","P60709","P30462","P30475","P30479")
colnames(temp)<-c("Exp1","Exp2")
temp
convertId(temp,filters="uniprot_swissprot",keepMultipleId=TRUE)










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
nrow(m.filtered) ## 623 if using >2 fold, FDR<0.01
unique(substr(row.names(m.filtered),start=1,stop=6))  ### 251 proteins of those have at least one significant change in acylation site

uniprot<-substr(row.names(m.filtered),start=1,stop=6)
m<-m.filtered
m[is.na(m)]<-0
m.dist<-dist(m,method = "manhattan")
#m.dist<-dist(m,method = "maximum")

m.clust<-hclust(m.dist)
?dist
m[m==0]<-NA
library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white","red"))(n = 11)
?heatmap.2
m.hm<-heatmap.2(as.matrix(m),scale="none", trace="none", 
                tracecol="black", colsep=c(5,10,16,21,26), 
                sepcolor = "black",
                Rowv=as.dendrogram(m.clust), Colv=FALSE, 
                na.col="white",
                row.names<-"", 
                dendrogram="none",col=my_palette,
                margins = c(10,5),
                labRow = c(""))
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
write.table(t[m.hm$rowInd,],sep="\t",file="PIQED_10v2_sKfix_manh.txt")
m<-read.delim(file="PIQED_hm2w10w_sKfix_manh.txt")
head(m)
na.omit(m)
lapply(na.omit,m)
for(x in 1:ncol(m)){
  print(length(na.omit(m[,x])))
}




uniprot<-substr(m[,1],start=1,stop=6)
uniprot<-substr(row.names(m),start=1,stop=6)

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
bOx<-read.delim("C:/users/jmeyer/documents/R24/dietstudy/sharedGenes.txt",skip = 1,head=T,stringsAsFactors = F)
bOx<-read.delim("AAmetab.txt",skip = 1,head=T,stringsAsFactors = F)
toupper(unlist(bOx))
tmp.index<-is.na(match(geneids,toupper(unlist(bOx))))
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

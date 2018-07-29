


setwd("C:/users/jmeyer/documents/KAT/")


list.files()

### load in gluC and trypsin data
read.delim(file="analysis_output_trypsin.txt")->k
read.delim(file="analysis_output_GluC.txt")->g
head(k)
##### Load the acK-regulated sites 

read.delim(file="acka_reg.txt")->a
head(a)
a<-a[,1:6]
a<-a[1:592,]

uniprot<-substr(a[,1],start=4,stop=9)
uniprot_site<-paste(uniprot,a[,5],sep="_")

### genes are better for finding overlap
gene_site<-paste(toupper(a[,"Gene.Name"]),a[,5],sep="_")
aa<-cbind(gene_site,uniprot_site,uniprot,a)

head(aa)
aa.unique.genes<-toupper(unique(aa[,"Gene.Name"]))
aa.unique.uniprot<-unique(aa[,"uniprot"])

kg.unique.genes<-unique(geneids)

intersect(aa.unique.uniprot, unique(uniprot))
intersect(aa.unique.genes, kg.unique.genes)

nrow(aa)
aa[,"Acetyl.Site"]


library(reshape2)
### cast acetyl
k.castFC <- dcast(data = k, formula = Protein ~ Label2,value.var=c("log2FC"))
k.castFDR <- dcast(data = k, formula = Protein ~ Label2,value.var=c("FDR"))
g.castFC <- dcast(data = g, formula = Protein ~ Label2,value.var=c("log2FC"))
g.castFDR <- dcast(data = g, formula = Protein ~ Label2,value.var=c("FDR"))
head(k.castFC)
head(g.castFC)

#### compare the uniprot values with each

#aa.unique.uniprot<-unique(aa[,"uniprot"])

#k.unique.uniprot<-unique(substr(k.castFC[,1],start=1,stop=6))
#g.unique.uniprot<-unique(substr(g.castFC[,1],start=1,stop=6))

#kg.unique.uniprot<-unique(c(k.unique.uniprot,g.unique.uniprot))
#kg.unique.uniprot<-uniqu
### get genes



#### order the columns so they are the order 2w vs control, 10w vs control, 10w vs 2w

### Birgit wants: Yfiq, YjaB, YiaC, RimI, PhnO
order1k<-c(4,6,5,3,2)
order1g<-c(3,5,4,2)



k.fc.order<-k.castFC[,order1k]
k.fdr.order<-k.castFDR[,order1k]
g.fc.order<-g.castFC[,order1g]
g.fdr.order<-g.castFDR[,order1g]

# make the row names match
row.names(k.fc.order)<-gsub(k.castFC[,1],pattern=".K\\+42.011",replacement = "")
row.names(k.fdr.order)<-gsub(k.castFC[,1],pattern=".K\\+42.011",replacement = "")
row.names(g.fc.order)<-gsub(g.castFC[,1],pattern=".K\\+42.011",replacement = "")
row.names(g.fdr.order)<-gsub(g.castFC[,1],pattern=".K\\+42.011",replacement = "")

#row.names(sl)<-gsub(sl[,1],pattern=".K\\+42.0105",replacement = "")

m.fc<-merge(k.fc.order,g.fc.order,all=T,by=0)
m.fdr<-merge(k.fdr.order,g.fdr.order,all=T,by=0)

head(m.fc)
head(m.fdr)

row.names(m.fc)<-m.fc[,1]
row.names(m.fdr)<-m.fdr[,1]
m.fc<-m.fc[,c(2:ncol(m.fc))]
m.fdr<-m.fdr[,c(2:ncol(m.fdr))]

nrow(m.fc) ###1367 sites

row.names(k.fc.order)
uniprot<-unique(substr(row.names(m.fc),start=1,stop=6))  ### 585 unique proteins with at least one mod quantified
unique(uniprot)
nrow(k.fc.order) ###1105 sites
unique(substr(row.names(k.fc.order),start=1,stop=6))  ### 586 unique proteins with at least one mod quantified
nrow(g.fc.order) ###381 sites
unique(substr(row.names(g.fc.order),start=1,stop=6))  ### 247 unique proteins with at least one mod quantified

#### k only significant
k.onlysig<-k.fc.order
k.onlysig[k.fdr.order>=0.01]<-NA
k.onlysig[k.onlysig<2]<-NA
k.onlysig[is.na(k.onlysig)]<-0   #### from NA to zero
k.onlysig[k.onlysig==0]<-NA  #### from zero to NA

### remove the rows that don't have all significant 
k.filtered<-k.onlysig[rowSums(is.na(k.onlysig))!=ncol(k.onlysig),]
nrow(k.filtered) ## 600 if using >2 fold, FDR<0.01
length(unique(substr(row.names(k.filtered),start=1,stop=6)))  ### 400 proteins of those have at least one significant change in acylation site


#### g only significant
g.onlysig<-g.fc.order
g.onlysig[g.fdr.order>=0.01]<-NA
g.onlysig[g.onlysig<2]<-NA
g.onlysig[is.na(g.onlysig)]<-0   #### from NA to zero
g.onlysig[g.onlysig==0]<-NA  #### from zero to NA

### remove the rows that don't have all significant 
g.filtered<-g.onlysig[rowSums(is.na(g.onlysig))!=ncol(g.onlysig),]
nrow(g.filtered) ## 169 if using >2 fold, FDR<0.01
length(unique(substr(row.names(g.filtered),start=1,stop=6)))  ### 426 proteins of those have at least one significant change in acylation site

head(m.fc)
#### combined only significant
m.onlysig<-m.fc
m.onlysig[m.fdr>=0.01]<-NA
m.onlysig[m.onlysig<2]<-NA
m.onlysig[is.na(m.onlysig)]<-0   #### from NA to zero
m.onlysig[m.onlysig==0]<-NA  #### from zero to NA

### remove the rows that don't have all significant 
m.filtered<-m.onlysig[rowSums(is.na(m.onlysig))!=ncol(m.onlysig),]
nrow(m.filtered) ## 706 if using >2 fold, FDR<0.01
length(unique(substr(row.names(m.filtered),start=1,stop=6)))  ### 434 proteins of those have at least one significant change in acylation site

###########################3

##### TO DO: here we need to setup the k/g merged table for merging with acka-regulated sites
row.names(m.filtered)


head(m)
#### Put here, naming the filtered rows by gene_site



############################

uniprot<-substr(row.names(m.filtered),start=1,stop=6)
m<-m.filtered
m[is.na(m)]<-0
m.dist<-dist(m,method = "manhattan")
#m.dist<-dist(m,method = "maximum")
m.clust<-hclust(m.dist)
m[m==0]<-NA

head(m)
colnames(m)

library(gplots)
library("RColorBrewer")
my_palette <- colorRampPalette(c("blue", "white","red"))(n = 51)

my_palette <- colorRampPalette(c("red","black"))(n = 21) ### with <0 excluded

#?heatmap.2

#############3 gene names from uniprot

#availableUniprotSpecies(pattern="K12")

#### annotate with gene names and protein descriptions

library(UniProt.ws)
split_sp <- function(x)unlist(strsplit(x, " "))[1]

up <- UniProt.ws(taxId=83333)
taxId(up) <- 83333



### genes
genemap <- select(up,unique(uniprot), "GENES")
genes <- unlist(lapply(genemap[,2], split_sp))
geneid <- genes[match(uniprot, genemap[,"UNIPROTKB"])]
geneids <- toupper(geneid)



### protein descriptions
protein.description.map<-select(up,keys=unique(uniprot),columns="PROTEIN-NAMES")
prot.descrs<-unlist(protein.description.map[,2])
descr.id <- prot.descrs[match(uniprot, protein.description.map[,"UNIPROTKB"])]


### get sequence protein descriptions
protein.sequence.map<-select(up,keys=unique(uniprot),columns="SEQUENCE")
protein.sequence.map[[2]][1]

prot.seq<-unlist(protein.sequence.map[,2])

#### loop throught the hits and retrieve their position
resplit<-strsplit(row.names(m), split="_")
lengths<-unlist(lapply(resplit,length))
uniprot<-substr(row.names(m),start=1,stop=6)
## need uniprot and site
uniprot

add.buffers=function(x)paste("XXXXXXXXXX",x,"XXXXXXXXXX",collapse="",sep="")
prot.seq.buf<-lapply(prot.seq,add.buffers)

split_num=function(x)unlist(strsplit(x, "_"))[2]
first.site<-as.numeric(unlist(lapply(row.names(m),split_num)))
seq.windows<-rep(NA,times=nrow(m))
for(i in 1:nrow(m)){
  if(lengths[i]==2){
    tmp.uni<-uniprot[i]
    tmp.site<-first.site[i]+10
    window<-c(tmp.site-10, tmp.site+10)
    
    tmp.seq.win<-substr(prot.seq.buf[match(tmp.uni, protein.sequence.map[,"UNIPROTKB"])],start=window[1],stop=window[2])
    seq.windows[i]<-tmp.seq.win
    }
  
}


write.table(cbind(seq.windows,m), file="merged.sig.motif.txt",sep="\t", quote=F)
cbind(seq.windows,m)
for(x in uniprot){
  prot.seq[match(x, protein.sequence.map[,"UNIPROTKB"])]

}


row.names(m)

ncol(protein.description.map)

row.names(m)<-paste(geneids,row.names(m),sep="_")
tmp.names<-strsplit(row.names(m),split="_")

odd<-seq(from=1,7,2)
new.names<-c()
for(i in 1:length(tmp.names)){
  new.names<-c(new.names,paste(na.omit(tmp.names[[i]][odd]),collapse="_",sep="_"))
  #typeof(na.omit(tmp.names[[i]][odd]))
  
}


row.names(m)<-new.names
colnames(m)<-c("YfiQ", "YjaB","YiaC", "RimI", "PhnO", "YfiQ.gluC", "YjaB.gluC","YiaC.gluC", "RimI.gluC")

intersect(row.names(m),aa[,1])

head(m)
colnames(m)<-


row.names(aa)<-aa[,1]
head(aa)
a2<-data.frame(aa[,"FC"])

row.names(a2)<-row.names(aa)
colnames(a2)<-c("acka")
head(a2)
ma<-merge(m,a2,by=0,all=T)
head(ma)

#acka.overlap.index<-match(row.names(m),row.names(aa))
#nrow(m)
#m2<-cbind(m, acka=rep(NA, times=nrow(m)))
#ncol(m2)
### loop through the acka match index, add the FC value to positions where there was overlap
#collector<-c()
#for( i in 1:nrow(m)){
#  if(is.na(acka.overlap.index[i])==FALSE){
#    print(row.names(m)[i])
 #   m2[i,ncol(m2)]<-aa[acka.overlap.index[i],"FC"]
#    collector<-c(collector,acka.overlap.index[i])
#  }
  
  
#}

#### finally add the sites that were not in the other set

#all.rows<-1:nrow(aa)
#collector
#unique.rows<-all.rows[-collector]

#natemplate<-rep(NA,times=ncol(m2)-1)
#m3<-m2
#for(i in 1:length(unique.rows)){
#  m3<-rbind(m3,c(natemplate,unique.rows[i]))
#}


#row.names(m)[37]
#row.names(aa)[534]

row.names(ma)<-ma[,1]
ma<-ma[,2:11]
ma[is.na(ma)]<-0
ma.dist<-dist(ma,method = "manhattan")
ma.clust<-hclust(ma.dist)
ma[ma==0]<-NA


head(m)
row.names(m)

nrow(m)


ma.hm<-heatmap.2(as.matrix(ma),scale="none", trace="none", 
                tracecol="black", colsep=c(5,9), 
                sepcolor = "black",
                Rowv=as.dendrogram(ma.clust), Colv=NA, 
                na.col="grey",
                dendrogram="none",col=my_palette,
                margins = c(10,10),cexRow = 0.2)




t<-cbind(descr.id,uniprot,m)
head(t)

t[m.hm$rowInd,]
output.table="Exclude_Negatives_include_acka_FC.txt"
write.table(ma[ma.hm$rowInd,],sep="\t",file=output.table)


########## Count the significant values found from each 
for( i in 1:ncol(m)){
  print(colnames(m)[i])
  print(length(na.omit(m[,i])))
}

head(t)

### function to count the proteins and sites for arbitrary list of sets
by.kat<-list(c(3,8),
             c(4,9),
             c(5,10),
             c(6,11),
             c(7)
               )


count.sitesandprots=function(table=t,setlist=by.kat){
  n.sets<-length(setlist)
  print("number of sets")
  print(n.sets)
  for( i in 1: n.sets){
    tmpset<-setlist[[i]]
    if(length(tmpset)>1){
      tmp.table<-table[rowSums(is.na(table[,tmpset]))!=ncol(table[,tmpset]),]
    }
    if(length(tmpset)==1){
      tmp.table<-table[is.na(table[,tmpset])==FALSE,]
    }
    gene_sites<-row.names(tmp.table)
    gene_sites_split<-strsplit(gene_sites,"_")
    line.n<-length(gene_sites_split)
    prots<-c(rep(0,times=line.n))
    sites<-c()
    for(j in 1:line.n){
      prots[j]<-gene_sites_split[[j]][1]
      sites<-c(sites,gene_sites_split[[j]][2:length(gene_sites_split[[j]])])
    }
    #print(i)
    print(colnames(table[,tmpset]))
    print("unique proteins")
    print(length(unique(prots)))
    print("total sites")
    print(length(sites))
  }
}

### set memberships for counting

by.ALL<-list(3:11)
by.protease<-list(c(3:7),c(8:11))
by.singles<-list(c(3),
                 c(4),
                 c(5),
                 c(6),
                 c(7),
                 c(8),c(9),c(10),c(11))
by.singles<-list(c(3,8),
                 c(4,9),
                 c(5,10),
                 c(6,11),
                 c(7))

### count significant sites + proteins
count.sitesandprots(setlist=by.ALL)
count.sitesandprots(setlist=by.singles)
count.sitesandprots(setlist=by.protease)

### count total identified 


k.fc.order

tmp.t<-k.fc.order
tmp.t<-g.fc.order
tmp.t<-m.fc

row.names(k.fc.order)
tmp.uni<-substr(row.names(tmp.t),start=1,stop=6)
row.names(tmp.t)<-paste(tmp.uni,row.names(tmp.t),sep="_")

tmp.names<-strsplit(row.names(tmp.t),split="_")

odd<-seq(from=1,7,2)
new.names<-c()
for(i in 1:length(tmp.names)){
  new.names<-c(new.names,paste(na.omit(tmp.names[[i]][odd]),collapse="_",sep="_"))
  #typeof(na.omit(tmp.names[[i]][odd]))
  
}


row.names(tmp.t)<-new.names
#colnames(tmp.t)<-c("YfiQ", "YjaB","YiaC", "RimI", "PhnO", "YfiQ.gluC", "YjaB.gluC","YiaC.gluC", "RimI.gluC")

gene_sites<-new.names
gene_sites_split<-strsplit(gene_sites,"_")
line.n<-length(gene_sites_split)
prots<-c(rep(0,times=line.n))
sites<-c()
for(j in 1:line.n){
  prots[j]<-gene_sites_split[[j]][1]
  sites<-c(sites,gene_sites_split[[j]][2:length(gene_sites_split[[j]])])
}
#print(i)
print(colnames(table[,tmpset]))
print("unique proteins")
print(length(unique(prots)))
print("total sites")
print(length(sites))

#### counting done
#######################################

### subsets for motifs
colnames(t)
setn<-c(3,8)
setn<-c(4,9)
setn<-c(5,10)
setn<-c(6,11)
nrow(t)
yfiq.only<-t[rowSums(is.na(t[,setn]))!=ncol(t[,setn]),]
nrow(yfiq.only)
#second

nrow(t[rowSums(is.na(t[,setn]))!=ncol(t[,setn]),])

#m<-read.delim(file="PIQED_hm2w10w_sKfix_manh.txt")
head(m)
head(yfiq.only)
na.omit(m)
lapply(na.omit,m)
for(x in 1:ncol(m)){
  print(length(na.omit(m[,x])))
}

write.table(file="rimI.subset.table.txt",yfiq.only, quote=F,sep="\t")


uniprot<-substr(m[,1],start=1,stop=6)
uniprot<-substr(row.names(m),start=1,stop=6)

m<-m[,2:ncol(m)]
head(m)




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

list.files()
yj.mot<-read.delim(file="yjab_motifs.txt",head=T,stringsAsFactors = F)
yi.mot<-read.delim(file="yiaC_motifs.txt",head=T,stringsAsFactors = F)
yf.mot<-read.delim(file="yfiq_motifs.txt",head=T, stringsAsFactors = F)

nchar(yj.mot[2,])

rm.k=function(motif.list=yj.mot){
  p1<-lapply(FUN=substr, motif.list,start=1,stop=10)
  p2<-lapply(FUN=substr, motif.list,start=12,stop=21)
  #lapply(FUN=paste, unlist(p1), unlist(p2))
  fixed<-paste(unlist(p1), unlist(p2),sep="")
  fixed
}

yi.fix<-rm.k(motif.list = yi.mot)
write.table(yi.fix,file="yiac.motif.krm.txt",sep="\t", row.names = F,quote=F)

### change mapDIA output to include gene name and protein description


setwd("C:/users/jmeyer/documents/KAT/")


### load in gluC and trypsin data
read.delim(file="analysis_output_trypsin.txt")->t
read.delim(file="analysis_output_GluC.txt")->g
head(t)

t<-cbind(uniprot.gene=gsub(t[,1],pattern=".K\\+42.011",replacement = ""),t)
g<-cbind(uniprot.gene=gsub(g[,1],pattern=".K\\+42.011",replacement = ""),g)

head(g)
nrow(g)
head(t)
t[1,]

t<-cbind(substr(t[,1],start=1, stop=6), t)
g<-cbind(substr(g[,1],start=1, stop=6), g)

head(g)
nrow(g)
head(t)


uniprot.t<-as.character(t[,1])
uniprot.g<-as.character(g[,1])
uniprot<-c(as.character(t[,1]),as.character(g[,1]))

library(UniProt.ws)
split_sp <- function(x)unlist(strsplit(x, " "))[1]
up <- UniProt.ws(taxId=83333)
taxId(up) <- 83333

### genes
?select
unique(uniprot)[1:50]
length(unique(uniprot))

genemap1 <- select(up,unique(uniprot)[1:100], "GENES")
genemap2 <- select(up,unique(uniprot)[101:200], "GENES")
genemap3 <- select(up,unique(uniprot)[201:300], "GENES")
genemap4 <- select(up,unique(uniprot)[301:400], "GENES")
genemap5 <- select(up,unique(uniprot)[401:500], "GENES")
genemap6 <- select(up,unique(uniprot)[501:668], "GENES")

genemap<-rbind(genemap1,genemap2, genemap3, genemap4, genemap5, genemap6)


typeof(genemap)
length(genemap)

genes <- unlist(lapply(genemap[,2], split_sp))
geneid.t <- genes[match(uniprot.t, genemap[,"UNIPROTKB"])]
geneid.g <- genes[match(uniprot.g, genemap[,"UNIPROTKB"])]

protein.description.map1<-select(up,keys=unique(uniprot)[1:100],columns="PROTEIN-NAMES")
protein.description.map2<-select(up,keys=unique(uniprot)[101:200],columns="PROTEIN-NAMES")
protein.description.map3<-select(up,keys=unique(uniprot)[201:300],columns="PROTEIN-NAMES")
protein.description.map4<-select(up,keys=unique(uniprot)[301:400],columns="PROTEIN-NAMES")
protein.description.map5<-select(up,keys=unique(uniprot)[401:500],columns="PROTEIN-NAMES")
protein.description.map6<-select(up,keys=unique(uniprot)[501:668],columns="PROTEIN-NAMES")

protein.description.map<-rbind(protein.description.map1,protein.description.map2,protein.description.map3, protein.description.map4, protein.description.map5, protein.description.map6)


prot.descrs<-unlist(protein.description.map[,2])
descr.id.t <- prot.descrs[match(uniprot.t, protein.description.map[,"UNIPROTKB"])]
descr.id.g <- prot.descrs[match(uniprot.g, protein.description.map[,"UNIPROTKB"])]

head(g)

g.out<-cbind(geneid.g,descr.id.g, g)
head(g.out)
write.table(file="GluC.allQuant.txt",g.out,quote=F, sep="\t")

t.out<-cbind(geneid.t, descr.id.t, t)

write.table(file="Trypsin.allQuant.txt",t.out,quote=F, sep="\t")





### change mapDIA output to include gene name and protein description


setwd("C:/BAT/")


### load in from written file data
read.delim(file="analysis_output_trypsin.txt")->regulated
reg<-s.set[hm1$rowInd,]

r.u<-cbind(substr(rownames(reg),start=1, stop=6),reg)

head(r.u)
nrow(r.u)


uniprot<-as.character(r.u[,1])


library(UniProt.ws)
split_sp <- function(x)unlist(strsplit(x, " "))[1]
?UniProt.ws
availableUniprotSpecies()

up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
genemap<-select(up,unique(uniprot), "GENES")


genes <- unlist(lapply(genemap[,2], split_sp))
geneid.reg <- genes[match(uniprot, genemap[,"UNIPROTKB"])]



protein.description.map<-select(up,keys=unique(uniprot),columns="PROTEIN-NAMES")

prot.descrs<-unlist(protein.description.map[,2])
descr.id.reg <- prot.descrs[match(uniprot, protein.description.map[,"UNIPROTKB"])]


reg.out<-cbind(geneid.reg,descr.id.reg, r.u)
head(reg.out)
getwd()
write.table(file="suK.reg.annotated.txt",reg.out,quote=F, sep="\t")


head(sites)
uniprot<-substr(sites[,1],start=1, stop=6)


genemap<-select(up,unique(uniprot), "GENES")


genes <- unlist(lapply(genemap[,2], split_sp))
geneid.sites <- genes[match(uniprot, genemap[,"UNIPROTKB"])]



protein.description.map<-select(up,keys=unique(uniprot),columns="PROTEIN-NAMES")

prot.descrs<-unlist(protein.description.map[,2])
descr.id.sites <- prot.descrs[match(uniprot, protein.description.map[,"UNIPROTKB"])]


sites.out<-cbind(geneid.sites,descr.id.sites, sites)
head(sites.out)
getwd()
write.table(file="BAT.acKsuK.unfiltered.allsites.annotated.txt",sites.out,quote=F, sep="\t")




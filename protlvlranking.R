


list.files()
setwd("Z:/sucla2/")


fibro<-read.delim(file="20171031_141337_fibro_qsparse_protlvl_pivotProtReport.xls",stringsAsFactors = F)
myo<-read.delim(file="20171031_141228_myo_protlvl_qsparse_pivotProtReport.xls",stringsAsFactors = F)

head(myo)
head(fibro)
fibro[,"PG.UniProtIds"]
range(fibro[,"PG.Qvalue"])


uniprot<-unique(substr(c(fibro[,"PG.UniProtIds"], myo[,"PG.UniProtIds"]),start=1,stop=6))
library(UniProt.ws)

up <- UniProt.ws(taxId=9606)
taxId(up) <- 9606
?select

uniprot.lengths <- select(up,uniprot, "LENGTH")
#uniprot.lengths<-genemap

ncol(uniprot.lengths)

### compute normalized quantity across all measures = sum(intensity)/length
m.quant<-grep(names(myo),pattern="Quantity")


norm.prot=function(table=myo, length.table=uniprot.lengths){
  quant.cols<-grep(names(table),pattern="Quantity")
  table.rows<-nrow(table)
  tmp<-data.frame(uniprot=as.character(substr(table[1,"PG.UniProtIds"],start=1,stop=6)),sum=as.character(sum(table[1,quant.cols])),stringsAsFactors = F)
  for(i in 1:table.rows){
    tmp<-rbind(tmp,c(as.character(substr(table[i,"PG.UniProtIds"],start=1,stop=6)),sum(table[i,quant.cols])))
    #tmp[i,]<-c(substr(table[i,"PG.UniProtIds"],start=1,stop=6),sum(table[i,quant.cols]))
  }
  head(tmp)
  ### divide each protein intensity by the length of that protein
  tmp.cor<-tmp
  for(i in 1:table.rows){
    tmp.len<-as.numeric(length.table[grep(length.table[,1],pattern=tmp[i,1]),2])
    tmp.cor[i,2]<-as.numeric(tmp[i,2])/tmp.len
  }
  hist(log(as.numeric(tmp.cor[,2]),base=10))
  tmp.cor<-tmp.cor[order(as.numeric(tmp.cor[,2]),decreasing=T),]
  tmp.cor
  
  
}


m<-norm.prot()
f<-norm.prot(table=fibro)


hist(log(as.numeric(myo.cor.abundance[,2]),base=10),main="myo.log.abundance",xlab=c("log10(corr. abund)"),breaks=20)
hist(log(as.numeric(f[,2]),base=10),main="fibro.log.abundance",xlab=c("log10(corr. abund)"),breaks=20)

myo.cor.abundance[,1]

sort(as.numeric(m[,2]),decreasing=T,index.return=T)



?sort
m.genes<- select(up,m[,1], "GENES")
f.genes<- select(up,f[,1], "GENES")
split_sp <- function(x)unlist(strsplit(x, " "))[1]

m.genes.first <- unlist(lapply(m.genes[,2], split_sp))
f.genes.first <- unlist(lapply(f.genes[,2], split_sp))

cbind(m.genes.first,f.genes.first)
write.table(m.genes.first,row.names = F, quote=F,file="myo.proteins.ordered.txt",sep="\t")
write.table(f.genes.first,row.names = F, quote=F,file="fibro.proteins.ordered.txt",sep="\t",col.names = F)


m.genes



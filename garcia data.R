
setwd("~/R24/dietstudy/revision data/isotopetracing/")


loc.tab<-read.delim("human_uniprot_subcellular.tab",stringsAsFactors = F)
head(loc.tab)

loc.tab<-cbind(loc.tab[,1],loc.tab[,ncol(loc.tab)-1])
head(loc.tab)

mito.list<-loc.tab[grep(loc.tab[,2],pattern="mitochondria"),]

mito.list[,1]
mito.list[782,]

garcia.data<-read.delim("garcia.txt",stringsAsFactors = F)
head(garcia.data)

is.mito<-match(garcia.data[,1],mito.list[,1])


garcia.data<-cbind(is.mito, garcia.data)
head(garcia.data)
mito.subset.garcia<-garcia.data[is.na(garcia.data[,1])==FALSE,]

mito.subset.garcia[,"MK.test.Glucose"]<=0.01

mito.sig.ss<-mito.subset.garcia[mito.subset.garcia[,"MK.test.Glucose"]<=0.01,]

mito.sig.ss
write.table(mito.sig.ss,file="mito.mksig.txt",sep="\t",quote=F, row.names = F)
write.table(mito.subset.garcia,file="mito.subset.all.txt",sep="\t",quote=F, row.names = F)

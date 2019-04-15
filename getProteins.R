setwd("P:/JGM_DI2A/20190405/FAIMS")
list.files(pattern="txt")

res.tab<-read.delim("1da10ppm_MCF7_FAIMS_18_2.txt",stringsAsFactors = F, header=T)
sorted.results<-res.tab[order(-res.tab$IonCount),]
head(results_byIonCount)

sorted.results<-res.tab[order(-res.tab[,"cosine"]),]
#peptide.factors<-as.factor(sorted.results[,"Peptide"])
t.first <- sorted.results[match(unique(sorted.results$Peptide), sorted.results$Peptide),]

head(t.first)
peptides<-t.first$Peptide

peptides_cleaned<-gsub(x=peptides, pattern="+[0-9]*.[0-9]", replacement = "")


getProteinList = function(sequences = peptides_cleaned,
                       fasta= "P:/JGM_DI2A/MSPLIT-DIAv1.0/2019-03-14-td-UP000005640.fasta"){
  require(seqinr)
  require(Biostrings)
  fas_ss<-readAAStringSet(filepath=fasta, format="fasta",
                       nrec=-1L, skip=0L, seek.first.rec=FALSE,
                       use.names=TRUE, with.qualities=FALSE)
  ### clean peptide string

  proteins<-list()
  i=1
  for( x in sequences){
    print(i)
    proteins[[x]] <- names(unlist(vmatchPattern(subject=fas, pattern=x, fixed=TRUE)))
    i=i+1
  }
  return(proteins)
}


### list apply form of assigning proteins
findpep=function(X, fas){names(unlist(vmatchPattern(subject=fas, pattern=X, fixed=TRUE)))}

op <- lapply(FUN=findpep, X=sequences, fas=fas_ss)

unlist_proteins <- unlist(op)


op[[1]]
op[[2]]
op[[3]]
op[[45]]

op

proteins < -getProteinList()
proteins[[1]]
proteins
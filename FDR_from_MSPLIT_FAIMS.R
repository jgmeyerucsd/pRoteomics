setwd("D:/20180816_FIA_DIA/results/")

list.files()

### gives error if 1% FDR not reached
peplvlfdr=function(msplitresults="msplit_out_pt8da.txt", fdrlevel=0.01){
  res.tab<-read.delim(msplitresults,stringsAsFactors = F, header=T)
  sorted.results<-res.tab[order(-res.tab[,"cosine"]),]
  #peptide.factors<-as.factor(sorted.results[,"Peptide"])
  t.first <- sorted.results[match(unique(sorted.results$Peptide), sorted.results$Peptide),]
  maxlines<-nrow(t.first)
  decoylines<-grep(t.first[,"Name"], pattern="DECOY")
  n.decoylines<-length(decoylines)
  fdr=0
  #i = 1
  lastdecoy<-c()
  ### loop through the decoy lines
  if(n.decoylines==0){lastdecoy<-NULL
  } else { 
    for(i in 1:n.decoylines){
      fdr<-i/decoylines[i]
      if(fdr>fdrlevel){
        lastdecoy<-c(lastdecoy,i)
      }
    }
  }
  if(is.null(lastdecoy)==TRUE){
    print("not enough decoys, accept all hits @ FDR=")
    print(n.decoylines/maxlines)
    print(paste("peptides =",maxlines-n.decoylines))
  }
  
  if(is.null(lastdecoy)==FALSE){
    i=min(lastdecoy)-1
    fdr<-i/decoylines[i]
    print(paste("fdr", round(length(decoylines[1:i])/decoylines[i], digits = 4)))
    print("score cutoff")
    print(t.first[decoylines[i],"cosine"])
    print(paste("peptide hits=", decoylines[i]-i))
  }
  ### make output
  #pep.output<-t.first[1:decoylines[i-2],]
  #pep.output
  
}





list.files()
peplvlfdr(msplitresults="msplitout_IT_1da_500ppmFrag.txt")
peplvlfdr(msplitresults="msplitout_IT_pt8da_500ppmFrag.txt")
pt4<-peplvlfdr(msplitresults="msplitout_IT_pt4da_500ppmFrag.txt")

peplvlfdr(msplitresults="msplit_out_pt4da.txt")


setwd("D:/20180816_FIA_DIA/results/")
oneda<-peplvlfdr(msplitresults="msplitout_IT_1da_500ppmFrag_pt8prec.txt")
pt8set<-peplvlfdr(msplitresults="msplitout_IT_pt8_500ppmFrag_pt4prec.txt")
pt4set<-peplvlfdr(msplitresults="msplitout_IT_pt4_500ppmFrag_pt4prec.txt")
p8p8<-peplvlfdr(msplitresults="msplitout_IT_pt8_500ppmFrag_pt4prec.txt")






getwd()
seqinr
library(seqinr)
require(Biostrings)


add.protein=function(msplitresults="D:/20180816_FIA_DIA/results/msplitout_IT_pt8da_500ppmFrag.txt",
                     fasta= "D:/20180816_FIA_DIA/2018-08-14-td-UP000002311.fas"){
  fas<-readAAStringSet(filepath=fasta, format="fasta",
                       nrec=-1L, skip=0L, seek.first.rec=FALSE,
                       use.names=TRUE, with.qualities=FALSE)
  
  ### loop through peptides, add the protein name to a vector
  table<-read.delim(msplitresults,stringsAsFactors = F)
  head(table)
  pepvec<-table[,"Peptide"]
  pepvec.cleaned<-gsub(x=pepvec, pattern="+[0-9]*.[0-9]", replacement = "")
  
  
  proteins<-rep(0, times= length(pepvec.cleaned))
  
  for(i in 1:length(pepvec.cleaned)){
    tmp.prot<-names(unlist(vmatchPattern(subject=fas, pattern=pepvec.cleaned[i], fixed=TRUE)))
    if(length(tmp.prot)==0){
      proteins[i]<-"DECOY"
    }
    if(length(tmp.prot)>0){
      proteins[i]<-tmp.prot
    }
    #proteins[i]<-names(unlist(vmatchPattern(subject=fas, pattern=pepvec.cleaned[i], fixed=TRUE)))
    print(i)
  }
  pep.pos <- lapply(file, function(x) {
    string <- BString(paste(x, collapse = ""))
    matchPattern(pattern, string, "naive exact")
  })
  peppos <- lapply(fas, function(x) {
    string <- BString(paste(x, collapse = ""))
    matchPattern(pattern, string)
  })
}

prot.fdr=function(fdrlevel=0.01,
                  msplitresults="D:/20180816_FIA_DIA/results/msplitout_IT_pt8da_500ppmFrag.txt",
                  fasta= "D:/20180816_FIA_DIA/2018-08-14-td-UP000002311.fas"){
  res.tab<-read.delim(msplitresults,stringsAsFactors = F, header=T)
  sorted.results<-res.tab[order(-res.tab[,"cosine"]),]
  #peptide.factors<-as.factor(sorted.results[,"Peptide"])
  t.first <- sorted.results[match(unique(sorted.results$Peptide), sorted.results$Peptide),]
  decoylines<-grep(t.first[,"Name"], pattern="DECOY")
  fdr=0
  i = 1
  while(fdr<fdrlevel){
    fdr<-length(decoylines[1:i])/decoylines[i]
    #cutoffscore <- t.first[decoylines[i],"cosine"]
    i=i+1
    print(fdr)
  }
  ### the single-pep version of the table
  pep.output<-t.first[1:decoylines[i-2],]
  
  ##### add protein names, determine protein-level FDR
  fas<-readAAStringSet(filepath=fasta, format="fasta",
                       nrec=-1L, skip=0L, seek.first.rec=FALSE,
                       use.names=TRUE, with.qualities=FALSE)
  
  ### loop through peptides, add the protein name to a vector
  
  pepvec<-pep.output[,"Peptide"]
  pepvec.cleaned<-gsub(x=pepvec, pattern="+[0-9]*.[0-9]", replacement = "")
  
  
  proteins<-rep(0, times= length(pepvec.cleaned))
  
  for(i in 1:length(pepvec.cleaned)){
    tmp.prot<-names(unlist(vmatchPattern(subject=fas, pattern=pepvec.cleaned[i], fixed=TRUE)))
    if(length(tmp.prot)==0){
      proteins[i]<-"DECOY"
    }
    if(length(tmp.prot)>0){
      proteins[i]<-tmp.prot
    }
    #proteins[i]<-names(unlist(vmatchPattern(subject=fas, pattern=pepvec.cleaned[i], fixed=TRUE)))
    print(i)
  }
  new.table<-cbind(proteins, pep.output)
  
  
  pdecoylines<-grep(new.table[,"proteins"], pattern="DECOY")
  
  pfdr=0
  pi = 1
  
  unique.proteins<-unique(new.table$proteins)
  #proteins
  
  while(pfdr<fdrlevel){
    n.unique.prot<-length(unique(proteins[seq(from=1, to=pdecoylines[pi])]))
    
    pfdr<-length(pdecoylines[1:pi])/n.unique.prot
    
    cutoffscore <- t.first[decoylines[pi],"cosine"]
    pi=pi+1
    print(pfdr)
  }
  
  
  
  
  ############### print stuff
  # peptide level
  print(paste("fdr", round(length(decoylines[1:(i-2)])/decoylines[i-2], digits = 4)))
  print("score cutoff")
  print(t.first[decoylines[i-2],"cosine"])
  print(paste("peptide hits=", decoylines[i-2]))
  
  # protein level
  print(paste("fdr", round((pi-2)/length(unique(proteins[seq(from=1, to=pdecoylines[pi-2])])), digits = 4)))
  print("prot score cutoff")
  print(new.table[pdecoylines[pi-2],"cosine"])
  print(paste("protein hits=", length(unique(proteins[seq(from=1, to=pdecoylines[pi-2])]))))
  
  ### make output
  return(list(  pep.output, new.table[1:pdecoylines[pi-2],]))
  
}
write.table(new.table[1:pdecoylines[pi-2],],file="protein.fdr.txt", sep="\t", row.names = F)
getwd()




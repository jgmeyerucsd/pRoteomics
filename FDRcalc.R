setwd("P:/JGM_DI2A/20190402_FAIMS_boudica/DI2A_conditionGrid/FAIMS/")

f<-list.files(pattern=".txt")

### gives error if 1% FDR not reached
peplvlfdr=function(msplitresults="20190402_MCF7_FAIMS_17_1.txt", fdrlevel=0.01){
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
    print(paste("not enough decoys, accept all",maxlines-n.decoylines,  "peptides @ FDR=", n.decoylines/maxlines, sep=" "))
    #print(n.decoylines/maxlines)
    #print(paste("peptides =",maxlines-n.decoylines))
    pep.output<-t.first[-decoylines,]
    return(pep.output)
  }
  
  if(is.null(lastdecoy)==FALSE & (min(lastdecoy)-1) !=0){
    i=min(lastdecoy)-1
    fdr<-i/decoylines[i]
    print(paste("fdr", round(length(decoylines[1:i])/decoylines[i], digits = 4)))
    print("score cutoff")
    print(t.first[decoylines[i],"cosine"])
    print(paste("peptide hits=", decoylines[i]-i))
    pep.output<-t.first[1:max(decoylines),]
    return(pep.output)
  }
  if(is.null(lastdecoy)==FALSE & (min(lastdecoy)-1) ==0){
    print(paste("FDR over", fdrlevel))
  }
  ### make output
  #pep.output
  }

f

peplvlfdr(msplitresults = f[13],fdrlevel = 0.01)

library(seqinr)
require(Biostrings)

f

msplitresults=f[5]

msplitresults=paste(getwd(), f[5], sep="/")

prot.fdr=function(fdrlevel=0.01,
                  msplitresults="D:/20180816_FIA_DIA/results/msplitout_IT_pt8da_500ppmFrag.txt",
                  fasta= "P:/JGM_DI2A/MSPLIT-DIAv1.0/2019-03-14-td-UP000005640.fasta"){
  res.tab<-read.delim(msplitresults,stringsAsFactors = F, header=T)
  sorted.results<-res.tab[order(-res.tab[,"cosine"]),]
  #peptide.factors<-as.factor(sorted.results[,"Peptide"])
  t.first <- sorted.results[match(unique(sorted.results$Peptide), sorted.results$Peptide),]
  decoylines<-grep(t.first[,"Name"], pattern="DECOY")
  fdr=0
  i = 1
  length(decoylines)
  if()
    while(fdr<fdrlevel){
      fdr<-length(decoylines[1:i])/decoylines[i]
      #cutoffscore <- t.first[decoylines[i],"cosine"]
      i=i+1
      print(fdr)
    }
  
  ### the single-pep version of the table
  pep.output<-t.first[1:decoylines[i-2],]
  ### if there are not enough decoys
  pep.output<-t.first
  
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
      print(length(tmp.prot))
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
  while(pfdr<fdrlevel){
    n.unique.prot<-length(unique(proteins[seq(from=1, to=pdecoylines[pi])]))
    pfdr<-length(pdecoylines[1:pi])/n.unique.prot
    cutoffscore <- t.first[decoylines[pi],"cosine"]
    pi=pi+1
    print(pfdr)
  }
  # protein level
  print(paste("fdr", round((pi-2)/length(unique(proteins[seq(from=1, to=pdecoylines[pi-2])])), digits = 4)))
  print("prot score cutoff")
  print(new.table[pdecoylines[pi-2],"cosine"])
  print(paste("protein hits=", length(unique(proteins[seq(from=1, to=pdecoylines[pi-2])]))))
  ### make output
  return(list(  pep.output, new.table[1:pdecoylines[pi-2],]))
}






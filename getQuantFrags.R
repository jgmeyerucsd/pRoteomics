getQuantFrags = function(sequence="A+42.01057AQGEPQVQFK"){
  require(MSnbase)
  hK<-8.014199 ## mass of heavy lysine
  hR<-10.00827 ## mass of heavy arg
  ### look for modificaiton masses in the peptides
  modvec<-c(C=57.02146)
  seq_mod<-gsub(x=sequence, pattern="[A-Z]", replacement= "")
  if(seq_mod=="+42.01057"){
    modvec<-c(C=57.02146, Nterm=42.01057)
  }
  seq_cleaned<-gsub(x=sequence, pattern="+[0-9]*.[0-9]", replacement = "")
  frags_light<-calculateFragments(seq_cleaned, modifications=modvec, neutralLoss=NULL)
  yfrags<-frags_light[frags_light$type=="y",]
  nfrags<-nrow(yfrags)
  z=1 ### only works for singly-charged fragments right now
  for(i in 1:nfrags){
    nk<- nchar(yfrags$seq[i]) - nchar(gsub("K", "", yfrags$seq[i])) # how many lysine?
    nr<- nchar(yfrags$seq[i]) - nchar(gsub("R", "", yfrags$seq[i])) # how many arg?
    new_mz <- yfrags$mz[i] + (hR/z)*nr+(hK/z)*nk
    yfrags<-rbind(yfrags, data.frame(mz=new_mz, 
                                     ion=yfrags$ion[i], 
                                     type="yheavy",
                                     pos=yfrags$pos[i],
                                     z=yfrags$z[i],
                                     seq=yfrags$seq[i]))
  }
  return(yfrags)
}
#source("https://bioconductor.org/biocLite.R")
#biocLite("mzR")
library(mzR)
#library(msdata)

#### gets a single spectra
###     to implement whole spectra matching process in R
####    change to get all spectra as R object

mgf.lib<-readLines(con="D:/FAIMS/msfragger_decoys_fixed.mgf")

get.mgf.spec=function(mgf=mgf.lib, 
                      specID="TITLE=MSfragger1.45180.45180")
  {
  tmp.index<-grep(pattern=specID, mgf)[1]+4
  mgf[tmp.index]#tmp.index
  first.line<-tmp.index+1
  while(mgf[tmp.index]!="END IONS")  {tmp.index=tmp.index+1}
  last.line<-tmp.index-1
  spectra <- as.numeric(unlist(strsplit(c(mgf[first.line:last.line]),split=" ")))
  mz<-spectra[seq(from=1, to=length(spectra), by=2)]
  int<-spectra[seq(from=2, to=length(spectra), by=2)]
  #plot(mz,int,type="h",lwd=1, xlim=c(200, 1200))  
  return(data.frame(mz,int))
}

### feed in output from spectra
filter.dia=function(spec=s18560,
                    lspec=lib.45180,
                    tol=0.3,
                    tol.type="ppm"){
  ### define the tolerances
  if(tol.type=="ppm"){
    low.mz=lspec[,"mz"]-lspec[,"mz"]*(tol/1000000)
    high.mz=lspec[,"mz"]+lspec[,"mz"]*(tol/1000000)
  }
  if(tol.type=="da"){
    low.mz=lspec[,"mz"]-tol
    high.mz=lspec[,"mz"]+tol
  }
  btw<-matrix(spec[spec[,1]>=low.mz[1] & spec[,1]<=high.mz[1]],byrow=F, ncol=2)
  #rbind(btw,btw)
  for(i in 2:length(low.mz)){
    btw<-rbind(btw,spec[spec[,1]>=low.mz[i] & spec[,1]<=high.mz[i]])
  }
 return(btw)
}

filtered.rawspec<-filter.dia()


#peaks(ms1da)
#s18560<-spectra(ms1da, scans=18560)


### plot all 3 ---- pt 4
### spec lib
mgf.lib<-readLines(con="C:/Users/jmeyer/Documents/msfragger_decoys_fixed.mgf")
### IT pt4
ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug22_JGM_FDIA2_IT_pt4_ol.mzXML")
### IT pt 8 
ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug15_JGM_FDIA2_IT_pt8_ol.mzXML")
### OT 1 Th
### 50k OT
ms1da<-openMSfile(filename="D:/FAIMS/20180822_DI2A/201808aug22_JGM_DIA2_OT_50k_3i1o_1.mzXML")


rawspec<-spectra(ms1da, scans=4984)
secondpep<-get.mgf.spec(mgf=mgf.lib,specID="MSfragger1.14753.14753")
### IT, 0.3 Da
#f.secondpep<-filter.dia(spec=rawspec, lspec=secondpep, tol=0.3, tol.type="da")
### OT, 10ppm
f.secondpep<-filter.dia(spec=rawspec, lspec=secondpep, tol=20, tol.type="ppm")

#dev.off()
#par(mfcol=c(3,1),cex.lab=1.2, cex.axis=1.2, mai=c(0.5,0.5,0.5,0))
#plot(rawspec[,1],rawspec[,2],type="h",lwd=1, xlim=c(200, 1200), xlab="mz", ylab="int", main="raw.dia")
#plot(secondpep[,1],secondpep[,2],type="h",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int", main="library spec")  
#plot(f.secondpep[,1],f.secondpep[,2],type="h",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int", main="filtered raw")  

# make mirrored library spectra and projected spectra
dev.off()
par(mfcol=c(2,1),cex.lab=1.2, cex.axis=1.2, mai=c(0.5,0.5,0.5,0.5))
plot(rawspec[,1],rawspec[,2],type="h",lwd=1, xlim=c(200, 1200), xlab="mz", ylab="int", main="raw.dia")
### normalize heights
normheights1<-secondpep[,2]/max(secondpep[,2])
normheights2<-f.secondpep[,2]/max(f.secondpep[,2])
plot(secondpep[,1],-normheights1,type="h",col="red",lwd=1, xlim=c(200, 1200),ylim=c(-1,1),xlab="mz", ylab="int", main="projected(top) vs. library(bot)")  
lines(f.secondpep[,1],normheights2,type="h",col="black",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int")  
abline(h=0)

library(OrgMassSpecR)
dev.off()

SpectrumSimilarity
SpectrumSimilarity(spec.top=data.frame(secondpep,normheights1), 
                   spec.bottom=data.frame(f.secondpep, normheights2), 
                   t = 0.11, b = 0.0001, 
                   top.label = NULL, bottom.label = "Library", 
                   xlim = c(200, 1200))


cos.sim <- function(ix) 
{
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
} 

ix<-matrix(c(secondpep[,1],normheights1),
           c(f.secondpep[,1],normheights2))
### compute the angle between the 2 spectra
theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )


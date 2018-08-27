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
ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug15_JGM_FDIA2_IT_pt4_ol.mzXML")
### IT pt 8 
ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug15_JGM_FDIA2_IT_pt8_ol.mzXML")
### OT 1 Th
ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug14_JGM_FDIA2_OT_1da.mzXML")


rawspec<-spectra(ms1da, scans=28405)
secondpep<-get.mgf.spec(mgf=mgf.lib,specID="MSfragger1.32000.32000")
### IT, 0.3 Da
f.secondpep<-filter.dia(spec=rawspec, lspec=secondpep, tol=0.3, tol.type="da")
### OT, 10ppm
f.secondpep<-filter.dia(spec=rawspec, lspec=secondpep, tol=10, tol.type="ppm")

dev.off()
par(mfcol=c(3,1),cex.lab=1.2, cex.axis=1.2, mai=c(0.5,0.5,0.5,0))
plot(rawspec[,1],rawspec[,2],type="h",lwd=1, xlim=c(200, 1200), xlab="mz", ylab="int", main="raw.dia")
plot(secondpep[,1],secondpep[,2],type="h",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int", main="library spec")  
plot(f.secondpep[,1],f.secondpep[,2],type="h",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int", main="filtered raw")  

# make mirrored library spectra and projected spectra
dev.off()
par(mfcol=c(2,1),cex.lab=1.2, cex.axis=1.2, mai=c(0.5,0.5,0.5,0.5))
plot(rawspec[,1],rawspec[,2],type="h",lwd=1, xlim=c(200, 1200), xlab="mz", ylab="int", main="raw.dia")
### normalize heights
normheights1<-secondpep[,2]/max(secondpep[,2])
normheights2<-f.secondpep[,2]/max(f.secondpep[,2])
plot(secondpep[,1],-normheights1,type="h",col="red",lwd=1, xlim=c(200, 1200),ylim=c(-1,1),xlab="mz", ylab="int", main="projected(top) vs. library(bot)")  
lines(f.secondpep[,1],normheights2,type="h",col="black",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int")  

# Lower minimum temperature in invasive range
tMin <- dat$Native_Temp_Min - dat$Invasive_Temp_Min
tMin[tMin < 0] <- 0
dat$TempMinIncrease <- tMin

# Higher maximum temperature in invasive range
tMax <- dat$Invasive_Temp_Max - dat$Native_Temp_Max
tMax[tMax < 0] <- 0
dat$TempMaxIncrease <- tMax


lib.45180

  









### first, add the lines with peptide to the mgf spec lib
##
?readLines
#mgf<-readLines(con="D:/jgm/MSfragger1.mgf")
#ms2<-readLines(con="D:/jgm/MSfragger1.ms2")

d.index.ms2<-grep(ms2, pattern="D")


modseq.index.ms2<-d.index.ms2[3:length(d.index.ms2)][seq(from=2,to= length(d.index.ms2), by =2)]

seqlines<-ms2[modseq.index.ms2]
seqlines
?gsub
seqlines<-gsub(seqlines, pattern="\\[", replacement="")
seqlines<-gsub(seqlines, pattern="\\]", replacement="")
seqlines<-gsub(seqlines, pattern="D\tmodified seq\t", replacement="SEQ=")


seqlines[1]


### where in mgf file are the 'PEPMASS' lines?
pepmass.mgf.index<-grep(mgf, pattern="PEPMASS")


unlink(outfile)

pepmass.mgf.index



txt <- "ampl.tab 2 1"
dat <- read.table(text = "
A  B
2   3
4   6
2   0
", header = TRUE)
tmp <- "tmp.txt"
getwd()
cat(txt, "\n", file = tmp) #" Don't forget the newline "\n"




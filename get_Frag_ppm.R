

f<-list.files(pattern=".txt")
f

fn<-"20190402_MCF7_FAIMS_17_1.txt"

ms1da<-openMSfile(filename="P:/JGM_DI2A/20190402_FAIMS_boudica/DI2A_conditionGrid/FAIMS/20190402_MCF7_FAIMS_17_1.mzXML")

po<-peplvlfdr(msplitresults = "20190402_MCF7_FAIMS_17_1.txt",fdrlevel = 0.01)


po$Name
scans<-unlist(lapply(FUN=function(x) x[[1]], strsplit(po$Name, split=" ")) )
scans

allppm<-c()
for(i in 1:nrow(po)){
  print(po$Scan.[i])
  rawspec<-spectra(ms1da, scans=po$Scan.[i])
  ls<-get.mgf.spec(mgf=mgf.lib,specID=scans[i])
  filtered<-filter.dia(spec=rawspec, lib_spec=ls, tol=30, tol.type="ppm")
  allppm<-c(allppm, filtered[[2]])
}
  


mean(allppm)
abline(v=mean(allppm))
hist(allppm[allppm>-10 & allppm<10])
mean(allppm[allppm>-10 & allppm<10])
sd(allppm[allppm>-10 & allppm<10])*3

rawspec<-spectra(ms1da, scans=677)
secondpep<-get.mgf.spec(mgf=mgf.lib,specID="TITLE=human.faims.61905.61905.")

### IT, 0.3 Da
#f.secondpep<-filter.dia(spec=rawspec, lspec=secondpep, tol=0.3, tol.type="da")
### OT, 10ppm
f.secondpep<-filter.dia(spec=rawspec, lib_spec=secondpep, tol=30, tol.type="ppm")


### first read the ID results from the 1:16 triplicate
f1 <- peplvlfdr(msplitresults = "2da10ppm_DI2A_1to16_tMS2_n1.txt",fdrlevel = 0.01)
f2 <- peplvlfdr(msplitresults = "2da10ppm_DI2A_1to16_tMS2_n2.txt",fdrlevel = 0.01)
f3 <- peplvlfdr(msplitresults = "2da10ppm_DI2A_1to16_tMS2_n3.txt",fdrlevel = 0.01)


### filter those for z=2

f1<-f1[f1$z.1==2,]
f2<-f2[f2$z.1==2,]
f3<-f3[f3$z.1==2,]



### get overlap between replicates
f1$Peptide

pepintersect<-intersect(intersect(f1$Peptide,f2$Peptide),f3$Peptide)
pepintersect
pf1<-f1[na.omit(match(f1$Peptide,pepintersect)),]
pf2<-f2[na.omit(match(f2$Peptide,pepintersect)),]
pf3<-f3[na.omit(match(f3$Peptide,pepintersect)),]

nrow(pf2)

head(pf1)
pf1$Scan.[1:10]
pf2$Scan.[1:10]

pf1$Peptide[1:10]
pf2$Peptide[1:10]

### do the scans also overlap?
tdf<-pf3
sp<-c()
for( i in 1:nrow(tdf)){
  sp<-c(sp,paste(tdf$Scan.[i],tdf$Peptide[i]))
}

sp1<-sp
sp2<-sp
sp3<-sp

pf1[pf1$Peptide=="LLADQAEAR",]
pf2[pf2$Peptide=="LLADQAEAR",]
pf3[pf3$Peptide=="LLADQAEAR",]



#####  Compute the singly-charged light y-ions for each peptide ###### 


??MSnbase
BiocManager::install("MSnbase")
a
library(MSnbase)


frags_light<-calculateFragments("LLADQAEAR", modifications=c(C=57.02146))
yfrags_light<-frags_light[frags_light$type=="y",]
yfrags_heavy




setdiff(sp1,sp2)
intersect(intersect(sp1,sp2),sp3)





list.files(pattern=".mzXML")


### get library spectra
#mgf.lib<-readLines(con="P:/JGM_DI2A/MSPLIT-DIAv1.0/human.faims.fixed.decoy.mgf")

get.mgf.spec=function(mgf=mgf.lib, 
                      specID="TITLE=MSfragger1.45180.45180")
{
  tmp.index<-grep(pattern=specID, mgf)[1]+5
  first.line<-tmp.index+1
  while(mgf[tmp.index]!="END IONS")  {tmp.index=tmp.index+1}
  last.line<-tmp.index-1
  spectra <- as.numeric(unlist(strsplit(c(mgf[first.line:last.line]),split=" ")))
  mz<-spectra[seq(from=1, to=length(spectra), by=2)]
  int<-spectra[seq(from=2, to=length(spectra), by=2)]
  #plot(mz,int,type="h",lwd=1, xlim=c(200, 1200))  
  return(data.frame(mz,int))
}



filtered[1,"Scan."]
filtered$Peptide[1]
unlist(strsplit(filtered$Name[1],split=" "))[1]


msfile<-openMSfile(filename="20190412_DI2A_1to16_tMS2_plus4p5_agc1e6_n1.mzXML")

rawspec<-spectra(msfile, scans=filtered[1,"Scan."])
secondpep<-get.mgf.spec(mgf=mgf.lib,specID=unlist(strsplit(filtered$Name[1],split=" "))[1])

specmatched<-filter.dia(spec=rawspec, lib_spec=secondpep, tol=30, tol.type="ppm")
f.secondpep<-specmatched[[1]]


plot(rawspec)
plot(rawspec[,1],rawspec[,2],type="h",lwd=1, xlim=c(200, 2000), xlab="mz", ylab="int", main="raw.dia")
normheights1<-secondpep[,2]/max(secondpep[,2])
normheights2<-f.secondpep[,2]/max(f.secondpep[,2])
plot(secondpep[,1],-normheights1,type="h",col="red",lwd=1, xlim=c(200, 2000),
     ylim=c(-1,1),xlab="mz", ylab="int", main="projected(top) vs. library(bot)")  
lines(f.secondpep[,1],normheights2,type="h",col="black",lwd=1, xlim=c(200, 2000),xlab="mz", ylab="int")  
abline(h=0)

par(mfcol=c(3,1))
plot(data.frame(mz,int), type="h", lwd=1, xlim=c(200, 2000),xlab="mz", ylab="int", main="library spec")
plot(secondpep[,1],-normheights1,type="h",col="red",lwd=1, xlim=c(200, 2000),
     ylim=c(-1,1),xlab="mz", ylab="int", main="projected(top) vs. library(bot)")  
plot(rawspec[,1],rawspec[,2],type="h",col="red",lwd=1, xlim=c(200, 2000),
    xlab="mz", ylab="int", main="rawspec")  
plot(f.secondpep[,1],normheights2,type="h",col="black",lwd=1, xlim=c(200, 2000),xlab="mz", ylab="int")














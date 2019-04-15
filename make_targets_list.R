
getwd()
setwd("P:/JGM_DI2A/SILAC")
list.files(pattern="txt")
### mass of heavy labels

filtered <- peplvlfdr(msplitresults = "2da10ppm_1to16_n1b.txt",fdrlevel = 0.01)

head(filtered)

filtered$Peptide[1]

filtered$Mz.1[1]
filtered$z.1[1]

peps<-filtered$Peptide


getHeavyMz = function(p=filtered$Peptide[1], 
                      mz = filtered$Mz.1[1], 
                      z = filtered$z.1[1]){
  hK<-8.014199 ## mass of heavy lysine
  hR<-10.00827 ## mass of heavy arg
  sk <- gsub("K","",p) ## replace any lysine with nothing
  sr <- gsub("R","",p) # replace any arginine with nothing
  nk <- nchar(p) - nchar(sk) # how many lysine?
  nr <- nchar(p) - nchar(sr) 
  new_mz <- mz+(hR/z)*nr+(hK/z)*nk
  return(new_mz)
}

#need columns: compound ==peptide_light, formula== NA, Adduct == H+, m/z, z, MSX ID (==i)

line1<-data.frame(Compound=filtered$Peptide[1], Formula="light", Adduct="H+", "m/z"= filtered$Mz.1[1], 
                  z= filtered$z.1[1], MSXID=1)
line2<-data.frame(Compound=filtered$Peptide[1], 
                  Formula="heavy", 
                  Adduct="H+", 
                  "m/z"= getHeavyMz(p=filtered$Peptide[1],
                                    mz = filtered$Mz.1[1],
                                    z = filtered$z.1[1]), 
                  z= filtered$z.1[1], MSXID=1)
df<-rbind(line1, line2)

#cur_line<-3
for( i in 2:nrow(filtered)){
  df<-rbind(df, data.frame(Compound=filtered$Peptide[i], Formula="light", Adduct="H+", "m/z"= filtered$Mz.1[i], 
                           z= filtered$z.1[i], MSXID=i))
  df<-rbind(df, data.frame(Compound=filtered$Peptide[i], 
                           Formula="heavy", 
                           Adduct="H+", 
                           "m/z"= getHeavyMz(p=filtered$Peptide[i],
                                             mz = filtered$Mz.1[i],
                                             z = filtered$z.1[i]), 
                           z= filtered$z.1[i], MSXID=i))
  print(i)
}


lightlines<-df[df$Formula=="light",]
heavylines<-df[df$Formula=="heavy",]


lightmz<-df[df$Formula=="light",]$m.z
heavymz<-df[df$Formula=="heavy",]$m.z

heavymz-lightmz

hist(as.numeric(lightlines$m.z), breaks=600)
hist(heavymz-lightmz, breaks=50)
unique(as.numeric(heavymz-lightmz))

mzdiff<-heavymz-lightmz

hist(unique(mzdiff), breaks=60)
typeof(mzdiff[0])

mzdiff
mzdiff0<-which(mzdiff==0)

lightlines[mzdiff0,]

#### make table where the co-isolated fragment is always +4, +4.5, or +5
diff=4.5
i=401
j=i+diff
msxid =1
df1<-data.frame(Compound="light", Formula="", Adduct="(no adduct)", m.z=i, z=2, msx_ID=msxid )
df1<-rbind(df1, data.frame(Compound="heavy", Formula="", Adduct="(no adduct)", m.z=j, z=2, msx_ID=msxid))

df1
while(i<1000){
  msxid=msxid+1
  i=i+1.5
  j=i+diff
  df1<-rbind(df1, data.frame(Compound="light", Formula="", Adduct="(no adduct)", m.z=i, z=2, msx_ID=msxid))
  df1<-rbind(df1, data.frame(Compound="heavy", Formula="", Adduct="(no adduct)", m.z=j, z=2,  msx_ID=msxid))
}


df1
write.table(df1, file="generaltargets_plus4p5.csv", sep=",",quote=FALSE, row.names = F)

write.table(df, file="generaltargets_plus5.csv", sep=",",quote=FALSE)


msfile<-openMSfile(filename=list.files(pattern="mzXML"))
??openMSfile
header1<-header(msfile, filtered$Scan.[1])
header1$

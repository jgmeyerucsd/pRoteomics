


setwd("C:/Users/jmeyer/Documents/MSfragger/")
list.files()

fr<-read.csv("frag_rep.csv", stringsAsFactors = F)


head(fr,10)

cvs<-c("30","40","50","60","70","80","90","100","110")


### list of column indexes for each CV
cvs.col.ind<-list()
for( x in cvs){
  cvs.col.ind[[x]]<-grep(x, colnames(fr))
}
cvs.col.ind


### list of which rows are IDed in each CV fraction
cvs.row.id.index<-list()
for( x in cvs){
  cvs.row.id.index[[x]]<-which(fr[,cvs.col.ind[[x]][2]]==TRUE)
}


colnames(fr)


### make list of the identified precursor m/z per fraction
cvs.prec.mz<-list()
for(x in cvs){
  cvs.prec.mz[[x]]<-fr[cvs.row.id.index[[x]],"Precursor.Mz"]
}


#y<-unique(cvs.prec.mz[[x]])[1]
### list of row indexes where those precursors are
### instead grab these on the fly per spec gen
cvs.prec.mz.row.index<-list()
for(x in cvs){
  tmp.index<-c()
  print(x)
  for(y in unique(cvs.prec.mz[[x]])){
    tmp.index<-c(tmp.index,which(fr[,"Precursor.Mz"]==y))
    }
  cvs.prec.mz.row.index[[x]]<- tmp.index
  #cvs.prec.mz.row.index[[x]]<- which(fr[,"Precursor.Mz"]==unique(cvs.prec.mz[[x]]))
}

cvs.prec.mz[[1]]

### make list of the fragment ions
cvs.frag.index<-which(fr[,"Fragment.Ion.Type"]!="precursor")

 


### make list of precursor isotope areas
cvs.prec.area<-list()
for(x in cvs){
  cvs.prec.area[[x]]<-fr[cvs.row.id.index[[x]],cvs.col.ind[[x]][1]]
}

#cvs.prec.area

n=0
dev.off()
pdf(file=paste("plot",n,".pdf",sep="", collapse = ""),width=9, height=24)
par(mfcol=c(9,1),cex=0.8)
for(x in cvs){
  #if(is.integer(n/3)){
  #  dev.off()
  #  tiff(filename=paste("plot",n,".tiff",sep="", collapse = ""),width=600, height=1440)
  #  }
  print(x)
  print(min(cvs.prec.mz[[x]]))
  print(max(cvs.prec.mz[[x]]))
  hist(cvs.prec.mz[[x]],
       xlab="m/z",
       breaks=seq(from=300,to=1250, by=0.7),
       main=x)
  n=n+1
  }
dev.off()


### smallest possible windows, overlap by 1/2 their width
q.bottom.slices<-seq(from=299.5, to= 1199.8, by =0.35)
q.top.slices<-seq(from=300.2, to = 1200.5, by = 0.35)

windows<-list()
for(i in 1:length(q.bottom.slices)){
  windows[[i]]<-c(q.bottom.slices[i],q.top.slices[i])
}
windows

dev.off()

### which precursors are fragments

frags.index<-list()
for(x in cvs){
  frags.index[[x]]<-intersect(cvs.prec.mz.row.index[[x]],cvs.frag.index)
}



### find exact examples
unique(cvs.prec.mz[[1]])[1:5]
length(windows)/20

x=4 ### CV==60

tmp.unique.prec<-unique(cvs.prec.mz[[x]]) ### CV==60
tmp.window<-windows[[1000]]  ### window is 649.15 to 649.85


tmp.prec.mz<-tmp.unique.prec[which( tmp.unique.prec >= tmp.window[1] & tmp.unique.prec <= tmp.window[2] )]


### get the frag masses for those precursors


tmp.frags<-c()

# y=649.3774
for(y in tmp.prec.mz){
  print(y)
  tmp.frag.indexes<-intersect(which(fr[,"Precursor.Mz"]==y), cvs.frag.index)
  print(length(tmp.frag.indexes))
  tmp.frags<-c(tmp.frags,tmp.frag.indexes)
  print(length(tmp.frags))
  }

sort(fr[tmp.frags,"Product.Mz"])
hist(fr[tmp.frags,"Product.Mz"],
     xlab="m/z",
     breaks=seq(from=200,to=1202, by=0.1),
     main=x)

colnames(fr)

length(tmp.prec.mz)*5


1200*0.000001




cvs.prec.mz.row.index[[1]][1:15]

fr[45:49,]






for(x in windows){
  
}





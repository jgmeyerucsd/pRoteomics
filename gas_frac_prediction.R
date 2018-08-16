


??s4




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



###### make ids master list with all the lists
pepids<-setClass("pepids", slots=c(skyline.all="ANY", subsets="list", subsets.col.index="ANY"))
pi<-new("pepids")
pi
pi@skyline.all<-fr
pi@subsets.col.index<-cvs.col.ind

#pi@skyline.all

#pi@subsets

### list of which rows are IDed in each CV fraction
#cvs.row.id.index<-list()

colnames(fr)
head(fr,20)

for(x in names(pi@subsets.col.index)){
  pi@subsets[[x]][["prec.id.indexes"]] <- which(pi@skyline.all[,cvs.col.ind[[x]][2]]==TRUE)
  pi@subsets[[x]][["prec.mzs"]] <- na.omit(pi@skyline.all[ pi@subsets[[x]] [["prec.id.indexes"]] , "Precursor.Mz" ] [seq(from=1,to=300000,by=3)])
  pi@subsets[[x]][["peptides"]] <- na.omit(pi@skyline.all[ pi@subsets[[x]] [["prec.id.indexes"]] , "Peptide" ] [seq(from=1,to=300000,by=3)])
  
  ### get the fragments now per peptide for those subsets
  ### by looping through the length of those prec.mz/peptide pairs, 
  ### taking all those column values from 4:length(retrieved rows)
  pi@subsets[[x]][["frag.mzs"]]<-list()
  for( i in 1: length(pi@subsets[[x]][["peptides"]]) ){
    tmp.index<-which(pi@skyline.all$Precursor.Mz == pi@subsets[[x]][["prec.mzs"]][i] & pi@skyline.all$Peptide == pi@subsets[[x]][["peptides"]][i])
    pi@subsets[[x]][["frag.mzs"]][[i]] <- pi@skyline.all$Product.Mz[tmp.index][4:length(tmp.index)]
    ### assign new list slot with same index as other two but called frag.mzs
    print(i)
  }
  
  print(x)
    #pi@subsets[[x]][["peptides"]]
}



### since there are groups of 3 precursor masses for each precursor followed by fragments
### the above takes all the triplets, then only every 3rd value out to 3 million
### and uses NA omit to clean off the trailing excess values

#length( pi@subsets[[x]][["peptides"]] )
#length( pi@subsets[[x]][["prec.mzs"]] )
#pi@subsets[[x]][["prec.mzs"]] [10000]
#pi@subsets[[3]][["dat"]]




colnames(fr)


### make list of the identified precursor m/z per fraction 
#cvs.prec.mz<-list()


  ### make list of the identified precursor m/z per fraction 

  
}

#pi@subsets[["110"]][["dat"]]


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

 
#### add list for all precursor m/z values for contrast
all.prec.z<-c()
for(x in names(pi@subsets.col.index)){
  all.prec.z<-c(all.prec.z, pi@subsets[[x]][["prec.mzs"]])
}
pi@subsets[["all"]][["prec.mzs"]]<-unique(all.prec.z)


### make list of precursor isotope areas
cvs.prec.area<-list()
for(x in cvs){
  cvs.prec.area[[x]]<-fr[cvs.row.id.index[[x]],cvs.col.ind[[x]][1]]
}

#cvs.prec.area

n=0
dev.off()
pdf(file=paste("all precursors pt35 iso",n,".pdf",sep="", collapse = ""),width=9, height=36)
par(mfcol=c(10,1),cex=0.8)
for(x in names(pi@subsets.col.index)){
  print(x)
  print(min(cvs.prec.mz[[x]]))
  print(max(cvs.prec.mz[[x]]))
  hist(pi@subsets[[x]][["prec.mzs"]],
       xlab="m/z",
       breaks=seq(from=300,to=1250, by=0.35),
       main=paste("CV=",x),
       ylim=c(0,100))
  n=n+1
}
hist(pi@subsets[["all"]][["prec.mzs"]],
     xlab="m/z",
     breaks=seq(from=300,to=1250, by=0.35),
     main="all precursors from all CVs",
     ylim=c(0,100))
dev.off()


### smallest possible windows, overlap by 1/2 their width
q.bottom.slices<-seq(from=399.5, to= 1199.8, by =0.35)
q.top.slices<-seq(from=400.2, to = 1200.5, by = 0.35)

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
length(windows)/60

x=4 ### CV==60
cvs[4]
tmp.unique.prec<-unique(cvs.prec.mz[[x]]) ### CV==60
tmp.window<-windows[[950]]  ### window is 649.15 to 649.85


tmp.prec.mz<-tmp.unique.prec[which( tmp.unique.prec >= tmp.window[1] & tmp.unique.prec <= tmp.window[2] )]


### get the frag masses for those precursors


tmp.frags<-c()

# y=649.3774

#### this gets all the precursors that match that window, which over estimates 
#### fix to only sample the list of current 'x' group (CV)
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
     breaks=seq(from=200,to=1202, by=0.3),
     main=paste("CV slice = ", cvs[x], ", DIA window = ", tmp.window[1], "-", tmp.window[2],sep=""))

colnames(fr)

length(tmp.prec.mz)*5


1200*0.000001




cvs.prec.mz.row.index[[1]][1:15]

fr[45:49,]






for(x in windows){
  
}





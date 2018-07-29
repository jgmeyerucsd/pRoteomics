getwd()

head(pl)
proteinlevels="proteinlevels.txt"


pl<-read.delim(proteinlevels,header = T,stringsAsFactors = F)
protlvl.proteins<-substr(pl[,1],start=1, stop=6)

#### match the site level columns to the correct protein level columns
pl.colnames<-names(pl)


namemap<-read.delim("name_mapping.txt",header = T,stringsAsFactors = F,colClasses = "character")
n.groups<-ncol(namemap)
groups<-colnames(namemap)
group.n<-c()
groupcol.index<-c()

for(x in groups){
  print(x)
  tmp.grp.nm<-namemap[,x][nchar(namemap[,x])>1]
  grp.cols<-grep(paste(tmp.grp.nm,sep="|",collapse = "|"), colnames(pl))
  i=1
  groupcol.index<-c( groupcol.index,grp.cols)
  for(y in grp.cols){
    
    colnames(pl)[y]<-paste(x,i,sep = ".")
    i=i+1
  }
  
  group.n<-c(group.n,i-1)
}


head(pl)[,32]
pl[,32]
finalout<-pl
### reorder to appropriate order
pl<-pl[,c(1,groupcol.index)]
names(pl)[32]
ncol(pl)
pl.10w<-pl[,c(32:61)]
rownames(pl.10w)<-pl[,1]

head(pl.10w)

pl.filt<-read.delim("Z:/Jesse/R24/protlvl/281candidates.onlysig.txt",header = T,stringsAsFactors = F)




pl.10w.filt<-pl.10w[na.omit(match(substr(rownames(pl.10w),start=1,stop=6),rownames(pl.filt))),]

nrow(pl.10w.filt)
nrow(pl.filt)
head(pl.10w.filt)
rownames(pl.10w.filt)

excludecontam<-pl.10w.filt[c(1:4,7,10,13:281),]

head(aaoa)
allgroups<-c("con2w","fr2w","gl2w","hfd2w","hfd_fr2w","hfd_gl2w","con10w","fr10w","gl10w","hfd10w","hfd_fr10w","hfd_gl10w")
row.def<-list(1:8,9:16,17:24,25:32,33:40,41:48,49:54,55:60,61:66,67:72,73:78,79:84)

row.names(table)
group.aves=function(table=aaoa, row.def=row.def, group.names=allgroups){
  if(length(row.def)<1){
    groups<-unique(substr(names(table),start=2,stop=nchar(names(table))-2))
  }
  if(length(row.def)>1){
    n.groups<-length(group.names)
    for(x in 1:n.groups){
      row.names(table)[row.def[[x]]]<-paste(group.names[x], row.names(table)[row.def[[x]]],sep="_")
    }
    table<-t(table)
    groups<-group.names
  }
  
  averaged<-data.frame(rep(0,times=nrow(table)))
  rownames(averaged)<-rownames(table)
  stdevs<-data.frame(rep(0,times=nrow(table)))
  rownames(stdevs)<-rownames(table)
  i=1
  for(x in groups){
    tmp.grp<-grep(colnames(table),pattern=x)
    if(i==1){
      averaged[,1]<-rowMeans(table[,tmp.grp])
      stdevs[,1]<-apply(table[,tmp.grp],1, sd)
      names(averaged)[i]<-groups[i]
      names(stdevs)[i]<-groups[i]
    }
    if(i>=2){
      averaged<-cbind(averaged,rowMeans(table[,tmp.grp]))
      stdevs<-cbind(stdevs,apply(table[,tmp.grp],1, sd))
      names(averaged)[i]<-groups[i]
      names(stdevs)[i]<-groups[i]
    }
    i=i+1
  }
  dev.off()
  par(mar = c(5, 6, 4, 5) + 0.1)
  row.names(averaged)
  n=7
  plotTop <- max(averaged[n,]) +
    stdevs[n,][averaged[n,] == max(averaged[n,])] * 3
  
  typeof(averaged)
  barplot(height = as.matrix(averaged),beside=T)
  
  barCenters <- barplot(height = as.matrix(averaged[n,]),
                        names.arg = colnames(averaged[n,]),
                        beside = T, las = 2,
                        ylim = c(0, plotTop),
                        cex.names = 0.75, xaxt = "n",
                        main = "citrate",
                        ylab = "umol",
                        border = "black", axes = TRUE)
  text(x = barCenters, y = par("usr")[3] - 0.05, srt = 45,
       adj = 1, labels = colnames(averaged[1,]), xpd = TRUE)
  
  segments(barCenters, as.numeric(averaged[n,] - stdevs[n,]), barCenters,
           as.numeric(averaged[n,] + stdevs[n,]), lwd = 1.5)
  
  arrows(barCenters, as.numeric(averaged[n,] - stdevs[n,]), barCenters,
         as.numeric(averaged[n,] + stdevs[n,]), lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
  
  head(averaged)
  head(stdevs)
  scaled<-t(scale(t(averaged)))
  head(scaled)
  sd(scaled[1,])
  library(gplots)
  library(RColorBrewer)
  ?brewer.pal
  mycolors<-colorRampPalette(c("blue","white","red"))(n=50)
  
  heatmap.2(scaled,scale="none",col = mycolors,trace="none",Colv = F,Rowv = T,dendrogram = "none",tracecol = "black",key.title = "")
  return(scaled)
}


group.aves(table=excludecontam)->subset







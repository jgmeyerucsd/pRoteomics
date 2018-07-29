prepMapDIAin(ptmProphName = "", 
                      skyline.output= "2016_0826_mapDIA.csv", 
                      ptm.score=0.99,
                      modstring= "K:42.0105",
                      wd="C:/users/jmeyer/documents/R24/piqed_sitelvl/",
                      namemapping=TRUE,
                      protlvl.correction=FALSE)
  
  
  


library(qvalue)

c.s3<-read.csv("C:/users/jmeyer/documents/R24/sirt35overlap/coonSirt3.csv",header = T)

c.s3[1,]
names(c.s3)

qtest<-qvalue(c.s3[,10])
qvalue(c.s3[,10])

plot(qvalue(c.s3[,10]))

plot(qtest$qvalue)

plot(c.s3[,11],qtest$qvalue)

c.s3[which(c.s3[,11]<=0.05),]->coon.p05
c.s3[which(c.s3[,11]<=0.01),]->coon.p01

nrow(coon.p05)
nrow(coon.p01)

head(coon.p01)
write.table(coon.p01,file="sirt3.coon.q01.tsv",sep="\t",quote=F,row.names=F)

getwd()
g.s3<-read.csv("C:/users/jmeyer/documents/R24/sirt35overlap/gibson.sirt3.csv",header = T)
head(g.s3)
filtered<-cbind(g.s3[,c(1:11)],qvalue(g.s3[,"p.Value"])$qvalue)
gibson.p01<-filtered[which(filtered[,12]<=0.01),]
write.table(gibson.p01,file="C:/users/jmeyer/documents/R24/sirt35overlap/sirt3.gibson.q01.tsv",sep="\t",quote=F,row.names=F)

site.n.gibson<-gsub(x=gibson.p01[,"Acetyl.Site"],", ",replacement = "_")

acc_site.gibson<-paste(gibson.p01[,"Accession."],site.n.gibson,sep="_")

write.table(unique(acc_site.gibson),file="C:/users/jmeyer/documents/R24/sirt35overlap/gibson_siteonly.tsv",sep="\t",row.names = F)

### make the sites from coon same format
pos.n<-gsub(coon.p01[,"AcetylSite"],pattern="K",replacement = "")



pos.n2<-gsub(pos.n,pattern="|",replacement=" ")
cbind(paste(coon.p01[,2],pos.n,sep = "_"),coon.p01)


coon.sites<-paste(coon.p01[,2],pos.n,sep = "_")
gib.sites<-acc_site.gibson


intersect(coon.sites,gib.sites)
both.sites<-unique(c(coon.sites,gib.sites))

getwd()
setwd()
write.table(both.sites,file="sirt3.coon.gibson.tsv",sep="\t",quote=F,row.names = F)
both.sites<-read.csv("C:/users/jmeyer/documents/R24/dietstudy/sirt3.coon.gibson.tsv",header = F,stringsAsFactors = F)
both.sites
setwd("C:/users/jmeyer/documents/R24/dietstudy/")

ac<-read.delim(file="mergedsites_PIQED_stringent_4overlap.txt",row.names = 1)
head(ac)
ac<-m
?gsub
gsub("[:punct:]","_",as.character(both.sites),fixed = T)

typeof(both.sites)
unlist(both.sites)
overlap<-intersect(unlist(both.sites),row.names(ac))
overlap
overlap.rows<-match(overlap,row.names(ac))
sirt3ol<-rep(0,times=nrow(ac))
sirt3ol[match(overlap,row.names(ac))]<-1


s5<-read.delim(file="sirt5_succinyl.txt",stringsAsFactors = F)
head(s5)
s5qtest<-qvalue(s5[,"pvalue"])
hist(s5qtest$qvalues)
hist(s5[,"pvalue"])
s5q<-cbind(s5,qvalue(s5[,"pvalue"])$qvalues)
head(s5q)
colnames(s5q)[ncol(s5q)]<-"q"
s5q.sig<-s5q[s5q[,"q"]<=0.01,]
s5p.sig<-s5q[s5q[,"pvalue"]<=0.01,]
head(s5p.sig)
head(s5q.sig)

nrow(s5p.sig)
nrow(s5q.sig)

write.table(paste(s5q.sig[,1],s5q.sig[,"site"],sep="_"),file="sirt5reg_siteonly.tsv",sep="\t",row.names = F,quote=F)
s5reg<-read.table("sirt5reg_siteonly.tsv",stringsAsFactors = F)

overlap<-intersect(unlist(s5reg),row.names(ac))
overlap
overlap.rows<-match(overlap,row.names(ac))
sirt5ol<-rep(0,times=nrow(ac))
sirt5ol[match(overlap,row.names(ac))]<-1
sirt3ol
sirt5ol
both35ol.tmp<-sirt3ol+sirt5ol

both35ol<-rep(0,times=nrow(ac))
both35ol[both35ol.tmp==2]<-1
sirt3ol[both35ol.tmp==2]<-0
sirt5ol[both35ol.tmp==2]<-0

sum(both35ol)
sum(sirt3ol)
sum(sirt5ol)

both35ol+sirt3ol+sirt5ol



#kac.onlysig[641,]
#overlap[1]
kac.ol<-cbind(ac[,1:10],sirt3ol,sirt5ol,both35ol)
kac.ol<-cbind(ac[,11:20],sirt3ol,sirt5ol,both35ol)

head(kac.ol)

#### loop through columns

up<-c()
down<-c()
s3.up<-c()
s3.down<-c()
s5.up<-c()
s5.down<-c()
b.up<-c()
b.down<-c()
kac.ol[is.na(kac.ol)]<-0


for(x in 1:10){
  up.tmp<-0
  down.tmp<-0
  s3.up.tmp<-0
  s3.down.tmp<-0
  s5.up.tmp<-0
  s5.down.tmp<-0
  b.up.tmp<-0
  b.down.tmp<-0
  
  for(y in 1:nrow(kac.ol)){
    if(kac.ol[y,x]>0 & sum(kac.ol[y,11:13])==0){up.tmp=up.tmp+1} 
    if(kac.ol[y,x]<0 & sum(kac.ol[y,11:13])==0){down.tmp=down.tmp+1} 
    if(kac.ol[y,x]>0 & kac.ol[y,11]==1){s3.up.tmp=s3.up.tmp+1} 
    if(kac.ol[y,x]<0 & kac.ol[y,11]==1){s3.down.tmp=s3.down.tmp+1}
    if(kac.ol[y,x]>0 & kac.ol[y,12]==1){s5.up.tmp=s5.up.tmp+1} 
    if(kac.ol[y,x]<0 & kac.ol[y,12]==1){s5.down.tmp=s5.down.tmp+1}
    if(kac.ol[y,x]>0 & kac.ol[y,13]==1){b.up.tmp=b.up.tmp+1} 
    if(kac.ol[y,x]<0 & kac.ol[y,13]==1){b.down.tmp=b.down.tmp+1}
    
  }
  up<-c(up,up.tmp)
  down<-c(down,down.tmp)
  s3.up<-c(s3.up,s3.up.tmp)
  s3.down<-c(s3.down,s3.down.tmp)
  s5.up<-c(s5.up,s5.up.tmp)
  s5.down<-c(s5.down,s5.down.tmp)
  b.up<-c(b.up,b.up.tmp)
  b.down<-c(b.down,b.down.tmp)
  #check that sum is the same for components
  #table(kac.ol[,x]>0)
  #table(kac.ol[,x]<0)
  
  }
groups<-c("2wFr",
          "2wGl",
          "2wHFD",
          "2wHFD+Fr",
          "2wHFD+Gl",
          "10wFr",
          "10wGl",
          "10wHFD",
          "10wHFD+Fr",
          "10wHFD+Gl")

### in order of layer drawing where each layer is the sum of the remaining layetrs

totalup<-up+s3.up+s5.up+b.up
t2.up<-s3.up+s5.up+b.up
t3.up<-s3.up+s5.up
t5.up<-s5.up

totaldown<-down+s3.down+s5.down+b.down
t2.d<-s3.down+s5.down+b.down
t3.d<-s3.down+s5.down
t5.d<-s5.down

groups<-c(rep("up",times=10),rep("down",times=10))
s3.groups<-c(rep("s3.up",times=10),rep("s3.down",times=10))
s5.groups<-c(rep("s5.up",times=10),rep("s5.down",times=10))
b.groups<-c(rep("b.up",times=10),rep("b.down",times=10))

dat<-data.frame(groups,up,down,b.up,b.down,s3.up,s3.down,s5.up,s5.down)

dat<-data.frame(c(totalup,-totaldown,t2.up,-t2.d,t3.up,-t3.d,t5.up,-t5.d))
dat<-cbind(c(groups,b.groups,s3.groups,s5.groups),dat)
dat<-cbind(c(1:10,1:10,1:10,1:10),dat)
dat<-cbind(c(11:20,11:20,11:20,11:20),dat)

head(dat)
colnames(dat)<-c("x","groups","y")
ds<-dat
da<-dat

datb<-rbind(da,ds)

s.down
library(ggplot2)
plot(height=matrix(up,down,s.up,s.down))
group.colors<-c("blue", "#33CC33", "magenta","red","blue", "#33CC33", "magenta","red")
group.colors<-c("blue", "blue", "orange","red", "red","#33CC33", "#33CC33","orange")
?svg
dev.off()
tiff(filename = "Overlap_sfix1.tif",
     width = 4, height = 5, units = "in", pointsize = 1,
     compression = c("none"),
     bg = "white", res = 600, family = "", restoreConsole = TRUE,
     type = c("windows"))
svg(filename = "Overlap_sfix1.svg",
   width = 7, height = 7, pointsize = 12,
   onefile = FALSE, family = "sans", bg = "white",
   antialias = c("default", "none", "gray", "subpixel"))


ggplot(datb, aes(x=x, y=y, fill=groups)) + 
  geom_bar(stat="identity", position="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.y=element_blank()) +
  geom_abline(slope=0,intercept = 0, lwd=2) +
  scale_fill_manual(values=group.colors)
dev.off()
?geom_bar


colnames(kac.ol)
kac.ol[,1]
getwd()

getwd()

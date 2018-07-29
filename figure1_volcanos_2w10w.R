


setwd("C:/users/jmeyer/documents/R24/dietstudy/")


#### read site levels and remove anything but time dep.  in excel
read.delim(file="ak_analysis_output_2w10wonly.txt")->k
read.delim(file="sk_analysis_output_2w10only.txt")->s


head(k)

### read protlvl and flip upside down rows
p<-read.delim(file="Z:/R24/resource/proteinlvl/2w10w_prot_candidates_unfilt.txt")
head(p)
p.flip<-p
#p[,"Condition.Numerator"]=="10w_chow"
p[,"Condition.Denominator"]!="control"

p.flip[p[,"Condition.Denominator"]!="control","AVG.Log2.Ratio"]<-p[p[,"Condition.Denominator"]!="control","AVG.Log2.Ratio"]*-1
head(p)
head(p.flip)
### split into 2w/10w



p2<-p.flip[grep(p[,"Comparison..group1.group2."],pattern="2w"),]
p10<-p.flip[grep(p[,"Comparison..group1.group2."],pattern="10w"),]



#### p.flip should be fixed
plot.p2<-cbind(p2[,"AVG.Log2.Ratio"],-log(p2[,"Qvalue"],base=10))
plot.p10<-cbind(p10[,"AVG.Log2.Ratio"],-log(p10[,"Qvalue"],base=10))
colnames(plot.p2)<-c("FC","FDR")
colnames(plot.p10)<-c("FC","FDR")

p2.blue<-which(plot.p2[,1]<(-0.58) & plot.p2[,2]>2)
p2.red<-which(plot.p2[,1]>0.58 & plot.p2[,2]>2)
p10.blue<-which(plot.p10[,1]<(-0.58) & plot.p10[,2]>2)
p10.red<-which(plot.p10[,1]>0.58 & plot.p10[,2]>2)

### need sites and separately the significant ones

head(k)
k2<-k[grep(k[,"Label2"], pattern="2w"),]
k10<-k[grep(k[,"Label2"], pattern="10w"),]
# acetyl 2w
plot.k2<-cbind(k2[,"log2FC"],-log(k2[,"FDR"],base=10))
colnames(plot.k2)<-c("FC","FDR")
k2.blue<-which(plot.k2[,1]<(-1 )& plot.k2[,2]>2)
k2.red<-which(plot.k2[,1]>1 & plot.k2[,2]>2)

# acetyl 10w
plot.k10<-cbind(k10[,"log2FC"],-log(k10[,"FDR"],base=10))
colnames(plot.k10)<-c("FC","FDR")
k10.blue<-which(plot.k10[,1]<(-1 )& plot.k10[,2]>2)
k10.red<-which(plot.k10[,1]>1 & plot.k10[,2]>2)



s2<-s[grep(s[,"Label2"], pattern="2w"),]
s10<-s[grep(s[,"Label2"], pattern="10w"),]
plot.s2<-cbind(s[,"log2FC"],-log(s[,"FDR"],base=10))
plot.s10<-cbind(s[,"log2FC"],-log(s[,"FDR"],base=10))

colnames(plot.s2)<-c("FC","FDR")
colnames(plot.s10)<-c("FC","FDR")

s2.red<-which(abs(plot.s2[,1])>1 & plot.s2[,2]>2)
s2.blue<-which(plot.s2[,1]<(-1) & plot.s2[,2]>2)
s10.red<-which(abs(plot.s10[,1])>1 & plot.s10[,2]>2)
s10.blue<-which(plot.s10[,1]<(-1) & plot.s10[,2]>2)

###############################################33
######## plot everything
###########################################
dev.off()
par(mfcol=c(3,2))
#prot 2w
plot(plot.p2[,"FC"],plot.p2[,"FDR"],pch=20,
     xlim=c(-5,5),bty="n",xlab="log2FC",ylab="-log10(FDR)")
points(plot.p2[p2.red,],col="red",pch=20)
points(plot.p2[p2.blue,],col="blue",pch=20)
#acetyl 2w
plot(k2[,"log2FC"],-log(k2[,"FDR"],base=10),pch=20,
     xlim=c(-5,5),bty="n",xlab="log2FC",ylab="-log10(FDR)")
points(plot.k2[k2.red,],col="red",pch=20)
points(plot.k2[k2.blue,],col="blue",pch=20)
# succinyl 2w
plot(s2[,"log2FC"],-log(s2[,"FDR"],base=10),
     xlim=c(-5,5),pch=20,bty="n",xlab="log2FC",ylab="-log10(FDR)")
points(plot.s2[s2.red,],col="red",pch=20)
points(plot.s2[s2.blue,],col="blue",pch=20)
#10w prot
plot(plot.p10[,"FC"],plot.p10[,"FDR"],pch=20,
     xlim=c(-5,5),bty="n",xlab="log2FC",ylab="-log10(FDR)")
points(plot.p10[p10.red,],col="red",pch=20)
points(plot.p10[p10.blue,],col="blue",pch=20)
### 10w acetyl
plot(k10[,"log2FC"],-log(k10[,"FDR"],base=10),pch=20,
     xlim=c(-5,5),bty="n",xlab="log2FC",ylab="-log10(FDR)")
points(plot.k10[k10.red,],col="red",pch=20)
points(plot.k10[k10.blue,],col="blue",pch=20)
#10w succinyl
plot(s10[,"log2FC"],-log(s10[,"FDR"],base=10),
     xlim=c(-5,5),pch=20,bty="n",xlab="log2FC",ylab="-log10(FDR)")
points(plot.s10[s10.red,],col="red",pch=20)
points(plot.s10[s10.blue,],col="blue",pch=20)

plot.p<-

#####

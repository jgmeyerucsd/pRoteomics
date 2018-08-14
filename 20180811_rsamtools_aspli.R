
###### event poiinter

source("http://www.bioconductor.org/biocLite.R")
biocLite("EventPointer")

source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment")
biocLite("spliceR")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocGenerics", version = "devel")
BiocManager::install("GenomicRanges", version = "devel")
BiocManager::install("S4Vectors", version = "devel")

a
BiocManager::install("GenomicRanges", version = "devel",update=T,lib="C:/Users/jmeyer/Documents/")
?BiocManager::install
library(GenomicRanges)
library(SGSeq)

source("https://bioconductor.org/biocLite.R")

#biocLite("Rsamtools")
library(Rsamtools)

getwd()

setwd("D:/20180718_RNAseq/BAM/")
files<-list.files(,pattern="bam")
files

x<-files[1]
for(x in files){
  print(x)
}
for(x in files){
  print(x)
  sortBam(x, paste(x,".sort",collapse = "",sep=""),maxMemory=8000)
  indexBam(paste(x,".sort.bam",collapse = "",sep=""),paste(x,".sort.bam.i",collapse = "",sep=""))
}

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")

library(GenomicFeatures)

annFile <- aspliExampleGTF()
aTxDb <- makeTxDbFromGFF(annFile)
> # extract features from annotation

> # Accesors of ASpliFeatures object
geneCoord <- featuresg( features )
binCoord <- featuresb( features )
junctionCoord <- featuresj( features )
# Extract metadata annotation, such as bins names, event and feature type
binMetadata <- mcols( binCoord )


### splicing analysis
source("https://www.bioconductor.org/biocLite.R")
biocLite("ASpli")
a
library(ASpli)
getwd()
setwd("D:/20180718_RNAseq/BAM/sorted/")
setwd("C:/RNAseq/")
files<-list.files(pattern="bam")
files
odd

#### get  features
TxDb <- makeTxDbFromGFF(
  file="D:/RNAseq/Bam_Data/Caenorhabditis_elegans.refGene.ce10.gtf",
  format="gtf")
features <- binGenome( TxDb )




bamFiles<-files[seq(1,23, by=2)]


bam <- loadBAM(targets, cores=1)
?readCounts

sit<-getBamInfo(sample_info=targets,yieldSize=NULL)

##### this function gives an ERROR!!!
counts <- readCounts (
  features,
  bam,
  cores = 1,
  readLength = 125,
  maxISize = 50000,
  minAnchor = NULL )

getwd()
?getBamInfo

setwd("D:/RNAseq/Bam_Data/")





importTranscripts()

#### remapped reads using hisat2
setwd("C:/RNAseq/")
si<-data.frame(sample_name=c("G1","G2","G3","G4","G5","N1","N2","N3","N4","N5"),file_bam=list.files(pattern="bam"),stringsAsFactors = F)

sit<-getBamInfo(sample_info=targets,yieldSize=NULL)

?importTranscripts
importTranscripts(file="Caenorhabditis_elegans.WBcel235.90.gff3")


sgf_ucsc <- convertToSGFeatures(txf_ucsc)
head(sgf_ucsc)

sgfc_ucsc <- analyzeFeatures(sit, features = txf_ucsc)
library(GenomicFeatures)
samplefile <- system.file("C:/RNAseq/ce10.sqlite",
                          package="GenomicFeatures")
samplefile <- system.file("C:/RNAseq/ce10.sqlite")


samplefile

txdb <- loadDb("C:/RNAseq/ce10.sqlite")



seqlevelsStyle(txdb)<-"NCBI"
txf_ucsc <- convertToTxFeatures(txdb)
txf2 <- convertToTxFeatures(TxDb)

head(txf_ucsc)
head(txf2)

sgf_ucsc <- convertToSGFeatures(txf2)

sgfc_ucsc <- analyzeFeatures(sit, features = txf2)
sgvc_ucsc <- analyzeVariants(sgfc_ucsc)
sgvc_ucsc



colData(sgfc_ucsc)
rowRanges(sgfc_ucsc)
head(counts(sgfc_ucsc))
head(FPKM(sgfc_ucsc))

df <- plotFeatures(sgfc_ucsc, geneName = "set-16")
dev.off()
par(mfrow = c(11, 1), mar = c(1, 3, 2, 1))
plotSpliceGraph(rowRanges(sgfc_ucsc), geneName = "set-16", toscale = "none", color_novel = "red")
for (j in 1:10) {
  plotCoverage(sgfc_ucsc[, j], geneName = "set-16", toscale = "none")
}

plotCoverage(geneName="set-16",sgfc_ucsc[,1])


sgfc_ucsc
head(sgf_ucsc)
seqlevels(param)
features

#### testing for diff usage, ce10 definitions
sgv <- rowRanges(sgvc_ucsc)
sgvc <- getSGVariantCounts(sgv, sample_info = sit)
sgvc

#### de novo prediction of splicing alterations
sgfc_pred <- analyzeFeatures(sit)
head(rowRanges(sgfc_pred))


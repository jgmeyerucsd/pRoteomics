
sortBam("D:/RNAseq/Bam_Data/2.G4_2.bam","D:/RNAseq/Bam_Data/2.G4_2.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/3.G4_3.bam","D:/RNAseq/Bam_Data/3.G4_3.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/4.G4_4.bam","D:/RNAseq/Bam_Data/4.G4_4.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/5.G4_5.bam","D:/RNAseq/Bam_Data/5.G4_5.sort",maxMemory=6000)

indexBam("D:/RNAseq/Bam_Data/2.G4_2.sort.bam","D:/RNAseq/Bam_Data/2.G4_2.indexed")
indexBam("D:/RNAseq/Bam_Data/3.G4_3.sort.bam","D:/RNAseq/Bam_Data/3.G4_3.indexed")
indexBam("D:/RNAseq/Bam_Data/4.G4_4.sort.bam","D:/RNAseq/Bam_Data/4.G4_4.indexed")
indexBam("D:/RNAseq/Bam_Data/5.G4_5.sort.bam","D:/RNAseq/Bam_Data/5.G4_5.indexed")


sortBam("D:/RNAseq/Bam_Data/N2_1.bam","D:/RNAseq/Bam_Data/N2_1.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/N2_2.bam","D:/RNAseq/Bam_Data/N2_2.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/N2_3.bam","D:/RNAseq/Bam_Data/N2_3.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/N2_4.bam","D:/RNAseq/Bam_Data/N2_4.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/N2_5.bam","D:/RNAseq/Bam_Data/N2_5.sort",maxMemory=6000)

indexBam("D:/RNAseq/Bam_Data/N2_1.sort.bam","D:/RNAseq/Bam_Data/N2_1.indexed")
indexBam("D:/RNAseq/Bam_Data/N2_2.sort.bam","D:/RNAseq/Bam_Data/N2_2.indexed")
indexBam("D:/RNAseq/Bam_Data/N2_3.sort.bam","D:/RNAseq/Bam_Data/N2_3.indexed")
indexBam("D:/RNAseq/Bam_Data/N2_4.sort.bam","D:/RNAseq/Bam_Data/N2_4.indexed")
indexBam("D:/RNAseq/Bam_Data/N2_5.sort.bam","D:/RNAseq/Bam_Data/N2_5.indexed")


sortBam("D:/RNAseq/Bam_Data/3.G4_3.bam","D:/RNAseq/Bam_Data/3.G4_3.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/4.G4_4.bam","D:/RNAseq/Bam_Data/4.G4_4.sort",maxMemory=6000)
sortBam("D:/RNAseq/Bam_Data/5.G4_5.bam","D:/RNAseq/Bam_Data/5.G4_5.sort",maxMemory=6000)

indexBam("D:/RNAseq/Bam_Data/2.G4_2.sort.bam","D:/RNAseq/Bam_Data/2.G4_2.indexed")
indexBam("D:/RNAseq/Bam_Data/3.G4_3.sort.bam","D:/RNAseq/Bam_Data/3.G4_3.indexed")
indexBam("D:/RNAseq/Bam_Data/4.G4_4.sort.bam","D:/RNAseq/Bam_Data/4.G4_4.indexed")
indexBam("D:/RNAseq/Bam_Data/5.G4_5.sort.bam","D:/RNAseq/Bam_Data/5.G4_5.indexed")









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
library(ASpli)
getwd()
setwd("D:/RNAseq/Bam_Data/")
setwd("C:/RNAseq/")

bamFiles <- c( "N2.1.bam", "N2.2.bam", "N2.3.bam","N2.4.bam","N2.5.bam",
               "G4.1.bam","G4.2.bam","G4.3.bam","G4.4.bam","G4.5.bam")
targets <- data.frame( row.names = c("N2_rep1","N2_rep2","N2_rep3","N2_rep4","N2_rep5","G4_rep1", "G4_rep2", "G4_rep3", "G4_rep4", "G4_rep5"),
                       bam = bamFiles,
                       genotype = c("N2","N2", "N2","N2","N2","G4","G4", "G4","G4", "G4") ,
                       stringsAsFactors = FALSE )
bam <- loadBAM(targets)
#### get  features
TxDb <- makeTxDbFromGFF(
  file="D:/RNAseq/Bam_Data/Caenorhabditis_elegans.refGene.ce10.gtf",
  format="gtf")
features <- binGenome( TxDb )




counts <- readCounts (
  features,
  bam,
  l=51)

?readCounts




list.files()
bamFiles <- c( "N2_1.sort.bam", "N2_2.sort.bam", "N2_3.sort.bam","N2_4.sort.bam","N2_5.sort.bam",
               "1.G4_1.sort.bam","2.G4_2.sort.bam","3.G4_3.sort.bam","4.G4_4.sort.bam","5.G4_5.sort.bam")
targets <- data.frame( row.names = c("N2_rep1","N2_rep2","N2_rep3","N2_rep4","N2_rep5","G4_rep1", "G4_rep2", "G4_rep3", "G4_rep4", "G4_rep5"),
                         bam = bamFiles,
                         genotype = c("N2","N2", "N2","N2","N2","G4","G4", "G4","G4", "G4") ,
                         stringsAsFactors = FALSE )
si<-data.frame(sample_name=as.character(row.names(targets)),file_bam=as.character(targets[,1]),stringsAsFactors = F)
bam <- loadBAM(targets)
?readCounts

##### this function gives an ERROR!!!
counts <- readCounts (
  features,
  bam,
  cores = 1,
  l = 51,
  maxISize = 50000,
  minAnchor = NULL )

###### event poiinter

source("http://www.bioconductor.org/biocLite.R")
biocLite("EventPointer")

source("https://bioconductor.org/biocLite.R")
biocLite("spliceR")
si

library(SGSeq)
getwd()
?getBamInfo

setwd("D:/RNAseq/Bam_Data/")





importTranscripts()

#### remapped reads using hisat2
setwd("C:/RNAseq/")
si<-data.frame(sample_name=c("G1","G2","G3","G4","G5","N1","N2","N3","N4","N5"),file_bam=list.files(pattern="bam"),stringsAsFactors = F)

sit<-getBamInfo(sample_info=si,yieldSize=NULL)

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


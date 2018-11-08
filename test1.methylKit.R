rm(list=ls())
# Set packages
devtools::source_gist("4839e615e2401d73fe51")
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("genomation", version = "3.8")
library(genomation)
update.packages("methylKit")
library(methylKit)

#my.File<-list.files("/Users/wwang/Desktop/bisout/bismarkCOV", pattern = ".cov")
#my.File

# Load data
my.test <- "/Users/wwang/Desktop/bisout/bismarkCOV/Recipient_merge.cov"
my.control <- "/Users/wwang/Desktop/bisout/bismarkCOV/Donor_merge.cov"
my.file=list(my.test,my.control)

myobj = readBismarkCoverage(my.file,sample.id=list("CASE1","ctrl1"), assembly="hg38",treatment=c(1,0))
# Since we read the methylation data now, we can check the basic stats about the methylation data such as coverage and percent methylation. We now have a methylRawList object which contains methylation information per sample. The following command prints out percent methylation statistics
getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
meth=unite(myobj, destrand=FALSE)
head(meth)
getCorrelation(meth,plot=TRUE)
PCASamples(meth, screeplot=TRUE)

## Tiling windows analysis
tiles1000=tileMethylCounts(myobj,win.size=1000,step.size=1000, mc.cores=4)
tiles200=tileMethylCounts(myobj,win.size=200,step.size=50, mc.cores=4)
tiles100s50=tileMethylCounts(myobj,win.size=100,step.size=50, mc.cores=4)
head(tiles[[1]],3)

myTileDiff=calculateDiffMeth(tiles, mc.cores=4)

#Finding differentially methylated bases or regions
# using Fisher’s Exact test and logistic regression
myDiff=calculateDiffMeth(meth,mc.cores=2)
head(myDiff)
myDiff10p=getMethylDiff(myDiff,difference=10,qvalue=0.01)
myDiff10p005=getMethylDiff(myDiff,difference=10,qvalue=0.05)
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
myDiff25p005=getMethylDiff(myDiff,difference=25,qvalue=0.05)
head(myDiff25p)
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
write.csv(myDiff10p005, "/Users/wwang/Desktop/bisout/bismarkCOV/MDS_diff10p005.txt")

#  visualize the distribution of hypo/hyper-methylated bases/regions per chromosome using the following function
diffMethPerChr(myDiff,plot=FALSE, qvalue.cutoff=0.01, meth.cutoff=25)

# read the gene BED file
refAnn<- "/Users/wwang/Desktop/gencode.v24.annotation.sorted.bed"
gene.obj=readTranscriptFeatures(refAnn)
# annotate differentially methylated CpGs with
# promoter/exon/intron using annotation data
annotateWithGeneParts(as(myDiff25p005,"GRanges"),gene.obj)
# read the shores and flanking regions and name the flanks as shores
# and CpG islands as CpGi
cpg.obj= readFeatureFlank(refAnn, feature.flank.name = c("CpGi", "shores"))
diffCpGann=annotateWithFeatureFlank(as(myDiff25p005,"GRanges"),cpg.obj$CpGi,cpg.obj$shores, feature.name = "CpGi", flank.name = "shores")
# Summarize methylation information over a set of defined regions such as promoters or CpG islands.
promoters=regionCounts(myobj,gene.obj$promoters)
head(promoters[[1]])
# Get the annotation of differentially methylated regions, we can get the distance to TSS and nearest gene name using the getAssociationWithTSS function from genomation package
diffAnn=annotateWithGeneParts(as(myDiff25p005,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnn))
# It is also desirable to get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
getTargetAnnotationStats(diffAnn, percentage=TRUE, precedence=TRUE)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")


class(meth)
as(meth,"GRanges")

class(myDiff)
as(myDiff,"GRanges")

class(myobj[[1]])
as(myobj[[1]],"methylRaw")

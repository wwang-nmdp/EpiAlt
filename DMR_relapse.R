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
##################################################################################################################################

#my.File<-list.files("/Users/wwang/Desktop/bisout/bismarkCOV", pattern = ".cov")
#my.File

# Load data
my.control <- "/Users/wwang/Desktop/MDS/bisout/bismarkCOV/bedGraph/corrected_merge/recipient_relapse_0_merge.cov"
my.test  <- "/Users/wwang/Desktop/MDS/bisout/bismarkCOV/bedGraph/corrected_merge/recipient_relapse_1_merge.cov"
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
myobj
# DonorTiles100s50=tileMethylCounts(myobj[[2]], win.size=100,step.size=50, mc.cores=4)
# RecipientTiles100s50=tileMethylCounts(myobj[[1]], win.size=100,step.size=50, mc.cores=4)
# tile.file = list(DonorTiles100s50,RecipientTiles100s50)
# tileobj = readBismarkCoverage(tile.file,sample.id=list("CASE1","ctrl1"), assembly="hg38",treatment=c(1,0))

tiles100s20=tileMethylCounts(myobj,win.size=100,step.size=20, mc.cores=4)
tiles30s20=tileMethylCounts(myobj,win.size=30,step.size=20, mc.cores=4)

##################################################################################################################################
# Calculate the DMR 100 frame with 20 step move forward
tile.meth=unite(tiles100s20, destrand=FALSE)
tile.meth2=unite(tiles30s20, destrand=FALSE)
head(tile.meth)

tile100s20Diff=calculateDiffMeth(tile.meth, mc.core = 4)
tile30s20Diff=calculateDiffMeth(tile.meth2, mc.core = 4)

tile100s20Diff50p005=getMethylDiff(tile100s20Diff,difference=50,qvalue=0.05)
tile30s20Diff50q005=getMethylDiff(tile30s20Diff,difference=50,qvalue=0.05)
tile30s20Diff90q005=getMethylDiff(tile30s20Diff,difference=90,qvalue=0.05)


head(tile100s20Diff50p005)
write.table(tile100s20Diff50p005,"/Users/wwang/Desktop/MDS/bisout/bismarkCOV/allDMR_tile100s20Diff50p005.txt")
write.table(tile30s20Diff50q005,"/Users/wwang/Desktop/MDS/bisout/bismarkCOV/allDMR_tile30s20Diff50q005.txt")
write.table(tile30s20Diff90q005,"/Users/wwang/Desktop/MDS/bisout/bismarkCOV/allDMR_tile30s20Diff90q005.txt")

annotateWithGeneParts(as(tile30s20Diff50q005,"GRanges"),gene.obj)
annotateWithGeneParts(as(tile30s20Diff90q005,"GRanges"),gene.obj)

tile30s20Diff50q005Promoter=regionCounts(tile30s20, gene.obj$promoters)

diffMethPerChr(tile100s50Diff50p005,plot=FALSE, qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr(tile30s20Diff90q005,plot=FALSE, qvalue.cutoff=0.01, meth.cutoff=25)



diffAnn_tile100s50Diff50p005=annotateWithGeneParts(as(tile30s20Diff50q005,"GRanges"),gene.obj)
diffAnn_tile100s50Diff90p005=annotateWithGeneParts(as(tile30s20Diff90q005,"GRanges"),gene.obj)

Ann_tile100s50Diff50p005<-getAssociationWithTSS(diffAnn_tile100s50Diff50p005)
Ann_tile100s50Diff90p005<-getAssociationWithTSS(diffAnn_tile100s50Diff90p005)
write.table(Ann_tile100s50Diff50p005,"/Users/wwang/Desktop/MDS/bisout/bismarkCOV/allDMR_tile30s20Diff50q005_feature.txt")
write.table(Ann_tile100s50Diff90p005,"/Users/wwang/Desktop/MDS/bisout/bismarkCOV/allDMR_Ann_tile100s50Diff90p005_feature.txt")



# Calculate the DMR based on promoter regions only
tile100s50.promoter=unite(tile100s50Diff50p005Promoter,destrand=FALSE)
tile100s50.promoterDiff=calculateDiffMeth(tile100s50.promoter, mc.core = 4)
tile100s50.promoterDiff25p005=getMethylDiff(tile100s50.promoterDiff,difference=25,qvalue=0.01)
write.table(tile100s50.promoterDiff25p005,"/Users/wwang/Desktop/bisout/bismarkCOV/DMR_tile100s50.promoterDiff25p005.txt")
tile100s50.promoterDiff25p005Ann = annotateWithGeneParts(as(tile100s50.promoterDiff25p005,"GRanges"), gene.obj)
tile100s50.promoterDiff25p005AnnFeature <- getAssociationWithTSS(tile100s50.promoterDiff25p005Ann)
write.table(tile100s50.promoterDiff25p005AnnFeature,"/Users/wwang/Desktop/bisout/bismarkCOV/DMR_tile100s50.promoterDiff25p005.features.txt")

# the result shows 62.12% DMR are happened in promoter regions
tile.promoters=regionCounts(tiles100s50, gene.obj$promoters)

head(tile.promoters[[1]])
tile100s50.promoter.ann = annotateWithGeneParts(as(tile.promoters[[1]],"GRanges"),gene.obj)
tile100s50Ann = annotateWithGeneParts(as(til100s50Diff,"GRanges"), gene.obj)
head(tile100s50Ann)
write.table(tile.promoters[[1]],"/Users/wwang/Desktop/bisout/bismarkCOV/allDMR_promoters_tileDiff100s50p005.txt")
##################################################################################################################################
# Filter the output
tileDiff50p005 <- getMethylDiff(til100s50Diff,difference=50,qvalue=0.05)
write.table(tileDiff50p005,"/Users/wwang/Desktop/bisout/bismarkCOV/allDMR_donor-recipient_tileDiff50p005.txt" )

##################################################################################################################################
#Finding differentially methylated bases or regions
# using Fisher’s Exact test and logistic regression
myDiff=calculateDiffMeth(meth,mc.cores=2)
head(myDiff)
allDMB <- myDiff
write.table(allDMB, "/Users/wwang/Desktop/bisout/bismarkCOV/allDMB_donor-recipient.txt")
myDiff10p=getMethylDiff(myDiff,difference=10,qvalue=0.01)
myDiff10p005=getMethylDiff(myDiff,difference=10,qvalue=0.05)
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
myDiff25p=getMethylDiff(myDiff,difference=10,qvalue=0.1)
myDiff25p005=getMethylDiff(myDiff,difference=25,qvalue=0.05)
# creat a dataset with broadly info of mehty status for circos plots
myDiff10p01=getMethylDiff(myDiff,difference=40,qvalue=0.49)
write.table(myDiff10p01, "/Users/wwang/Desktop/bisout/bismarkCOV/allDMB_myDiff5p04.txt")

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
diff25p005promoter <- promoters[[1]]
write.table(diff25p005promoter, "/Users/wwang/Desktop/bisout/bismarkCOV/MDS_diff10p005promoter.txt" )
# Get the annotation of differentially methylated regions, we can get the distance to TSS and nearest gene name using the getAssociationWithTSS function from genomation package
diffAnn=annotateWithGeneParts(as(myDiff25p005,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnn))
diff25p005feature <- getAssociationWithTSS(diffAnn)
write.table(diff25p005feature, "/Users/wwang/Desktop/bisout/bismarkCOV/MDS_diff10p005feature.txt")
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

?annotateWithGeneParts


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
diff25p005promoter <- promoters[[1]]
write.table(diff25p005promoter, "/Users/wwang/Desktop/bisout/bismarkCOV/MDS_diff10p005promoter.txt" )
# Get the annotation of differentially methylated regions, we can get the distance to TSS and nearest gene name using the getAssociationWithTSS function from genomation package
diffAnn=annotateWithGeneParts(as(myDiff25p005,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnn))
diff25p005feature <- getAssociationWithTSS(diffAnn)
write.table(diff25p005feature, "/Users/wwang/Desktop/bisout/bismarkCOV/MDS_diff10p005feature.txt")
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

?annotateWithGeneParts


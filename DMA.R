# Clear workspace
rm(list=ls())

# Install bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite()

# Install edageR: 
# Empirical Analysis of Digital Gene Expression Data in R
# biocLite("edgeR")

library("edgeR")

#The downloaded source packages are in
#	‘/tmp/RtmpMF9KNl/downloaded_packages’
# Tiny test for modeling:
# Test NB (Negative Binomial) linear model

#create a dataset
#counts <- matrix(c(2,12,11,0),1,4)
#dimnames(counts) <- list("Locus", c("A_Me","A_Un","B_Me","B_Un"))

# Check dataset
#counts

#design <- cbind(Sample1 = c(1,1,0,0),Sample2 = c(0,0,1,1),A_MvsU = c(1,0,1,0),BvsA_MvsU = c(0,0,1,0))
#design

#fit <- glmFit(counts, design, lib.size=c(100,100,100,100), dispersion= 0.0247)
# The dispersion parameter 
#   ->controls the degree of biological variability

#lrt <- glmLRT(fit, coef="BvsA_MvsU")
#topTags(lrt)
# In this analysis, the first two coefficients are used to model 
# the total number of reads (methylated or unmethylated) 
# for samples 1 and 2, respectively
# P-value for differential methylation (B vs A)


########################################################
#     Differential methylation analysis at CpG loci    #
########################################################
#experimental design


# Load the zipped coverage bed files produced by Bismark
#mcw283<- read.delim("MCW2018-001283_pe.bismark.cov.gz",header= FALSE)
#head(mcw283)

# The columns in the bed file represent: 
#   V1: chromosome number; 
#   V2: start position of the CpG site;
#   V3: end position of the CpG site;
#   V4: methylation proportion;
#   V5: number of methylated Cs;
#   V6: number of unmethylated Cs.

# Read the data into a list in R

targets <- read.delim("targets.txt", stringsAsFactors=FALSE)
targets
Sample <- targets$Sample
Group <- factor(targets$Population)
Group
fn <-paste0(Sample, ".bismark.cov.gz")
data<- list()

for(i in 1:length(Sample)){
  data[[i]] <- read.delim(file=fn[i], header=FALSE)[,-(3:4)]
  names(data[[i]]) <- c("Chr", "Position", "Meth", "Un")
}

position <- sapply(data, function(x) paste(x[,1], x[,2], sep="-"))
position_all <- unique(unlist(position))
  
counts <- matrix(0L, nrow=length(position_all), ncol=2*length(Sample))
for(i in 1:length(Sample)){
  m <- match(position[[i]], position_all)
  counts[m, c(2*i-1,2*i)] <- as.matrix(data[[i]][, 3:4])
}

#  The counts object is a matrix of integer counts with 12 columns,
#  two for each sample. The odd number of columns contain the numbers of methylated Cs,
#  whereas the even number of columns contain the numbers of unmethylated Cs.
#  The genomic positions are used as the row names of the count matrix.

rownames(counts) <- position_all
Sample2 <- rep(Sample, each=2)
Sample2 <- factor(Sample2)
Meth <- rep(c("Me","Un"), length(Sample))
Meth <- factor(Meth, levels=c("Un","Me"))
colnames(counts) <- paste(Sample2, Meth, sep="-")

# Proceed to the edgeR analysis of the methylation data.
# The edgeR package stores data in a simple list-based data object called a DGEList.
# We first create a DGEList object using the count matrix generated before.
# The information of CpG sites is converted into a data frame and
# stored in the genes component of the DGEList object.
library(edgeR)
options(digits=3)
Chr <- gsub("-.*$", "", position_all)
Position <- gsub("^.*-", "", position_all)
Genes <- data.frame(Chr=Chr, Position=Position)
y <- DGEList(counts, genes=Genes, group=rep(Group,each=2))

# Sum up the read counts of both methylated and unmethylated Cs
# Remove low counts
counts_total <- t(rowsum(t(counts), Sample2))
head(counts_total)

# CpG loci that have very low counts across all the samples 
# shall be removed prior to downstream analysis.
keep <- rowSums(counts_total >= 10) == 6
table(keep)
y <- y[keep,,keep.lib.sizes=FALSE]
  
  #keep.lib.sizes=FALSE causes the library sizes to be recomputed after the filtering

# Normalization
  # set the library sizes for each sample to be the average of the total read counts for the methylated and unmethylated libraries:

TotalReadCount <- colMeans(matrix(y$samples$lib.size, nrow=2, ncol=6))
y$samples$lib.size <- rep(TotalReadCount, each=2)
y$samples

  
# Exploring differences between samples
# For a particular CpG site in one sample, 
# the M-value can be computed by subtracting the log2 
# count-per-million (CPM) of the unmethylated Cs from 
# that of the methylated Cs
# This is equivalent to the calculation of the defined M-values
# as the library sizes are set to be the same for each pair of 
# methylated and unmethylated columns and they cancel each other 
# out in the subtraction. A prior count of 2 is added to the calculation 
# of log2-CPM to avoid undefined values and to reduce the variability of 
# M-values for CpG sites with low counts. The calculation of β-value is 
# straight-forward though a small offset may also be added to the calculation.

Beta <- y$counts[, Meth=="Me"] / counts_total[keep, ]
logCPM <- cpm(y, log=TRUE, prior.count=2)
M <- logCPM[, Meth=="Me"] - logCPM[, Meth=="Un"]
colnames(Beta) <- colnames(M) <- Sample
  
# The outputs Beta and M are numeric matrices with six columns, 
#  each of which contains the β-values or M-values calculated 
#  at each CpG site in one sample.

par(mfrow=c(1,2))
plotMDS(Beta, cex = .9,col=rep(1:3, each=2), main="Beta-values")
plot(Beta,pch = 16, cex = .9)
  #legend(x="topright", legend=Sample,col=rep(1:3, each=2))
plotMDS(M, cex = .9, col=rep(1:3, each=2), main="M-values")


# The MDS plots of the methylation levels of the data set. 
#  Methylation levels are measured in beta values (left) and M-values (right). 
#  Samples are separated by the cell population in the first dimension 
#  in both MDS plots.
  
  
#  Design matrix
  
  design <- model.matrix(~ Sample2 + Meth)
  colnames(design) <- gsub("Sample2","",colnames(design))
  colnames(design) <- gsub("Meth","",colnames(design))
  colnames(design)[1] <- "Int"
  design <- cbind(design,Me2=c(0,0,0,0,1,0,1,0,0,0,0,0),Me3=c(0,0,0,0,0,0,0,0,1,0,1,0))
  design
  
  y <- estimateDisp(y, design=design, trend="none")
  y$common.dispersion
  summary(y$prior.df)

# Test and graph
ftGML <- read.table("/Users/wwang/Desktop/filteredGML.txt", header = TRUE)
boxplot(ftGML$GML,ftGML$DnR, ylab= "GML(Global Methylation Level)", col= 'grey')
scatterplot(ftGML$GML~ftGML$age, col="navy", pch = 20)
barplot(table(ftGML$filterRate))

plot(ftGML$age,ftGML$GML, col="blue", pch = 20)
abline(lm(ftGML$GML~ftGML$age), col="red")
line(lowess(ftGML$age,ftGML$GML))

ftGML$DnR

testMDS <- read.table("~/Test_2018-001228.txt", header = TRUE)
summary.matrix(testMDS)
par(mar=c(2, 2, 2, 2) + 0.1)
barplot(table(testMDS$methylated), xlim=c(0,60), cex.name=0.4, col=rainbow(40))
barplot(table(testMDS$unmethylated), xlim=c(0,60),cex.name=0.4,col=rainbow(40))



MDSgml <- read.table("~/MDS_GML.txt",header = TRUE )
MDSgml
summary(MDSgml)

# Install ChAMP and dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("ChAMP")

# Install methylKit and dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("methylKit")

library(methylKit)
setwd("/home/wwang/RRBS")
mySaveFolder="/home/wwang/RRBS/result/methylKit-result"

my.methRaw=processBismarkAln(location = "MCW2018-001242_pe.sam.sorted.sam",
                              sample.id="MCW2018-001242", assembly="H. sapiens",
                              read.context="CpG",save.folder=mySaveFolder)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
## get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
### get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)


#test MethylKit
devtools::source_gist("4839e615e2401d73fe51")
library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE,
               repos=BiocManager::repositories(),
               dependencies=TRUE)
library(methylKit)
file.list=list(system.file("mdsdata", "/Users/wwang/Desktop/bisout/bismarkCOV/MCW2018-001228_pe.bismark.cov.gz", package = "methylKit"),system.file("mdsdata", "/Users/wwang/Desktop/bisout/bismarkCOV/MCW2018-001242_pe.bismark.cov", package = "methylKit"))
myobj = readBismarkCoverage(file.list,sample.id=list("CASE1","ctrl1"), assembly="hg38",treatment=c(1,0))

myobjDB=methRead(file.list,sample.id = list("test1","ctr1"),assembly = "hg38", treatment = c(1,0),context = "CpG", dbtype = "tabix")
myobjDB2=readBismarkCoverage("/Users/wwang/Desktop/bisout/bismarkCOV/MCW2018-001228_pe.bismark.cov", MCW2018-001228, assembly = "hg38", treatment = c(1,0))


# Accounting for covariates

covariates=data.frame(age=c(30,80,34,30,80,40))
sim.methylBase<-dataSim(replicates=6,sites=1000,
                        treatment=c(rep(1,3),rep(0,3)),
                        covariates=covariates,
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
                        )
my.diffMeth3<-calculateDiffMeth(sim.methylBase,
                                covariates=covariates,
                                overdispersion="MN",test="Chisq",mc.cores=1)
# library(Biobase)
library(rMAT)

pwd<-"./"

####################################################
#See vignette on how to get data files from Web site.
####################################################

setwd("/Users/gottarr/rglab/Papers/rMAT-paper2/")  # Set wroking directory where you are your CEL and bpmap files

bpmapFile<-paste(pwd,"P1_CHIP_A.Anti-Sense.hs.NCBIv35.NR.bpmap",sep="")
arrayFile<-paste(pwd,c("MCF_ER_A1.CEL","MCF_ER_A3.CEL","MCF_ER_A4.CEL","MCF_INP_A1.CEL","MCF_INP_A3.CEL","MCF_INP_A4.CEL"),sep="")
arrayFile<-paste(pwd,c("MCF_ER_A1.CEL"),sep="")
arrayFile<-paste(pwd,c("MCF_ER_A1.CEL","MCF_ER_A3.CEL"),sep="")

# Show the all the different sequences
# ReadBPMAPAllSeqHeader(bpmapFile)

# create a tiling Set from the corresponding data
# This will only grep the sequences with Sc
ScSet<-BPMAPCelParser(bpmapFile, arrayFile, verbose=TRUE,groupName="",seqName="chr21")   

# show the object
show(ScSet)

# summarize its content
summary(ScSet)

# Perform the MAT normalization
ScSetNorm<-NormalizeProbes(ScSet, method="MAT",robust=TRUE, all=FALSE, standard=TRUE, verbose=TRUE)

yN<-exprs(ScSetNorm)
y<-exprs(ScSet)

diag(cor(yN,y))

ScSetNorm<-NormalizeProbes(ScSet, method="MAT",robust=FALSE, all=FALSE, standard=TRUE, verbose=TRUE)

yN<-exprs(ScSetNorm)
y<-exprs(ScSet)

diag(cor(yN,y))

ScSetNorm<-NormalizeProbes(ScSet, method="PairBinned",robust=TRUE, all=FALSE, standard=TRUE, verbose=TRUE)

yN<-exprs(ScSetNorm)
y<-exprs(ScSet)

diag(cor(yN,y))

plot(y[,1],pch=".")

# show the object
show(ScSetNorm)

# summarize its content
summary(ScSetNorm)

# Compute the MAT scores and get the enriched regions
# Note that we only need to specify a unique name identifier for the control, here it is SUP
# So it's always a good idea to use unique identifiers for control and treatment
#ScScore<-MATScore(ScSetNorm, cName=NULL, dMax=600, nProbesMin=8, dMerge=300, method="score", threshold=5, verbose=TRUE, bedName="MyBedFile")

RD<-computeMATScore(ScSetNorm,cName="MCF_INP", dMax=600, verbose=TRUE) 

Enrich<-callEnrichedRegions(RD,dMax=600, dMerge=300, nProbesMin=8, method="pValue", threshold=1e-5, verbose=FALSE)  

# Show the object
show(ScScore)

# summarize its content
summary(ScScore)


# You have 2 possibilities to visualize data 

############################################
# 1) rtracklayer package
##########################################
library(rtracklayer)

#Export to BED format :
export(Enrich,"EnrichedData.bed")

genome(Enrich)<-"sacCer2"
names(Enrich)<-"chrI"

#Viewing the targets
session<- browserSession("UCSC")
track(session,"target") <- Enrich

#Get the first feature
subEnrich<-Enrich[2,]

#View with GenomeBrowser
view<- browserView(session,range(subEnrich) * -2)


#############################################
# 2) This is using the GenomeGraphs package
##########################################
library(GenomeGraphs)
mart<-useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genomeAxis<-makeGenomeAxis(add53 = TRUE,add35 = TRUE)

# Here I only look at chr1 between 1 and 50000

minbase<-15400000
maxbase<-15600000 

minbase<-14100000
maxbase<-15000000

genesplus<-makeGeneRegion(start = minbase, end = maxbase, strand = "+", chromosome = 21, biomart = mart) 
genesmin<-makeGeneRegion(start = minbase, end = maxbase, strand = "-", chromosome = "21", biomart = mart)

# Here I create a Generic array for chr1 ony
RD1<-RD[space(RD)=="chr21",]
Enrich1<-Enrich[space(Enrich)=="chr21",]
MatScore<-makeGenericArray(intensity=as.matrix(score(RD1)),  probeStart=start(RD1), dp=DisplayPars(size=1, color="black", type="l"))

# Here I create overlays to look at the enriched regions (above>threshold)
rectList<- makeRectangleOverlay(start = start(Enrich1), end = end(Enrich1), region = c(1, 4), dp = DisplayPars(color = "green", alpha = 0.1))

gdPlot(list("score" = MatScore, "Gene +" = genesplus, Position = genomeAxis, "Gene -" = genesmin), minBase = minbase, maxBase = maxbase, labelCex = 1,overlays=rectList) 


gdPlot(list("score" = MatScore,  Position = genomeAxis), minBase = minbase, maxBase = maxbase, labelCex = 1) 



### My test

# library(Biobase)
library(rMAT)

bpmapFile<-"P1_CHIP_A.Anti-Sense.hs.NCBIv35.NR.bpmap"
arrayFile<-c("MCF_ER_A1.CEL","MCF_ER_A3.CEL","MCF_ER_A4.CEL","MCF_INP_A1.CEL","MCF_INP_A3.CEL","MCF_INP_A4.CEL")

# Show the all the different sequences
ReadBPMAPAllSeqHeader(bpmapFile)

# create a tiling Set from the corresponding data
# This will only grep the sequences with Sc
ERset<-BPMAPCelParser(bpmapFile, arrayFile, seqName="chr21")

# Perform the MAT normalization
# ERsetNorm1<-NormalizeProbes(ERset, method="MAT",robust=TRUE)
ERsetNorm<-NormalizeProbes(ERset, method="PairBinned",robust=TRUE)

# Compute MAT scores
ERscore<-computeMATScore(ERsetNorm,cName="INP") 
ERscore

# Export a wig file
export(ERscore,con="ERscore.wig")



library(rtracklayer)

genome(Enrich)<-"sacCer2"
names(Enrich)<-"chrI"

#Viewing the targets
session<- browserSession("UCSC")
track(session,"target") <- Enrich

#Get the first feature
subEnrich<-Enrich[2,]

#View with GenomeBrowser
view<- browserView(session,range(subEnrich) * -2)



library(GenomeGraphs)
mart<-useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
genomeAxis<-makeGenomeAxis(add53 = TRUE,add35 = TRUE)

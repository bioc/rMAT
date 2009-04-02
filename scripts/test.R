# library(Biobase)
library(rMAT)
pwd<-"/rMAT/inst/doc/"


####################################################
#See vignette on how to get data files from Web site.
####################################################

setwd("/rMAT/inst/doc/")  # Set wroking directory where you are your CEL and bpmap files

bpmapFile<-paste(pwd,"Sc03b_MR_v04_10000.bpmap",sep="")
arrayFile<-paste(pwd,c("Swr1WTIP_Short.CEL"),sep="")

# Show the all the different sequences
ReadBPMAPAllSeqHeader(bpmapFile)

# create a tiling Set from the corresponding data
# This will only grep the sequences with Sc
ScSet<-BPMAPCelParser(bpmapFile, arrayFile, verbose=TRUE,groupName="Sc")     

# show the object
show(ScSet)

# summarize its content
summary(ScSet)

# Perform the MAT normalization
ScSetNorm<-NormalizeProbes(ScSet, method="MAT",robust=FALSE, all=FALSE, standard=TRUE, verbose=TRUE)

# show the object
show(ScSetNorm)

# summarize its content
summary(ScSetNorm)

# Compute the MAT scores and get the enriched regions
# Note that we only need to specify a unique name identifier for the control, here it is SUP
# So it's always a good idea to use unique identifiers for control and treatment
ScScore<-MATScore(ScSetNorm, cName=NULL, dMax=600, nProbesMin=8, dMerge=300, method="score", threshold=5, verbose=TRUE, bedName="MyBedFile")

# Show the object
show(ScScore)

# summarize its content
summary(ScScore)

# This is using the GenomeGraphs package

library(GenomeGraphs)
mart<-useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
genomeAxis<-makeGenomeAxis(add53 = TRUE,add35 = TRUE)

# Here I only look at chr1 between 1 and 200000

minbase<-1 
maxbase<-100000 
genesplus<-makeGeneRegion(start = minbase, end = maxbase, strand = "+", chromosome = "I", biomart = mart) 
genesmin<-makeGeneRegion(start = minbase, end = maxbase, strand = "-", chromosome = "I", biomart = mart)

# Here I create a Generic array for chr1 ony

MatScore<-makeGenericArray(intensity=as.matrix(ScScore@score[ScScore@featureChromosome=="chr1"]),  probeStart=ScScore@featurePosition[ScScore@featureChromosome=="chr1"], dp=DisplayPars(size=1, color="black", type="l"))

# Here I create overlays to look at the enriched regions (above>threshold)

featurePositionForRegion<-ScScore@featurePosition[ScScore@featureChromosome=="chr1" & ScScore@featurePosition< maxbase & ScScore@featurePosition> minbase]
regIndexForRegion<-ScScore@regIndex[ScScore@featureChromosome=="chr1" & ScScore@featurePosition< maxbase & ScScore@featurePosition> minbase]
RegionUnique<-unique(regIndexForRegion[regIndexForRegion>0])
rectList<-vector("list",length(RegionUnique))
for(i in 1:length(RegionUnique))
{
  ## Minimum
  m<-min(featurePositionForRegion[regIndexForRegion==RegionUnique[i]])
  ## Maximum
  M<-max(featurePositionForRegion[regIndexForRegion==RegionUnique[i]])
  rectList[i] <- makeRectangleOverlay(start = m, end = M, region = c(1, 4), dp = DisplayPars(color = "green", alpha = 0.1))
}

gdPlot(list("score" = MatScore, "Gene +" = genesplus, Position = genomeAxis, "Gene -" = genesmin), minBase = minbase, maxBase = maxbase, labelCex = 1, overlays=rectList) 




#getBM(attributes=c("chromosome_name","uniprot_swissprot","start_position","sequence_gene_chrom_start","end_position"), filters=c("chromosome_name","start","end"), values=list (Nom_chromosome,50679523,51195524), mart=ensembl)


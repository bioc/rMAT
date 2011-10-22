\name{computeMATScore}
\alias{computeMATScore}



\title{Detection of enriched regions}

\description{
   This function is used to compute the rMAT scores following normalization of expression values in order to locate putative enriched regions. 
  This function is now defunct now defunct and you should instead use 'computeMATScore'.
}

 \usage{
 computeMATScore(tilingSet, cName=NULL, dMax=600, verbose=FALSE)
        
 }
 
 
  \arguments{
     \item{tilingSet}{This object contains an ExpressionSet}
   \item{cName}{Unique identifier of control name}
   \item{dMax}{An integer value. The sliding window side of which the adjacent probes are to average upon in order to compute the rMAT score.}
   \item{verbose}{A logical value. If verbose is TRUE, progress information would be displayed.}

 }

 \details{
   For more details on the calculation of the rMAT score, pvalues, etc, please refer to the following paper:
   Johnson et al. Model-based analysis of tiling-arrays for ChIP-chip. Proc Natl Acad Sci USA (2006) vol. 103 (33) pp. 12457-62
 }
 
 \value{
   The rMAT Score, pValues, and regions. For the regions vector, let 0 denotes the unenriched region. If an enriched region is found, the interval of the region is labeled by a none 0 value. The first region detected is abeled 1 and the next regions are subsequently incremented. 
 }
 
 \author{
 Charles Cheung, \email{cykc@interchange.ubc.ca} and
 Raphael Gottardo, \email{rgottard@fhcrc.org}
 Arnaud Droit, \email{arnaud.droit@crchuq.ulaval.ca}
 }


\seealso{
 \code{NormalizeProbes}, \code{computeMATScore}, \code{callEnrichedRegions} for normalizing expression values before computing the rMAT enriched regions.
}
\examples{

####################################################
#The data are in inst/doc folder in rMAT package.
####################################################


pwd<-"" #INPUT FILES- BPMAP, ARRAYS, etc.
path<- system.file("doc/Sc03b_MR_v04_10000.bpmap",package="rMAT")

bpmapFile<-paste(pwd,path,sep="")

pathCEL<- system.file("doc/Swr1WTIP_Short.CEL",package="rMAT")
arrayFile<-paste(pwd,c(pathCEL),sep="")


# Show the all the different sequences
ReadBPMAPAllSeqHeader(bpmapFile)

# create a tiling Set from the corresponding data
# This will only grep the sequences with Sc
ScSet<-BPMAPCelParser(bpmapFile, arrayFile, verbose=FALSE,groupName="Sc")     

# show the object
show(ScSet)

# summarize its content
summary(ScSet)

ScSetNorm<-NormalizeProbes(ScSet, method="MAT",robust=FALSE, all=FALSE, standard=TRUE, verbose=FALSE)
	
RD<-computeMATScore(ScSetNorm,cName=NULL, dMax=600, verbose=TRUE) 
	
}


\keyword{file}
\keyword{IO}

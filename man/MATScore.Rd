\name{MATScore}
\alias{MATScore}



\title{Detection of enriched regions}

\description{
   This function is used to compute the rMAT scores following normalization of expression values in order to locate putative enriched regions. 
  This function is now defunct now defunct and you should instead use 'computeMATScore'.
}

 \usage{
 MATScore(tilingSet, cName="NULL", dMax=600,nProbesMin=8, dMerge=300,method="score",threshold=5,verbose=FALSE,bedName="NULL")
        
 }
 
 
 \arguments{
   \item{tilingSet}{This object contains an ExpressionSet}
   \item{cName}{Unique identifier of control name}
   \item{dMax}{An integer value. The sliding window side of which the adjacent probes are to average upon in order to compute the rMAT score.}
   \item{nProbesMin}{An integer value. The minimum number of probes to average upon. If the number of probes within the interval is less than nProbesMin, the rMAT score of the region will not be computed.}
   \item{dMerge}{An integer value. The maximum size to merge adjacent probes and categorize them as one region for scores of  adjacent probes uniformly above the input threshold.}
   \item{method}{A character string value equal to "score", "pValue" or "FDR". "score" denotes the method of calling enriched regions based sliding widow scores.
 "pValue" denotes the method of calling enriched regions based on p-values. Method "FDR" uses an FDR procedure to call regions. See Details below.}
   \item{threshold}{An integer value. The threshold of rMAT Score to be labeled as an enriched region. For method=1 or 3, the higher the score, the more confident we are about enriched regions. For method=2, the lower the score, the more confident we are about enriched regions.}
   \item{verbose}{A logical value. If verbose is TRUE, progress information would be displayed.}
   \item{bedName}{This file file includes columns "chromosome rMATScore region pValue" for each probe.}
 }
 
 \details{
   For more details on the calculation of the rMAT score, pvalues, etc, please refer to the following paper:
   Johnson et al. Model-based analysis of tiling-arrays for ChIP-chip. Proc Natl Acad Sci USA (2006) vol. 103 (33) pp. 12457-62
 }
 
 \value{
   The rMAT Score, pValues, and regions. For the regions vector, let 0 denotes the unenriched region. If an enriched region is found, the interval of the region is labeled by a none 0 value. The first region detected is labeled 1 and the next regions are subsequently incremented. 
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


#pwd<-"" #INPUT FILES- BPMAP, ARRAYS, etc.
#path<- system.file("extdata", "Sc03b_MR_v04_10000.bpmap",package="rMAT")

#bpmapFile<-paste(pwd,path,sep="")

#pathCEL<- system.file("extdata", "Swr1WTIP_Short.CEL",package="rMAT")
#arrayFile<-paste(pwd,c(pathCEL),sep="")


# Show the all the different sequences
#ReadBPMAPAllSeqHeader(bpmapFile)

# create a tiling Set from the corresponding data
# This will only grep the sequences with Sc
#ScSet<-BPMAPCelParser(bpmapFile, arrayFile, verbose=FALSE,groupName="Sc")     

# show the object
#show(ScSet)

# summarize its content
#summary(ScSet)


#ScSetNorm<-NormalizeProbes(ScSet, method="MAT",robust=FALSE, all=FALSE, standard=TRUE, verbose=FALSE)
	
#ScScore<- MATScore(ScSetNorm, cName=NULL, dMax=600,nProbesMin=8, dMerge=300,method="score",threshold=5,verbos=TRUE,bedName="MyBedFile")	
	
}

\keyword{file}
\keyword{IO}

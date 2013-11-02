\name{BPMAPCelParser}
\alias{BPMAPCelParser}

\title{BPMAP and CEL files Reader}

\description{
  One-step reading of BPMAP and CEL files, using Fusion SDK and affxparser.
}

\usage{
BPMAPCelParser(BPMAPFileName, CelFileNames, genomeName=NULL, verbose=FALSE,groupName="",seqName="")
}


\arguments{
  \item{BPMAPFileName}{String containing the full filename of the BPMAP file.}
  \item{CelFileNames}{Vector of strings containing full filenames of CEL files. i.e. c("F1.CEL", "F2.CEL"}
  \item{genomeName}{ String containing the genome name used.}
  \item{groupName}{ String containing the group of genome name used.}
  \item{seqName}{ String containing the group of sequence name (e.g. chromosome) used.}
  \item{verbose}{If verbose is selected, the progress and additional information will be displayed while the function is running}
}

\details{
  This function returns an object of class \code{tilingSet} containing all necessary information: probe sequences, genomic positions, chromosomes as well as the probe intensities. 
}

\value{
An object of class \code{tilingSet}.
}

\author{
  Charles Cheung, \email{cykc@interchange.ubc.ca} and
  Raphael Gottardo, \email{rgottard@fhcrc.org}
  Arnaud Droit, \email{arnaud.droit@crchuq.ulaval.ca}
}


\seealso{
  \code{affyTile} for information about the package.
}

\examples{

####################################################
#The data are in inst/doc folder in rMAT package
####################################################

pwd<-"" #INPUT FILES- BPMAP, ARRAYS, etc.
path<- system.file("extdata", "Sc03b_MR_v04_10000.bpmap",package="rMAT")

bpmapFile<-paste(pwd,path,sep="")

pathCEL<- system.file("extdata", "Swr1WTIP_Short.CEL",package="rMAT")
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
}

\keyword{file}
\keyword{IO}

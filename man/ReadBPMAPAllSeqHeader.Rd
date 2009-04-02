\name{ReadBPMAPAllSeqHeader}
\alias{ReadBPMAPAllSeqHeader}

\title{Reading All the BPMAP Sequence Header}

\description{
  Reading the header of a specified sequence in the BPMAP file. Several sequences could be stored in a single Affymetrix Tiling Array. For example, an array could contain probes from Chromosome 21 and Chromosome 22. The sequenceNum uniquely specifies a sequence. Information about this sequence could be determined in this function. The total number of sequences a tiling array contains can be determined in ReadBPMAPHeader(fileName). The sequenceNum indexes from 0 to (total number of sequences -1).
  
}

\usage{
ReadBPMAPAllSeqHeader(fileName)
}

\arguments{
  \item{fileName}{the full path of the BPMAP file to be read.}
}

\details{
  The BPMAP Sequence Header gives information about the design of the tiling array.
}

\value{
  A list of vectors containing SeqName, GroupName, version, npnrobeMapping, seqNum, and NumHits.
}

\author{
  Charles Cheung, \email{cykc@interchange.ubc.ca} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Arnaud Droit, \email{arnaud.droit@ircm.qc.ca}
}



\seealso{
  \code{BPMAPCelParser()} for an one-step BPMAP/CEL parser and \code{affyTile} for information about the package.
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

}

\keyword{file}
\keyword{IO}

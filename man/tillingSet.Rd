\name{tilingSet}
\alias{tilingSet}
\alias{tilingSet-class}
\alias{rbind,tilingSet-method}
\alias{\%in\%,ANY,tilingSet-method}


\title{This object contains an ExpressionSet } 

\description{
This object contains an ExpressionSet and has the following additional slots: 
genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber
}

\details{
Tiling set objects can also be combined using the \code{rbind} methods. This is particularly useful when several arrays span a genome/chromosome.
}

\usage{
 tilingSet('tilingSet', featureChromosome,featurePosition, featureCopyNumber,exprs, genomeName, featureSequence, experimentData)
}

\arguments{
\item{tilingSet}{This object contains an ExpressionSet}
\item{genomeName}{ String containing the genome name used (vector).}
\item{featureChromosome}{ Factor containing the name of chromosome used (vector).}
\item{featurePosition}{ String containing the Position of the sequences (vector).}
\item{featureCopyNumber}{ String containing the copy number of sequence (vector).}
\item{exprs}{ String containing the expresion data of enriched region (matrix with n column.}
\item{featureSequence}{String containing the sequence (vector).}
\item{experimentData}{String containing the type of experiments.}
}

\author{
Charles Cheung, \email{cykc@interchange.ubc.ca} and
Raphael Gottardo, \email{rgottard@fhcrc.org}
Arnaud Droit, \email{arnaud.droit@crchuq.ulaval.ca}
}

\references{
W. E. Johnson, Li, W., Meyer, C. A., Gottardo, R., Carroll, J. S., Brown, M., and Liu, X. S. (2006). Model-based analysis of tiling-arrays for ChIP-chip. PNAS 103:12457-12462.
}

\seealso{
\code{\link{BPMAPCelParser}}, \code{\link{NormalizeProbes}}
}


\examples{

featureChromosome=factor(c("chr1","chr1","chr1","chr1"))
featurePosition=c(as.integer(47193),as.integer(47197),as.integer(47201),
    as.integer(47205))
featureCopyNumber=c(as.integer(1),as.integer(1),as.integer(1),as.integer(1))
a=5.379897
exprs=matrix(a,nrow=4)
genomeName="Sc03b_MR_v04_10000"
featureSequence=c("TCATCAAGGGAAGAGAGTCTCTCAG","TGATCATCACGGGACTTCTGGTTTA","CGGGACTTCTGGTTTATGGAACTAT","ACTTCTGGTTTATGGAACTATAAGG")

newSet <- new('tilingSet', featureChromosome=featureChromosome,
    featurePosition=featurePosition,featureCopyNumber=featureCopyNumber,
    exprs=exprs, genomeName=genomeName, featureSequence=featureSequence)

}


\keyword{models}

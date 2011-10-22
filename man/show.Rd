\name{show,tilingSet-method}
\docType{methods}
\alias{show,tilingSet-method}
\title{show Method for tiling set object}

\description{
This methods show the content of tilinSet objets}

\usage{
\S4method{show}{tilingSet}(object)
}

\arguments{
\item{object}{Object returned of class \code{tilingSet}}
}

\details{
TilingSet contais an ExpressionSet and has the following additional slots: genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber, featureSequence.
}

\author{
Charles Cheung, \email{cykc@interchange.ubc.ca} and
Raphael Gottardo, \email{rgottard@fhcrc.org}
Arnaud Droit, \email{arnaud.droit@crchuq.ulaval.ca}
}


\seealso{
\code{\link{BPMAPCelParser}}, \code{\link{NormalizeProbes}}
}

\keyword{print}

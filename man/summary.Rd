\name{summary,MAT-method}
\docType{methods}
\alias{summary,MAT-method}
\alias{summary,tilingSet-method}
\title{Summary Method for MAT Object}

\description{
This methods summarize the objects MAT or tilinSet}

\usage{
\S4method{summary}{MAT}(object)
\S4method{summary}{tilingSet}(object)

}

\arguments{
\item{object}{Object returned from \code{\link{MAT}} .}
}

\details{
MAT is an object summarizing the output from the MAT function with the following slots : genomeName,featurePosition,featureChromosome,score,regIndex,method,threshold.

TilingSet is an object containing an ExpressionSet with the follonwing slots: 
genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber


}


\author{
  Charles Cheung, \email{cykc@interchange.ubc.ca} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Arnaud Droit, \email{arnaud.droit@ircm.qc.ca}
}

\seealso{
\code{\link{MAT}}
}

\keyword{print}

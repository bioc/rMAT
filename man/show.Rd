\name{show,MAT-method}
\docType{methods}
\alias{show,MAT-method}
\alias{show,tilingSet-method}
\title{show Method for MAT Object}

\description{
This methods show the objects MAT or tilinSet}

\usage{
\S4method{show}{MAT}(object)
\S4method{show}{tilingSet}(object)

}

\arguments{
\item{object}{Object returned from \code{\link{MAT}} .}
}

\details{
MAT or TilingSet contains an ExpressionSet and has the following additional slots: genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber, featureSequence. This methods shows the information.
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

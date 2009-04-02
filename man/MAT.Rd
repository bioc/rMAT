\name{MAT}
\alias{MAT}
\alias{MAT-class}
\alias{\%in\%,ANY,MAT-method}


\title{This object contains an ExpressionSet } 

\description{
An object summarizing the output from the MAT function with the following slots : genomeName,featurePosition,featureChromosome,score,regIndex,method,threshold}

\usage{

	
NewMAT<-new('MAT', genomeName="Sc", featurePosition, featureChromosome, score, regIndex, method, threshold)
}


\author{
  Charles Cheung, \email{cykc@interchange.ubc.ca} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Arnaud Droit, \email{arnaud.droit@ircm.qc.ca}
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
	
\code{\link{tilingSet}}
}

\examples{
	
#	genomeName="Sc03b_MR_v04"
#	featurePosition="4"
#	featureChromosome="chr1"
#	score="8.8623"
#	regIndex="177"
#	method=""
#	threshold="5"

#NewMAT<-new('MAT', genomeName="Sc", featurePosition, featureChromosome, score, regIndex, method, threshold)
	
}


\keyword{models}
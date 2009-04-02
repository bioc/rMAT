setClass("tilingSet", contains="ExpressionSet", representation(genomeName="character", featureSequence="character", featurePosition="integer", featureChromosome="character", featureCopyNumber="integer"),prototype(genomeName=character(0), featureSequence=character(0), featureChromosome=character(0), featurePosition=integer(0), featureCopyNumber=integer(0),  featureCopyNumber=integer(0)),

	validity=function(object){
		
		if(!((is.character(object@genomeName))) )
		{
		
		stop("The genomeName is not a caracter")	
		}
		
		if(!all(is.character(object@featureSequence)))
		{
		stop ("All featureSequence must be a character")	
		}
		
		if(length(object@featureSequence)!=length(object@featurePosition) | length(object@featureSequence)!=length(object@featureChromosome))		{		stop("All features (Sequence, Position, Chromosome) must have the same length!")		}
		
		if(!all(is.numeric(object@featurePosition)))
		{
		stop ("All featurePosition must be a numeric")	
		}
		
		if((!all(is.character(object@featureChromosome))))
		{
		stop ("All featureChromosome must be a character")	
		}
		
		if((!all(is.numeric(object@featureCopyNumber))))
		{
		stop ("All featureCopyNumber must be a numeric")	
		}

		
	})

setClass("MAT", representation(genomeName="character", featureChromosome="character", featurePosition="integer", score="numeric", pValue="numeric", regIndex="numeric",method="character",threshold="numeric"), prototype(score=rep(numeric(0),0), pValue=rep(numeric(0),0), regIndex=rep(numeric(0),0), method=character(0), threshold=numeric(0),featureChromosome=character(0), featurePosition=integer(0)))


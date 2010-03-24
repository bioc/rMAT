setClass("tilingSet", contains="ExpressionSet", representation(genomeName="character", featureSequence="character", featurePosition="integer", featureChromosome="factor", featureCopyNumber="integer"),prototype(genomeName=character(0), featureSequence=character(0), featureChromosome=factor(0), featurePosition=integer(0), featureCopyNumber=integer(0),  featureCopyNumber=integer(0)),

	validity=function(object){
		
		if(!((is.character(object@genomeName))) )
		{
		
		stop("The genomeName is not a caracter")	
		}
		
		if(!all(is.character(object@featureSequence)))
		{
		stop ("All featureSequence must be a character")	
		}
		
		if(length(object@featureSequence)!=length(object@featurePosition) | length(object@featureSequence)!=length(object@featureChromosome))
		{
		stop("All features (Sequence, Position, Chromosome) must have the same length!")
		}
		
		if(!all(is.numeric(object@featurePosition)))
		{
		stop ("All featurePosition must be a numeric")	
		}
		
		if((!is.factor(object@featureChromosome)))
		{
		stop ("featureChromosome must be a factor")	
		}
		
		if((!all(is.numeric(object@featureCopyNumber))))
		{
		stop ("All featureCopyNumber must be a numeric")	
		}
	})



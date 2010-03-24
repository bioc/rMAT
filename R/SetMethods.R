setMethod("show", "tilingSet",
function(object)
{
    cat("Object of class 'tilingSet'","\n")
    cat("This object contains an ExpressionSet and has the following additional slots: \n")
    cat("genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber, featureSequence\n")
}
)

setMethod("summary", signature("tilingSet"),
    function(object) {
      cat("   Genome interrogated: ",sort(object@genomeName)," \n")
      chr<-levels(object@featureChromosome)
      cat("   Chromosome(s) interrogated: ")
      cat(chr,sep=", ")
      cat(" \n")
      cat("   Sample name(s): ",sampleNames(object@phenoData)," \n")
      cat("   The total number of probes is: ",length(object@featureChromosome)," \n")
      cat("   Preprocessing Information \n")
      cat("     - Transformation:",preproc(object@experimentData)$transformation, "\n")
      cat("     - Normalization:",preproc(object@experimentData)$normalization, "\n")      
    }
)


## Bind several tiling sets together

setMethod("rbind", "tilingSet", function(..., deparse.level=1) {
  args <- list(...)

  names<-lapply(args,function(x){sampleNames(x@phenoData)})
  lnames<-unlist(lapply(names,"length"))
  # Check that the number of names is the same for all
  if(length(unique(lnames))!=1)
  {
    stop("Objects to be concatenated should have the same number of columns")    
  }

  featureChromosome<-unlist(lapply(args,function(x){x@featureChromosome}))
  featurePosition<-unlist(lapply(args,function(x){x@featurePosition}))
  featureCopyNumber<-unlist(lapply(args,function(x){x@featureCopyNumber}))
  featureSequence<-unlist(lapply(args,function(x){x@featureSequence}))
  
  ord<-order(featureChromosome,featurePosition)
  
  featureChromosome<-unlist(lapply(args,function(x){x@featureChromosome}))
  y<-do.call("rbind",lapply(args,function(x){exprs(x)}))
  
  ## For the sample name, I simply use the first set
  colnames(y)<-names[[1]]
  rownames(y)<-NULL
  newSet<-new('tilingSet', featureChromosome=featureChromosome[ord],featurePosition=featurePosition[ord],
  featureCopyNumber=featureCopyNumber[ord], exprs=y[ord,], genomeName=unlist(lapply(args,function(x){x@genomeName})),
  featureSequence=featureSequence[ord], experimentData=args[[1]]@experimentData)
  return(newSet)
})


#now we can only assume there is one chromosome per time
#Furthermore, the vectors must be sorted

computeMATScore<-function(tilingSet, cName=NULL, dMax=600, verbose=FALSE)
{
  if(class(tilingSet)!="tilingSet")
  {
    stop("tilingSet must be an object of class tilingSet (e.g. returned by BPMAPCelParser)!")
  }
  if(preproc(tilingSet@experimentData)$transformation!="log")
  {
    stop("The probe measurements need to be log transformed!")
  }
  if(preproc(tilingSet@experimentData)$normalization=="none")
  {
    warning("You should probably normalize your data before using this function")
  }
  
  #whether the Control array is available
  y<-exprs(tilingSet)
  
  if(!is.null(cName))
  {
    sNames<-sampleNames(tilingSet)       
    if(length(grep(cName,sNames))==0)
    {
      stop("The variable 'cName' must be a subset of the control sample name")
    }
    if(verbose)
    {
      cat("You are using: ",sNames[grep(cName,sNames)], " as the control\n" )
    }
    C<-as.matrix(y[,grep(cName,sNames)])
    I<-as.matrix(y[,-grep(cName,sNames)])
    nArraysC<-ncol(C)
  }
  else
  {
    C<-0
    nArraysC<-0
    I<-as.matrix(y)
  }

  nProbes<-nrow(I)
  nArraysI<-ncol(I)

  seqNum<-as.numeric(as.factor(tilingSet@featureChromosome))  

  obj<-.C("MATScore",
  as.double(t(C)),
  as.double(t(I)),
  as.integer(nProbes),
  as.integer(nArraysC),
  as.integer(nArraysI),
  as.integer(tilingSet@featurePosition),
  as.double(dMax),
  MATScore=double(nProbes),
  as.integer(seqNum),
  package="rMAT")

  if(verbose)
  {
    cat("** Finished processing ",nProbes," probes on ",nArraysC+nArraysI," arrays **\n");
  }

  # Here I assume that all the sequences have the same length
  ranges<-IRanges(as.integer(tilingSet@featurePosition),width=1)
  RD<-RangedData(ranges, score=obj$MATScore, space = tilingSet@featureChromosome)
  # Remove duplicated probes
  ind<-unlist(lapply(RD,function(x)duplicated(start(x))))
  RD[!ind,]
}

callEnrichedRegions<-function(MatScore, dMax=600, dMerge=300, nProbesMin=8, method="score", threshold=5, verbose=FALSE)
{
  if(!is(MatScore,"RangedData"))
  {
    stop("MatScore must be an object of class RangedData (e.g. returned by computeMATScore)!")
  }

  if(method=="score")
  {
    methodNum<-1
  }
  else if(method=="pValue")
  {
    if(threshold>1 | threshold<0)
    {
      stop("When the pValue method is selected, the threshold should be between 0 and 1")
    }
    methodNum<-2
  }
  else if(method=="FDR")
  {
    if(threshold >1 | threshold<0)
    {
      stop("When the FDR method is selected, the threshold should be between 0 and 1")
    }
    methodNum<-3
  }
  else
  {
    stop("Argument 'Method' must be either score or pValue or FDR")
  }
  chrAll<-space(MatScore)
  seqNum<-as.numeric(as.factor(chrAll))
  nProbes<-length(seqNum)
  startMat<-start(MatScore)
  scoreMat<-MatScore$score
  
  obj<-.C("callEnrichedRegions",
  as.double(scoreMat),
  as.integer(nProbes),
  as.integer(startMat),
  as.double(dMerge),
  as.double(dMax),  
  as.double(threshold),
  pValue=double(nProbes),
  as.integer(methodNum),
  regions=integer(nProbes), 
  as.integer(verbose),
  as.integer(seqNum),
  numRegions = integer(1),
  package="rMAT")

  pValue<-obj$pValue


  if(verbose)
  {
    cat("** Number of Enriched regions is ", obj$numRegions, " **\n")
  }
  numRegions<-obj$numRegions
  if(numRegions > 0)
  {
    # We have detected some regions
    obj<-.C("getIndices",
    as.integer(obj$regions),
    as.integer(nProbes),
    as.integer(numRegions),
    Start=integer(numRegions),
    End=integer(numRegions),    
    package="rMAT")
    
    center<-score<-start<-end<-chr<-rep(0,numRegions)
    
    for(i in 1:numRegions)
    {
      start[i]<-startMat[obj$Start[i]]
      end[i]<-startMat[obj$End[i]]
      score[i]<-max(scoreMat[obj$Start[i]:obj$End[i]])
      center[i]<-startMat[(obj$Start[i]:obj$End[i])[which.max(scoreMat[obj$Start[i]:obj$End[i]])]]
      chr[i]<-chrAll[obj$Start[i]]
    }
    ind<-obj$End-obj$Start>nProbesMin

    if(verbose)
    {
      cat("** ", sum(!ind), " enriched regions have less than ", nProbesMin, " probes and are filtered out**\n")
    }
    
    ranges<-IRanges(start=start[ind],end=end[ind])
    chr<-chr[ind]
    score<-score[ind]
    center<-center[ind]
  }
  else
  {
    if(verbose)
    {
      cat("** No regions to output **\n")
    }
    ranges<-IRanges(NULL)
    return(RangedData(ranges))
  }
  # Here I assume that all the sequences have the same length
  RD<-RangedData(ranges, space = chr, score=score, center=center)
  RD
}


MATScore<-function(tilingSet, cName=NULL, dMax=600, nProbesMin=8, dMerge=300, method="score", threshold=5, verbose=FALSE, bedName=NULL)
{
  .Defunct("computeMATScore",package="rMAT","This function is now Defunct. The MAT score function has been slipt into two seperated functions to facilitate export to wig/bed files, see the vignette for more details")
}

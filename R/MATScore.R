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
  RD
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
  print(obj$regions)


  if(verbose)
  {
    cat("** Number of Enriched regions is ", obj$numRegions, " **\n")
  }
  numRegions<-obj$numRegions
  print(numRegions)
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
    
    score<-start<-end<-chr<-rep(0,numRegions)
    print(obj$Start)
    print(obj$End)
    
    for(i in 1:numRegions)
    {
      start[i]<-startMat[obj$Start[i]]
      end[i]<-startMat[obj$End[i]]
      score[i]<-mean(scoreMat[obj$Start[i]:obj$End[i]])
      chr[i]<-chrAll[obj$Start[i]]
    }
    ind<-length(obj$End-obj$Start)>nProbesMin
    ranges<-IRanges(start=as.integer(start[ind]),end=as.integer(end[ind]))
    chr<-chr[ind]
    score<-score[ind]
  }
  else
  {
    if(verbose)
    {
      cat("** No regions to output **\n")
    }
    ranges<-NULL
  }
  # Here I assume that all the sequences have the same length
  RD<-RangedData(ranges, space = chr, score=score)
  RD
}


MATScore<-function(tilingSet, cName=NULL, dMax=600, nProbesMin=8, dMerge=300, method="score", threshold=5, verbose=FALSE, bedName=NULL)
{
  .Defunct("computeMATScore",package="rMAT","This function is now Defunct. The MAT score function has been slipt into two seperated functions to facilitate export to wig/bed files, see the vignette for more details")
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

  if(verbose)
  {
    cat("** Calculating MAT Scores **\n")
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

  # seqNum<-rep(0,nProbes)
  # tmp<-unique(tilingSet@featureChromosome)
  # for(i in 1:length(tmp))
  # {
  #   seqNum[tilingSet@featureChromosome==tmp[i]]<-i
  # }
  seqNum<-as.numeric(as.factor(tilingSet@featureChromosome))
  
  
  obj<-.C("MAT",
  as.double(t(C)),
  as.double(t(I)),
  as.integer(nProbes),
  as.integer(nArraysC),
  as.integer(nArraysI),
  as.integer(tilingSet@featurePosition),
  as.double(dMax),
  as.double(dMerge),
  as.integer(nProbesMin),
  as.double(threshold),
  MATScore=double(nProbes),
  pValue=double(nProbes),
  as.integer(methodNum),
  regions=integer(nProbes), 
  as.integer(verbose),
  as.integer(seqNum),
  numRegions = integer(1),
  package="rMAT")


  MATScore<-obj$MATScore
  pValue<-obj$pValue
  ### Replace the zeros (correponding to regions with too few probes) by NA:
  pValue[MATScore==0]<-NA
  MATScore[MATScore==0]<-NA    

  if(verbose)
  {
    cat("** Finished processing ",nProbes," probes on ",nArraysC+nArraysI," arrays **\n");
    cat("** Number of Enriched regions is ", obj$numRegions, " **\n")
  }
  numRegions<-obj$numRegions

  if(numRegions > 0)
  {
    if(!is.null(bedName))
    {
      cat("** Creating bed file **\n")
      bed<-lapply(1:numRegions,".extractRegion",obj$regions,obj$MATScore,tilingSet@featurePosition,tilingSet@featureChromosome,nProbesMin)
      bed<-unlist(bed)
      bed<-matrix(bed,length(bed)/4,4,byrow=TRUE)
      write.table(unlist(bed),file=bedName,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    }
    else
    {
      if(verbose)
      {
        cat("You have not specified an output. No output written.\n")
      }
    }
  }
  else
  {
    if(verbose)
    {
      cat("** No regions to output **\n")
    }
  }

  NewMAT<-new('MAT', genomeName=tilingSet@genomeName, featurePosition=tilingSet@featurePosition, featureChromosome=tilingSet@featureChromosome, score=obj$MATScore, pValue=obj$pValue, regIndex=obj$regions, method=method, threshold=threshold)
  NewMAT
}

.extractRegion<-function(index,regions,MATScore,featurePosition,featureChromosome,nProbesMin)
{
  ind<-regions==index
  if(sum(ind)>=nProbesMin)
  {
    chr<-unique(featureChromosome[ind])
    pos<-unique(featurePosition[ind])
    score<-max(MATScore[ind])
    return(list(chr=chr,start=min(pos),end=max(pos),score=score))
  }
  else
  {
    return(NULL)
  }
}

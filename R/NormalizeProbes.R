NormalizeProbes<-function(tilingSet, method="MAT",robust=FALSE, all=TRUE, standard=TRUE, verbose=FALSE)
{ 

  ### Sanity checks 
  if(class(tilingSet)!="tilingSet")
  {
    stop("tilingSet must be an object of class tilingSet (e.g. returned by BPMAPCelParser)!")
  }
  if(preproc(tilingSet@experimentData)$transformation!="log")
  {
    stop("The probe measurements need to be log transformed!")
  }

  if(verbose)
  {
    cat("** Initializing NormalizeProbes **\n")
  }

  if(method=="MAT")
  {
    methodNum<-1
  }
  else if(method=="PairBinned")
  {
    methodNum<-2
  }

  y<-exprs(tilingSet)
  nProbes<-nrow(y)
  nArrays<-ncol(y)    
  Seq<-tilingSet@featureSequence

  ### Get the indices of duplicates
  indRepeats<-duplicated(Seq)
  nRepeats<-sum(indRepeats)

  if(nRepeats>0 & verbose)
  {
    cat("** Removing ", nRepeats," in the estimation procedures **", "\n")
  }

  set.seed(20)

  # If there are too many probes, it's not necessary to use them all for efficiency  
  ordSeq<-c(sample((1:nProbes)[!indRepeats],nProbes-nRepeats),(1:nProbes)[indRepeats])

  ### Here we add some noise so that the design matrix is full rank if the copy number is not specified
  if(length(unique(tilingSet@featureCopyNumber))==1 & unique(tilingSet@featureCopyNumber)[1]==1)
  {
    copyNumber<-tilingSet@featureCopyNumber+rnorm(length(tilingSet@featureCopyNumber),0,0.001)
  }
  else
  {
 	 copyNumber<-tilingSet@featureCopyNumber
  }
  
  if(verbose)
  {
    cat("** Normalization is starting **\n")
  }

  # Calling C code
  obj<-.C("NormalizeProbes",
  as.character(Seq[ordSeq]),
  as.double(y[ordSeq,]),
  yNormalized=as.double(rep(0,nArrays*nProbes)),
  as.integer(nProbes),
  as.integer(nArrays),
  copyNumber=as.double(copyNumber[ordSeq]),
  as.integer(methodNum),
  as.integer(robust),
  adjRSquare=double(nArrays),
  RSquare=double(nArrays),
  BIC=double(nArrays),
  beta=double(100*nArrays),
  betaLength=integer(1),
  as.integer(all),
  as.integer(standard),
  NAOK=FALSE,
  DUP=TRUE,
  as.integer(verbose),
  package="rMAT")


  sampleNames<-sampleNames(tilingSet@phenoData)
  if(verbose)
  {
    cat("** Finished Normalizing ",nProbes," probes on ",nArrays," arrays **\n")
    cat("** Adjusted R-Squares values between fitted and observed values: **\n") 
    for(i in 1:nArrays)
    {
      cat("** ", sampleNames[i],": ", obj$adjRSquare[i]," **\n")
    }
  }

  normData<-matrix(obj$yNormalized, nProbes, nArrays)[order(ordSeq),]
  if(nArrays==1)
  {
    # Reconvert to a matrix, in case there is only one array
    normData<-matrix(normData,nProbes,1)
  }
  ### Set the normalized data as exprs

  exprs(tilingSet)<-normData

  ### Setting the normalization parameters
  normalization<-method
  if(robust)
  {
    normalization<-paste(normalization, "robust",sep=" ")
  }
  if(standard)
  {
    normalization<-paste(normalization, "standardized",sep=" ")
  }

  preproc(tilingSet@experimentData)$normalization<-normalization
  tilingSet
}

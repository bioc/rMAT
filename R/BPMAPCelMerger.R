".BPMAPCelMerger"<-function(BPMAPList, combinedCels, verbose=TRUE )
{
  if(verbose)
  {	
    cat("** Order BPMAP by Y and X coordinates\n")
  }
  BPMAPorder<-order(BPMAPList$Y, BPMAPList$X) #ranking vector ordered by Y and X coordinates
  outBPMAP<- vector("list", length(BPMAPList))

  # order each vector by BPMAPorder
  for(i in 1:length(BPMAPList))
  {
    outBPMAP[[i]]<-cbind(BPMAPList[[i]], deparse.level=2)[BPMAPorder,]
  }
  
  names(outBPMAP)<-names(BPMAPList)
  if(verbose)
  {	
    cat("** Ordering CEL by Y and X coordinates **\n")	
  }
  
  CELorder<-order(combinedCels$Y, combinedCels$X) 
  outCEL<-vector("list", length(combinedCels))
  
  #ranking vector ordered by Y and X coordinates
  for(i in 1:length(combinedCels))
  {
    outCEL[[i]]<-cbind(combinedCels[[i]],deparse.level=2)[CELorder,]
  }
  
  names(outCEL)<-names(combinedCels)
  
  if(verbose)
  {	
    cat("** C++ merging function **\n")
  }
  out<-.Call("BPMAPCelMerger", outBPMAP, outCEL, package="rMAT")	

  if(verbose)
  {	
    cat("** Merged finished **\n")	
  }
  out
}

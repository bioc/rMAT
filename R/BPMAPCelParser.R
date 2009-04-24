"BPMAPCelParser" <- function(BPMAPFileName, CelFileNames, genomeName=NULL, verbose=FALSE, groupName="",seqName="")
{

  lArray<-length(CelFileNames)

  ### Read the header file
  bpmapHeader<-ReadBPMAPAllSeqHeader(BPMAPFileName)

  ### Only read the sequence that we need, here I filter the sequences that contain chr
  seqToRead<-as.integer(bpmapHeader$seqNum[bpmapHeader$GroupName==groupName])
  seqToRead<-seqToRead[grep(seqName,bpmapHeader$SeqName[bpmapHeader$GroupName==groupName])]

  SeqChr<-.ReadBPMAPSeq(BPMAPFileName, seqToRead, readPM= TRUE, readMM=FALSE, readProbeLength=FALSE, readPMProbe=TRUE, readMatchScore=TRUE, readPosition=TRUE, readTopStrand=FALSE, verbose=verbose)

  #Reading Cel files using affxparser
  tmpCel<-vector("list",length(CelFileNames)+2)

  phenoData<-character(lArray)  
  for (i in 1:lArray)
  {
    tmp<-strsplit(CelFileNames[i], "/")
    tmp<-strsplit((tmp[[1]])[length(tmp[[1]])], ".[cC][eE][lL]")
    phenoData[i]<-(tmp[[1]])[1]
    if(i==1)
    {
      tmp<-readCel(CelFileNames[1], readXY=TRUE, readIntensities=TRUE, readOutliers=TRUE)
      tmpCel[1]<-list(tmp$x)
      tmpCel[2]<-list(tmp$y)
      tmpCel[3]<-list(tmp$intensities)
      names(tmpCel)[1] <- "X"
      names(tmpCel)[2] <- "Y"
      names(tmpCel)[3] <- "I1"
    }
    #for the other CEL files we don't have to read the X and Y coordinate, since they were already read.
    else
    {
      tmpCel[i+2]<-list(readCel(CelFileNames[i], readXY=FALSE, readIntensities=TRUE, readOutliers=TRUE)$intensities)
    }
    names(tmpCel)[i+2]<-paste("I", i, sep="")
  }
  names(SeqChr)[1]<-"X"
  names(SeqChr)[2]<-"Y"

  combineData<-.BPMAPCelMerger(SeqChr,tmpCel,verbose=verbose)

  seqNumNameIndex<-grep("SeqNum",names(combineData))

  ##order by SeqNum, Position now:
  if(verbose)
  {
    cat("** Sorting Data first by Sequence number, then by Position **\n")
  }

  seqNameVector<-as.character(combineData[[seqNumNameIndex]])

  for(i in seqToRead)
  {
    tmp<-strsplit(as.character(bpmapHeader$SeqName[bpmapHeader$seqNum==i]),"chr")    
    seqNameVector[seqNameVector==as.character(i)]<-(tmp[[1]])[length(tmp[[1]])]    
  }

  SeqNumorder<-order(seqNameVector,combineData$Position)
  combineDataBySP<-vector("list",length(combineData))
  for(i in 1:length(combineData))
  {
    combineDataBySP[[i]]<-cbind(combineData[[i]], deparse.level=2)[SeqNumorder,]       
  }

  names(combineDataBySP)<-names(combineData)  
  combineDataBySP

  Data<-combineDataBySP[grep("I",names(combineDataBySP))]
  Data<-as.data.frame(Data)

  names(Data)<-phenoData

  if(is.null(genomeName))
  {
    tmp<-strsplit(BPMAPFileName, "/")
    genomeName<-(tmp[[1]])[length(tmp[[1]])]
    genomeName<-strsplit(genomeName, ".bpmap")[[1]]
  }
  ### The BPMAP file need to contains copy numbers from xMAN
  copyNumber<-combineDataBySP$MatchScore/min(combineDataBySP$MatchScore)
  ### All copy numbers should be integer valued
    
  if(prod(copyNumber==round(copyNumber)))
  {
    copyNumber<-rep(1,length(combineDataBySP$Position))
  }

  myDesc <- new("MIAME")
  preproc(myDesc)<-list(transformation="log", normalization="none")

  newSet<-new('tilingSet', featureChromosome=paste("chr",seqNameVector[SeqNumorder],sep=""),featurePosition=combineDataBySP$Position, 
  featureCopyNumber=as.integer(copyNumber), exprs=as.matrix(log(Data)), genomeName=genomeName, 
  featureSequence=combineDataBySP$PMProbe, experimentData=myDesc)

  return(newSet)
}



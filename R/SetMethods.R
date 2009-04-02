setMethod("show", "tilingSet",
function(object)
{
    cat("Object of class 'tilingSet'","\n")
    cat("This object contains an ExpressionSet and has the following additional slots: \n")
    cat("genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber, featureSequence\n")
}
)

setMethod("show", "MAT",
function(object)
{
    cat("Object of class 'MAT'","\n")
    cat("This object contains an ExpressionSet and has the following additional slots: \n")
    cat("genomeName, featureSequence, featurePosition, featureChromosome, featureCopyNumber, score, pValue, regIndex, method, threshold \n")
}
)


setMethod("summary", signature("tilingSet"),
    function(object) {
      cat("   Genome interrogated: ",sort(object@genomeName)," \n")
      tmp<-unique(object@featureChromosome)
      tmp<-unlist((strsplit(tmp,"chr")))
      chr<-tmp[grep(".",tmp)]
      chr<-c(sort(as.integer(chr[grep("[1-9]",chr)])),sort(chr[grep("[a-zA-z]",chr)]))
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

setMethod("summary", signature("MAT"),
    function(object) {
      cat("   Genome interrogated: ", sort(object@genomeName)," \n\n")
      tmp<-unique(object@featureChromosome)
      tmp<-unlist((strsplit(tmp,"chr")))
      chr<-tmp[grep(".",tmp)]
      chr<-c(sort(as.integer(chr[grep("[1-9]",chr)])),sort(chr[grep("[a-zA-z]",chr)]))
      cat("   Chromosome(s) interrogated: ")
      cat(chr,sep=", ")
      cat(" \n")
      cat("   Regions selected by", object@method, "with a threshold of ", object@threshold, " \n")      
      cat("     - Total number of regions detected: ",sum(object@regIndex>0)," \n")
      cat("     - Number of regions detected by chromosomes: \n")
      for(i in 1:length(chr))
      {
        cat(chr[i],": ", sum(object@regIndex[object@featureChromosome==paste("chr",chr[i],sep="")]>0), "  ",sep="")
      }
      cat("\n")
    }
)

setMethod("combine", signature(x="tilingSet", y="tilingSet"),
    function(x, y) {
        if(class(x)!=tilingSet | class(y)!=tilingSet)
        {
          stop("The two objects need to be of class tilingSet")
        }
        xNames<-sampleNames(x@phenoData)
        yNames<-sampleNames(y@phenoData)
        xyNames<-xNames
        # try to find the root of the name file
        for(i in 1:length(xNames))
        {
          Name<-unique(strsplit(c(xNames[i],yNames[i]),"[.-]"))
          if(length(Name)>1)
          {
            stop("Are you sure you are merging the right objects? Sample names don't seem to have a common root.")
          }
          xyNames[i]<-Name
        }

        xgenomeName<-x@genomeName
        ygenomeName<-y@genomeName
        Name<-unique(strsplit(c(xgenomeName,ygenomeName),"[.-]"))
        if(length(Name)>1)
        {
          stop("Are you sure you are merging the right objects? Genome Names names don't seem to have a common root.")
        }
        xygenomeName<-Name
        
        if(preproc(x@experimentData)!=preproc(y@experimentData))
        {
          stop("Are you sure you are merging the right objects? They don't have the same ")
        }
        
        # We need to separate the numeric and alpha sequences
        xyfeatureChromosome<-c(x@featureChromosome,y@featureChromosome)
        tmp<-unlist((strsplit(xyfeatureChromosome,"chr")))
        indNum<-grep("[1-9]",tmp)
        NumChr<-tmp[indNum]
        indAlpha<-grep("[a-zA-z]",tmp)        
        AlphaChr<-tmp[indAlpha]
        
        xyfeaturePosition<-c(x@featurePosition,y@featurePosition)
        NumPos<-xyfeaturePosition[indNum]
        AlphaPos<-xyfeaturePosition[indAlpha]

        # Order by chromosome, then position
        NumOrd<-order(NumChr,NumPos)
        AlphaOrd<-order(AlphaChr,AlphaPos)
        
        ord<-c(indNum[NumOrd],indAlpha[AlphaOrd])
        xyfeatureCopyNumber<-c(x@featurePosition,y@featureCopyNumber)
        xyfeatureSequence<-c(x@featureSequence,y@featureSequence)
        xyData<-rbind(exprs(x),exprs(y))
        
      
        newSet<-new('tilingSet', featureChromosome=xyfeatureChromosome[ord],featurePosition=xyfeaturePosition[ord], 
        featureCopyNumber=xyfeatureCopyNumber[ord], exprs=xyData[ord,], genomeName=xygenomeName, 
        featureSequence=xyfeatureSequence[ord], experimentData=x@experimentData)
        return(newSet)        
    }
)



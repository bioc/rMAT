".ReadBPMAPSeq"<- function (fileName, range = NULL, readPM = TRUE, readMM = FALSE, readProbeLength = FALSE, readPMProbe =TRUE, readMatchScore = TRUE, readPosition = TRUE, readTopStrand = FALSE, verbose=TRUE){

  if(is.null(range))
  {
    maxRange <-as.integer(0)
    minRange <-as.integer(0)
  }
  else
  {	
    range<-as.integer(range)
    maxRange <-max(range)
    minRange<-min(range)
  }

  #borrow from affxparser - allow us to expand '~' pathnames to full pathnames.
  fileName <- file.path(dirname(fileName), basename(fileName))
  if (!file.exists(fileName)) 
  {
    stop("Cannot read BPMAP file. File not found: ", fileName)
  }
  #Calling the helper C++ function to read the BPMAP file
  content<-.Call("readBPMAPSeq", fileName, range, readPM, readMM, readProbeLength, readPMProbe, readMatchScore, readPosition, readTopStrand, verbose, maxRange, minRange, package="rMAT")
  content
}

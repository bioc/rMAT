".ReadBPMAPHeader" <-function (fileName){

    #borrow from affxparser
      # Expand '~' pathnames to full pathnames.
      fileName <- file.path(dirname(fileName), basename(fileName))
      if (!file.exists(fileName)) 
      {
        stop("Cannot read BPMAP file. File not found: ", fileName)
      }
      
	content<-.Call("readBPMAPFileHeader", fileName, package="rMAT")
	content
}


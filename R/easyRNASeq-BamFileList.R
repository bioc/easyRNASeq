##' Get a BamFileList from a list of filenames
##'
##' A utility function to create a \code{linkS4class{BamFileList-class}}{BamFileList}
##' object from a set of filenames. The filenames need to contain the file path if they
##' are not in the working directory.

##' @name getBamFileList
##' @rdname easyRNASeq-BamFileList
##' @aliases getBamFileList,character-method
##' @param filenames a character vector containing fully defined filenames
##' @return a \code{linkS4class{BamFileList-class}}{BamFileList}
##' @seealso \code{linkS4class{BamFileList-class}}{BamFileList}
##' \code{\link[base:list.files]{dir}}
##' @examples
##' 
##'   library("RnaSeqTutorial")
##' 
##' 	## creating a BamFileList using a directory and pattern
##' 	bfl <- getBamFileList(
##'                         dir(path=system.file("extdata",
##' 					                  package="RnaSeqTutorial"),
##' 					              pattern="[A,C,T,G]{6}\\.bam$",
##'   				              full.names=TRUE))
##' 
##'   ## using filenames
##'   filenames <- dir(path=
##'   	    			      system.file("extdata",
##' 					          package="RnaSeqTutorial"),
##' 					     pattern="[A,C,T,G]{6}\\.bam$",
##'   				     full.names=TRUE)
##'   bfl <- getBamFileList(filenames)
##'   
##'   ## get them recursively
##'   filenames <- dir(path=
##'         			      system.file("extdata",
##' 					          package="RnaSeqTutorial"),
##' 					     pattern="[A,C,T,G]{6}\\.bam$",
##'   				     full.names=TRUE,recursive=TRUE)
##'   bfl <- getBamFileList(filenames)
##'
setMethod(f="getBamFileList",
          signature="character",
          definition=function(filenames=character(0)){
            
            ## check the files
            if(length(filenames) == 0){
              stop("You need to provide a character vector of 'filenames'. ")
            }
            
            ## check that the file(s) exist(s)
            sel <- file.exists(filenames)
            if(any(!sel)){
              stop(paste("The file(s)",
                         paste(filenames[!sel],collapse = " and ")),"do(es) not exist.")
            }            
            
            ## check that the bai are not part of the filesList
            if(length(grep("\\.bam\\.bai$",filenames))>0){
                warning(paste("You either have provided BAM index files (.bai) as",
                              "part of the filenames argument or your pattern matched",
                              "BAM index files! Removing these from the files",
                              "list to process."))
                filesList <- filesList[-grep("\\.bam\\.bai$",names(filesList))]
            }

            ## check if we have index with bai
            indexes <- paste(filenames,"bai",sep=".")
            sel <- file.exists(indexes)
            if(any(!sel)){
              stop(paste("Index files (bai) are required. They are missing for the files: ",
                         paste(filesList[!sel],collapse = " and "),
                         ". Use the Rsamtools indexBam function or the samtools index command line utility to create them.",sep=""))
            }

            ## create the BamFileList
            ## as in Rsamtools the index is the filename, i.e. without the .bai extension
            bfl <- BamFileList(filenames,index=filenames)
            names(bfl) <- basename(filenames)

            ## validate them
            validate(bfl)

            ## return            
            return(bfl)
          })

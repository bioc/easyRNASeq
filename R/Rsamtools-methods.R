##' Extension of the Rsamtools package
##' 
##' Describes extensions to the Rsamtools package.
##' \itemize{
##' \item For \code{\linkS4class{BamFile}} and
##' \code{\linkS4class{BamFileList}} objects: 
##' \itemize{
##' \item \code{validate} validates a \code{\linkS4class{BamFile}} or
##' \code{\linkS4class{BamFileList}} object.
##' }
##' }
##' 
##' \code{validate} checks whether the BAM file exists and if a BAI index is present.
##' 
##' @rdname Rsamtools-methods
##' @aliases validate validate,BamFile-method validate,BamFileList-method
##' @param obj An object of the \code{\linkS4class{BamFile}} or 
##' \code{\linkS4class{BamFileList}} class
##' @param header a boolean to (de)activate the check for a BAM
##' header
##' @param cross.validation a boolean - only valid for
##' \code{\linkS4class{BamFileList}} objects - to
##' (de)activate the cross validation of all the BAM files
##' header
##' @return 
##' \code{validate} returns invisibly a vector of boolean.
##' Fails anyway if any file is missing.
##' @author Nicolas Delhomme
##' @seealso
##' \itemize{
##' \item \code{\linkS4class{BamFile}}
##' \item \code{\linkS4class{BamFileList}}
##' }
##' @keywords methods
##' @examples
##' library(Rsamtools)
##' filenames <- dir(path=system.file("extdata",
##'       		                  package="RnaSeqTutorial"),
##' 					       pattern="[A,C,T,G]{6}\\.bam$",
##'   				       full.names=TRUE)
##' 
##' bfl <-BamFileList(filenames,index=filenames)
##' 
##' validate(bfl)
##' 
setMethod(f="validate",
          signature="BamFile",
          definition=function(obj,header=TRUE,cross.validation=TRUE){
            
            if(!file.exists(obj)){
              stop(paste("The file:",path(obj),"does not exist."))
            }
            
            if(!file.exists(paste(index(obj),"bai",sep="."))){
              stop(paste("The index file: ",path(obj),".bai does not exist.",
                         "Use the Rsamtools indexBam function or the samtools index command line utility to create it.",
                         sep=""))
            }
            
            if(header){
              if(length(seqinfo(obj))==0){
                stop("Your BAM file has no sequence information in its header.")
              }

              ## sorted by coordinate
              if(scanBamHeader(obj)$text$`@HD`[2] != "SO:coordinate"){
                stop("Your BAM file needs to be sorted by coordinate; see (R)samtools sort.")
              }              
            }
            invisible(TRUE)
          })

setMethod(f="validate",
          signature="BamFileList",
          definition=function(obj,header=TRUE,cross.validation=TRUE){
            
            ## check that there are bam files
            if(length(obj)==0){
              stop("You have not provided any BAM files. Your BamFileList object is empty.")
            }
            
            ## check all the BamFile
            sapply(obj,validate,header)
            
            ## then cross.validate
            if(cross.validation){
              seqinfos <- lapply(obj,seqinfo)
              if(any(!sapply(seqinfos,"identical",seqinfos[[1]]))){
                stop("Your BAM files have different headers.")
              }
            }
            invisible(TRUE)
          })


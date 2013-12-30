##' Fetch genic annotation from a gff/gtf file or using biomaRt
##' 
##' The annotation can be retrieved in two ways
##' \itemize{
##' \item{biomaRt}{Use biomaRt and Ensembl to get organism specific annotation.}
##' \item{gff/gtf}{Use a gff or gtf local annotation file.}}
##' When using \pkg{biomaRt}, it is
##' important that the \code{organismName} slot of the
##' \code{\linkS4class{RNAseq}} object is set the prefix of one of the value
##' available using the \pkg{biomaRt}
##' \code{\link[biomaRt:listDatasets]{listDatasets}} function, e.g.
##' "Dmelanogaster".  When reading from a gff/gtf file, a version 3 formatted
##' gff (gtf are modified gff3 from Ensembl) is expected. The function
##' \pkg{genomeIntervals} \code{\link[genomeIntervals:readGff3]{readGff3}} is
##' used to read the data in.
##' Another annotation caveat is the reference names, \emph{i.e.} the
##' chromosome/scaffold names used in the alignment files and those fetched
##' when retrieving the genic annotation might differ. \pkg{easyRNASeq} tries
##' to be clever in this case and guess the correspondance. However, it is
##' not always obvious. Organisms were this has been checked can be listed with
##' the \emph{knownOrganisms} function.
##' 
##' \dots{} are for additional arguments, passed to the \pkg{biomaRt}
##' \code{\link[biomaRt:getBM]{getBM}} function or to the
##' \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
##' internal function that takes an optional arguments: annotation.type that
##' default to "exon". This is used to select the proper rows of the gff or gtf
##' file.
##' 
##' @aliases fetchAnnotation fetchAnnotation,RNAseq-method
##' knownOrganisms knownOrganisms,missing-method
##' @name easyRNASeq annotation methods
##' @rdname easyRNASeq-annotation-methods
##' @param obj An object of class \code{RNAseq}
##' @param annotationMethod one of biomaRt, gff, gtf
##' @param filename If the method is gff or gtf, the actual gtf, gff filename
##' @param ignoreWarnings set to TRUE (bad idea! they have a good reason to be
##' there) if you do not want warning messages.
##' @param ... See details
## TODO complain to roxygen2 about \dots 
##' @return A \code{\linkS4class{RangedData}} containing the fetched
##' annotations.
##' @author Nicolas Delhomme
##' @keywords connection data methods
##' @examples
##' 
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	obj <- new('RNAseq',
##' 		organismName="Dmelanogaster",
##' 		readLength=36L,
##' 		chrSize=as.list(seqlengths(Dmelanogaster))
##' 		)
##' 
##' 	obj <- fetchAnnotation(obj,
##' 				method="gff",
##'                                 filename=system.file(
##' 						"extdata",
##' 						"annot.gff",
##' 						package="RnaSeqTutorial"))
##' 	}
##' 
setMethod(
          f="fetchAnnotation",
          signature="RNAseq",
          definition=function(obj,
            annotationMethod=c("biomaRt","gff","gtf"),
            filename=character(1),
            ignoreWarnings=FALSE,...){
            
            ## check if method was provided
            annotationMethod <- match.arg(annotationMethod)
            
            ## check if method was provided the old way
            dots <- list(...)
            if ("method" %in% names(dots)){
              
              ## get the methods
              methods <- eval(formals("fetchAnnotation")$annotationMethod)
    
              ## check the provided one
              if(dots$method %in% methods){
                warning(
                  "The 'method' argument to fetchAnnotation is deprecated. Use 'annotationMethod' instead.")
                annotationMethod <- dots$method
              }
            }
            
            ## switch depending on the annotationMethod
            exon.range <- switch(EXPR=annotationMethod,
                                 "biomaRt"={.getBmRange(organismName(obj),ignoreWarnings=ignoreWarnings,...)},
                                 "gff"={.getGffRange(organismName(obj),filename=filename,ignoreWarnings=ignoreWarnings,...)},
                                 "gtf"={.getGtfRange(organismName(obj),filename=filename,ignoreWarnings=ignoreWarnings,...)}
                                 )
            
            ## update the obj
            genomicAnnotation(obj)<-exon.range
            
            ## return
            return(obj)
          })

##' @aliases knownOrganisms
##' @rdname easyRNASeq-annotation-methods
##' @return A vector containing the known organisms
##' @author Nicolas Delhomme
##' @keywords data methods
##' @examples
##' knownOrganisms()
##' 
setMethod(
          f="knownOrganisms",
          signature="missing",
          definition=function(){
            eval(formals(easyRNASeq:::.convertToUCSC)$organism)
          })

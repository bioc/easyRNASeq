##' Defunct annotation function
##'
##' The \code{fetchAnnotation} and \code{knownOrganisms} function are now
##' defunct. The \code{fetchAnnotation} function has been replaced by the
##' \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{getAnnotation}} method.
##'
##' @aliases fetchAnnotation-defunct knownOrganisms-defunct
##' @name easyRNASeq defunct annotation methods
##' @rdname easyRNASeq-defunct-annotation-methods
##' @author Nicolas Delhomme
NULL

##' Get genic annotation from a gff3/gtf file or using biomaRt
##'
##' The annotation can be retrieved in two ways
##' \itemize{
##' \item{biomaRt}{Use biomaRt and Ensembl to get organism specific annotation.}
##' \item{gff3/gtf}{Use a gff3 or gtf local annotation file.}
##' }
##' \itemize{
##' \item{When using \pkg{biomaRt}, it is
##' important that the \code{organism} argument to
##' \code{\linkS4class{AnnotParam}} is set the prefix of one of the value
##' available using the \pkg{biomaRt}
##' \code{\link[biomaRt:listDatasets]{listDatasets}} function, e.g.
##' "Dmelanogaster".}
##' \item{When reading from a gff3/gtf file, a version 3 formatted
##' gff or a gtf (an Ensembl defined gff2 version) is expected. The function
##' \pkg{genomeIntervals} \code{\link[genomeIntervals]{genomeIntervals-readGff3}} is
##' used to import the data.}
##' }
##'
##' \dots{} are for additional arguments, passed to the \pkg{biomaRt}
##' \code{\link[biomaRt:getBM]{getBM}} function or to the
##' \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
##' internal function that takes an optional arguments: annotation.type that
##' default to "exon". This is used to select the proper rows of the gff or gtf
##' file.
##'
##' @aliases getAnnotation getAnnotation,AnnotParam-method
##' @name easyRNASeq annotation methods
##' @rdname easyRNASeq-annotation-methods
##' @param obj An object of class \code{AnnotParam}
##' @param verbose a boolean to turn on verbosity
##' @param ... See details
##' @return A \code{\linkS4class{GRanges}} containing the fetched
##' annotations.
##' @author Nicolas Delhomme
##' @keywords connection data methods
##' @examples
##'
##'   \dontrun{
##' 	library("RnaSeqTutorial")
##'   getAnnotation(
##'     AnnotParam(
##'       organism="Dmelanogaster",
##'       datasource=system.file(
##'   		  "extdata",
##' 				"Dmel-mRNA-exon-r5.52.gff3",
##' 				package="RnaSeqTutorial"),
##'   		type="gff3"
##'   ))
##' }
##'
setMethod(
  f="getAnnotation",
  signature="AnnotParam",
  definition=function(obj,verbose=FALSE,...){

    ## first validate
    if(verbose){
      message("Validating the annotation source")
    }
    .validate(obj,verbose=verbose)

    ## then extract
    if(verbose){
      message("Fetching the annotation")
    }
    .extract(obj,verbose=verbose,...)
  }
)

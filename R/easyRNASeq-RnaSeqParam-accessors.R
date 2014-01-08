##' Accessors for RnaSeqParam class
##' 
##' These functions and generics define `accessors` (to get and set values) for
##' \code{\linkS4class{RnaSeqParam}} objects within the \pkg{easyRNASeq} package.
##' Implemented are:
##' \itemize{
##' \item annotParam
##' \item bamParam
##' \item countBy
##' \item datasource
##' \item paired
##' \item precision
##' \item stranded
##' \item yieldSize
##' }
##' 
##' @aliases RnaSeqParam-accessors
##' yieldSize,RnaSeqParam-method
##' bamParam annotParam countBy precision
##' countBy,RnaSeqParam-method precision,RnaSeqParam-method
##' bamParam,RnaSeqParam-method annotParam,RnaSeqParam-method
##' datasource,RnaSeqParam-method 
##' paired,RnaSeqParam-method stranded,RnaSeqParam-method
##' @name easyRNASeq RnaSeqParam accessors
##' @rdname easyRNASeq-RnaSeqParam-accessors
##' @param object An object derived from class \code{RnaSeqParam}.
##' @return
##' The value of the corresponding slot.
##' @author Nicolas Delhomme
##' @keywords manip
##' @seealso \itemize{
##'   \item The \code{\link[easyRNASeq:easyRNASeq-AnnotParam-class]{AnnotParam}} class
##'   \item The \code{\link[easyRNASeq:easyRNASeq-BamParam-class]{BamParam}} class
##'   \item The \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-class]{RnaSeqParam}} class
##' }
##' The \code{\link[easyRNASeq:easyRNASeq-BamParam-accessors]{BamParam yieldSize}} accessor
##' @examples
##' ## create the necessary AnnotParam
##' annotParam <- AnnotParam(
##'                 datasource=system.file(
##'                                 "extdata",
##'                                 "Dmel-mRNA-exon-r5.52.gff3",
##'                                 package="RnaSeqTutorial"))
##' 
##' ## create the RnaSeqParam
##' rsp <- RnaSeqParam(annotParam=annotParam)
##' ## get the yieldSize Parameter
##' ysize <-yieldSize(rsp)
##' 

setMethod(f="yieldSize",
          signature="RnaSeqParam",
          definition=function(object){
            return(yieldSize(bamParam(object)))
          })

setMethod(f="precision",
          signature="RnaSeqParam",
          definition=function(object){
            object@precision
          })

setMethod(f="countBy",
          signature="RnaSeqParam",
          definition=function(object){
            object@countBy
          })

setMethod(f="bamParam",
          signature="RnaSeqParam",
          definition=function(object){
            object@bamParam
          })

setMethod(f="annotParam",
          signature="RnaSeqParam",
          definition=function(object){
            object@annotParam
          })

setMethod(f="datasource",
          signature="RnaSeqParam",
          definition=function(object){
            datasource(annotParam(object))
          })

## not exported
setMethod(f="paired",
          signature="RnaSeqParam",
          definition=function(object){
            return(paired(bamParam(object)))
          })

setMethod(f="stranded",
          signature="RnaSeqParam",
          definition=function(object){
            return(stranded(bamParam(object)))
          })


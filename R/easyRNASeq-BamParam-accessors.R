##' Accessors for BamParam class
##' 
##' These functions and generics define `accessors` (to get and set values) for
##' \code{\linkS4class{BamParam}} objects within the \pkg{easyRNASeq} package.
##' 
##' @aliases BamParam-accessors yieldSize yieldSize,BamParam-method
##' paired stranded paired,BamParam-method stranded,BamParam-method
##' @name easyRNASeq BamParam accessors
##' @rdname easyRNASeq-BamParam-accessors
##' @param object An object derived from class \code{BamParam}.
##' @param ... Additional parameter inherited from the 
##' \code{\link[Rsamtools:RsamtoolsFile-class]{Rsamtools package yieldSize function}}.
##' Ignored here.
##' @usage
##' yieldSize(object,...)
##' paired(object)
##' stranded(object)
##' @return
##' The value of the corresponding slot.
##' @author Nicolas Delhomme
##' @keywords manip
##' @seealso The \code{\link[easyRNASeq:easyRNASeq-BamParam-class]{BamParam}} class
##' The \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam yieldSize}} accessor
##' @examples
##' 
##' bp <- BamParam()
##' ## get the yieldSize Parameter
##' ysize <-yieldSize(bp)
##' 

setMethod(f="yieldSize",
          signature="BamParam",
          definition=function(object,...){
            return(object@yieldSize)
          })

setMethod(f="paired",
          signature="BamParam",
          definition=function(object){
            return(object@paired)
          })

setMethod(f="stranded",
          signature="BamParam",
          definition=function(object){
            return(object@stranded)
          })


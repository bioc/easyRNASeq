##' Accessors for AnnotParam class
##' 
##' These functions and generics define `accessors` (to get and set values) for
##' \code{\linkS4class{AnnotParam}} objects within the \pkg{easyRNASeq} package.
##' Implemented are:
##' \itemize{
##' \item datasource
##' \item type
##' }
##' 
##' @aliases AnnotParam-accessors datasource
##' datasource,AnnotParam-method type,AnnotParam-method
##' @name easyRNASeq AnnotParam accessors
##' @rdname easyRNASeq-AnnotParam-accessors
##' @param x An object derived from class \code{AnnotParam}.
##' @param object An object derived from class \code{AnnotParam}.
##' @usage datasource(object)
##' \S4method{type}{AnnotParam}(x)
##' @return
##' The value of the corresponding slot.
##' @author Nicolas Delhomme
##' @keywords manip
##' @seealso The \code{\link[easyRNASeq:easyRNASeq-AnnotParam-class]{AnnotParam}} class.
##' The type and organism generics are imported from the \code{\link[BSgenome:BSgenome-class]{BSgenome}} and
##' \code{\link[Biostrings:XString-class]{Biostrings}} package, respectively.
##' @examples
##' 
##' annot <- AnnotParam(
##'   datasource=system.file(
##'                         "extdata",
##'                         "Dmel-mRNA-exon-r5.52.gff3",
##'                         package="RnaSeqTutorial"))
##' ## get the datasource Parameter
##' datasource(annot)
##' 

setMethod(f="datasource",
          signature="AnnotParam",
          definition=function(object){
            object@datasource
          })

setMethod(f="type",
          signature="AnnotParam",
          definition=function(x){
            switch(class(x),
                   "AnnotParamCharacter"=x@type,
                   NULL)
          })

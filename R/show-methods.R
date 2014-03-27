##' Display the content of classes from the easyRNASeq package.
##' 
##' Display the content of a \code{\linkS4class{RNAseq}}, 
##' \code{\linkS4class{AnnotParam}}, \code{\linkS4class{BamParam}} or 
##' \code{\linkS4class{RnaSeqParam}} object.
##' 
##' 
##' @name show methods
##' @rdname show-methods
##' @aliases show,RNAseq-method show,AnnotParam-method show,BamParam-method
##' show,RnaSeqParam-method
##' @docType methods
##' @section Methods: \describe{ 
##' \item{list("signature(object = \"RNAseq\")")}{
##' Display the values of the different slots of the
##' \code{\linkS4class{RNAseq}} object.  } 
##' \item{Annot/Bam/RnaSeqParam}{The respective object settings.}}
##' @keywords methods
setMethod(
          f="show",
          signature="RNAseq",
          definition=function(object){
            print(object)
          })

setMethod(
  f="show",
  signature="AnnotParam",
  definition=function(object){
    print(object)
  })

setMethod(
  f="show",
  signature="BamParam",
  definition=function(object){
    print(object)
  })

setMethod(
  f="show",
  signature="RnaSeqParam",
  definition=function(object){
    print(object)
  })

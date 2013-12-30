##' Display the content of a RNAseq object
##' 
##' Display the content of a \code{\linkS4class{RNAseq}} object.
##' 
##' 
##' @name show methods
##' @rdname show-methods
##' @aliases show,RNAseq-method
##' @docType methods
##' @section Methods: \describe{ \item{list("signature(object = \"RNAseq\")")}{
##' Display the values of the different slots of the
##' \code{\linkS4class{RNAseq}} object.  } }
##' @keywords methods
setMethod(
          f="show",
          signature="RNAseq",          
          definition=function(object){
            print(object)
          })

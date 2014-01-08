##' Extend the basename function to display Rsamtools BamFile class basename
##' 
##' Display the basename of the bam file represented by a \code{\linkS4class{BamFile}} object.
##' 
##' @name basename methods
##' @aliases basename basename,BamFile-method basename,BamFileList-method
##' @rdname basename-methods
##' @docType methods
##' @param path an object of class \code{\linkS4class{BamFile}} or \code{\linkS4class{BamFileList}}
##' @section Methods: \describe{ \item{list("signature(object = \"BamFile\")")}{
##' Display the basename of the bam file linked to by a
##' \code{\linkS4class{BamFile}} object.  } }
##' @keywords methods
setMethod(f="basename",
          signature="BamFile",
          definition=function(path){
            basename(path(path))
            })

setMethod(f="basename",
          signature="BamFileList",
          definition=function(path){
              sapply(path,function(p){
                  basename(path(p))
              })
          })

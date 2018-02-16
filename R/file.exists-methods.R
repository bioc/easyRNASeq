#' Extend the file.exists function to check the path slot of a Rsamtools BamFile class for existence
#'
#' Check if the bam file represented by a \code{\linkS4class{BamFile}} object exists.
#'
#' @name file.exists methods
#' @aliases file.exists file.exists,BamFile-method
#' @rdname file.exists-methods
#' @param ... a \code{\linkS4class{BamFile}} object
#' @docType methods
#' @section Methods: \describe{ \item{list("signature(object = \"BamFile\")")}{
#' Checkk if the bam file linked to by a
#' \code{\linkS4class{BamFile}} object exists.  } }
#' @keywords methods
setMethod(f="file.exists",
          signature="BamFile",
          definition=function(...){
            file.exists(path(...))
            })

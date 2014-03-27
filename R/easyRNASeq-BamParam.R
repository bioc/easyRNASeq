##' BamParam constructor
##' 
##' This constructs a \code{\linkS4class{BamParam}} object.
##' The default parameters are derived from the currently most 
##' common RNA-Seq experimental use-case and are detailed below:
##' \itemize{
##' \item paired is TRUE, \emph{i.e.} paired-end sequencing is expected.
##' \item stranded is FALSE \emph{i.e.} stranded sequencing is not expected.
##' \item yieldSize is set to 1,000,000. This is the amount of reads iteratively 
##' processed from the bam file stream. It is a compromise between speed, 
##' process-parallelization and memory usage.
##' }
##' 
##' Calling the constructor without argument result in the
##' default parameter described above to be returned. Calling
##' the constructor with any  parameter will affect the value
##' of the selected parameters, leaving the other parameters
##' unaffected.
##' 
##' @aliases BamParam BamParam,ANY-method
##' @name easyRNASeq BamParam constructor
##' @rdname easyRNASeq-BamParam
##' @param paired TODO 
##' @param stranded TODO
##' @param yieldSize TODO
##' @examples
##' ## the defaults
##' BamParam()
##' 
##' ## change the default
##' BamParam(paired=FALSE)
##' BamParam(stranded=TRUE,yieldSize=1L)
##' 
##' 
setMethod(f="BamParam",
          signature="ANY",
          definition=function(
                              paired=TRUE,
                              stranded=FALSE,
                              yieldSize=1e6L){
            new("BamParam",
                paired=paired,
                stranded=stranded,
                yieldSize=yieldSize)
          })

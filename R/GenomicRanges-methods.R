##' Extension of the GenomicRanges package
##' 
##' Describes extensions to the \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} package.
##' For \code{\linkS4class{GRanges}} and 
##' \code{\linkS4class{GRangesList}} objects:
##' \itemize{
##' \item \code{colnames} returns the column name of a \code{\linkS4class{GRanges}} or
##' \code{\linkS4class{GRangesList}} object.
##' \item \code{unsafeAppend} appends two \code{\linkS4class{GAlignments}}
##'   object together bypassing most sanity checks. Faster than the standard \code{c} or
##'   \code{append} function.
##' }
##' 
##' \itemize{
##' \item \code{colnames} returns the actual column names of the elementMetadata slot of the
##' \code{\linkS4class{GRanges}} or \code{\linkS4class{GRangesList}} object.
##' The elementMetadata contains a \code{\linkS4class{DataFrame}} object used
##' to store additional information provided by the user, such as exon ID in
##' our case.
##'   \item \code{unsafeAppend} appends two \code{\linkS4class{GAlignments}} objects.
##' }
##' 
##' @aliases colnames colnames,GRanges-method colnames,GRangesList-method
##' unsafeAppend unsafeAppend,GAlignments,GAlignments-method
##' @name easyRNASeq GenomicRanges package extension
##' @rdname GenomicRanges-methods
##' @param x An object of the \code{\linkS4class{GRanges}} or
##' \code{\linkS4class{GRangesList}} class
##' @param do.NULL see \code{\link[BiocGenerics:colnames]{colnames}} for details
##' @param prefix see \code{\link[BiocGenerics:colnames]{colnames}} for details
##' @param obj1 A \code{\linkS4class{GAlignments}} object
##' @param obj2 A \code{\linkS4class{GAlignments}} object
##' @usage colnames(x, do.NULL = TRUE, prefix = "col")
##' unsafeAppend(obj1,obj2)
##' @return \itemize{
##' \item \code{colnames}: A vector of column names.
##' \item \code{unsafeAppend}: A \code{\linkS4class{GAlignments}} object
##' }
##' @author Nicolas Delhomme
##' @seealso
##' \itemize{
##' \item \code{\linkS4class{DataFrame}}
##' \item \code{\linkS4class{GRanges}}
##' \item \code{\linkS4class{GRangesList}}
##' \item \code{\linkS4class{GAlignments}}
##' \code{\link[BiocGenerics:colnames]{colnames}}
##' }
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	## an example of a RangedData annotation
##' 	gAnnot <- RangedData(
##'                      IRanges(
##'                              start=c(10,30,100),
##'                              end=c(21,53,123)),
##'                           space=c("chr01","chr01","chr02"),
##'                           strand=c("+","+","-"),
##'                           transcripts=c("trA1","trA2","trB"),
##'                           gene=c("gA","gA","gB"),
##'                           exon=c("e1","e2","e3"),
##'                           universe = "Hs19"
##'                           )
##' 
##' 	## an example of a GRangesList annotation
##' 	grngs <- as(gAnnot,"GRanges")
##' 
##' 	## accessing the colnames
##' 	colnames(grngs)
##' 
##' 	## creating a GRangesList
##' 	grngsList<-split(grngs,seqnames(grngs))
##' 
##' 	## accessing the colnames
##' 	colnames(grngsList)
##' 	}
##'
##' ## For unsafeAppend
##' library(GenomicAlignments)
##' unsafeAppend(GAlignments(),GAlignments())
##'

## to extend GenomicRanges
##' @exportMethod colnames
setMethod(
  f="colnames",
  signature="GRanges",
  definition=function(x, do.NULL = TRUE, prefix = "col"){
    if (length(x) == 0){
      return(character(0))
    } else {              
      colnames(elementMetadata(x), do.NULL = do.NULL, prefix = prefix)
    }
  })

setMethod(
          f="colnames",
          signature="GRangesList",
          definition=function(x, do.NULL = TRUE, prefix = "col"){
            if (length(x) == 0){
              return(character(0))
            } else {              
              colnames(x[[1]], do.NULL = do.NULL, prefix = prefix)              
            }
          })

setMethod(f="unsafeAppend",
          signature=c("GAlignments","GAlignments"),
          definition=function(obj1,obj2){
            
            seqnames <- c(seqnames(obj1),seqnames(obj2))
            
            if(identical(seqlengths(obj1),seqlengths(obj2))){
              seqlengths <- seqlengths(obj1)
            } else {
              seqlengths <- c(seqlengths(obj1),seqlengths(obj2))
              seqlengths <- seqlengths[match(levels(seqnames),names(seqlengths))]
            }
            
            return(GAlignments(
              seqnames=seqnames,
              pos = c(start(obj1),start(obj2)),
              cigar = c(cigar(obj1),cigar(obj2)),
              strand = c(strand(obj1),strand(obj2)),
              names = NULL,
              seqlengths = seqlengths)
            )
          })

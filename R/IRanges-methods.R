# TODO this should be deprecated
#' Extension of the IRanges package
#'
#' Return the ranges of the genomic annotation.
#'
#' It retrieves the object stored in the genomicAnnotation slot of the RNAseq
#' object and apply the \code{ranges} function on it.
#'
#' @aliases ranges ranges,RNAseq-method
#' @name IRanges additional methods
#' @rdname IRanges-methods
#' @param x An object of the \code{\linkS4class{RNAseq}} class
#' @return An \code{\linkS4class{IRangesList}} object, where the split
#' is performed by seqnames (\emph{e.g.} chromosomes).
#' @author Nicolas Delhomme
#' @keywords methods
#' @examples
#'
#' 	\dontrun{
#' 	library("RnaSeqTutorial")
#'
#' 	obj <- getAnnotation(
#'             AnnotParam(
#'               organism="Dmelanogaster",
#'               datasource=system.file(
#'     	          "extdata",
#' 				        "Dmel-mRNA-exon-r5.52.gff3",
#' 				        package="RnaSeqTutorial"),
#'   		        type="gff3"
#'   ))
#'
#' 	ranges(obj)
#' 	}
#'
# extend IRanges
#' @exportMethod ranges
setMethod(
          f="ranges",
          signature="RNAseq",
          definition=function(x){
                     split(ranges(genomicAnnotation(x)),
                           seqnames(genomicAnnotation(x)))
          })


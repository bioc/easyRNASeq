##' Extension of the IRanges package
##' 
##' Return the ranges of the genomic annotation.
##' 
##' It retrieves the object stored in the genomicAnnotation slot of the RNAseq
##' object and apply the \code{ranges} function on it. The object retrieved can
##' be of the \code{\linkS4class{RangedData}} or
##' \code{\linkS4class{GRangesList}} class.
##' 
##' @aliases ranges ranges,RNAseq-method
##' @name IRanges additional methods
##' @rdname IRanges-methods
##' @param x An object of the \code{\linkS4class{RNAseq}} class
##' @return An \code{\linkS4class{IRangesList}} object, where the split
##' is performed by seqnames (\emph{e.g.} chromosomes).
##' @author Nicolas Delhomme
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	obj <- new('RNAseq',
##' 		organismName="Dmelanogaster",
##' 		readLength=36L,
##' 		chrSize=as.list(seqlengths(Dmelanogaster))
##' 		)
##' 
##' 	obj <- fetchAnnotation(obj,
##' 				method="gff",
##'                                 filename=system.file(
##' 						"extdata",
##' 						"annot.gff",
##' 						package="RnaSeqTutorial"))
##' 	ranges(obj)
##' 	}
##' 
## extend IRanges
##' @exportMethod ranges
setMethod(
          f="ranges",
          signature="RNAseq",
          definition=function(x){
            switch(class(genomicAnnotation(x)),
                   "GRanges"=split(ranges(genomicAnnotation(x)),
                                   seqnames(genomicAnnotation(x))), 
                   "RangedData"=ranges(genomicAnnotation(x))
            )
          })


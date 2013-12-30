##' Extension of the GenomicRanges package
##' 
##' Return the column name of a \code{\linkS4class{GRanges}} or
##' \code{\linkS4class{GRangesList}} object.
##' 
##' It returns the actual column names of the elementMetadata slot of the
##' \code{\linkS4class{GRanges}} or \code{\linkS4class{GRangesList}} object.
##' The elementMetadata contains a \code{\linkS4class{DataFrame}} object used
##' to store additional information provided by the user, such as exon ID in
##' our case.
##' 
##' @aliases colnames colnames,GenomicRanges-method colnames,GRangesList-method
##' @name GenomicRanges additional methods
##' @rdname GenomicRanges-methods
##' @param x An object of the \code{\linkS4class{GRanges}} or
##' \code{\linkS4class{GRangesList}} class
##' @param do.NULL see \code{\link[BiocGenerics:colnames]{colnames}} for details
##' @param prefix see \code{\link[BiocGenerics:colnames]{colnames}} for details
##' @usage colnames(x, do.NULL = TRUE, prefix = "col")
##' @return A vector of column names.
##' @author Nicolas Delhomme
##' @seealso
##' \code{\linkS4class{DataFrame}}
##' \code{\linkS4class{GRanges}}
##' \code{\linkS4class{GRangesList}}
##' \code{\link[BiocGenerics:colnames]{colnames}}
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
##'                           transcript=c("trA1","trA2","trB"),
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
## to extend GenomicRanges
##' @exportMethod colnames
setMethod(
          f="colnames",
          signature="GenomicRanges",
          definition=function(x, do.NULL = TRUE, prefix = "col"){
            if (length(x) == 0){
              return(character(0))
            } else {              
              return(switch(class(x),
                            "GRanges" = {
                              colnames(elementMetadata(x), do.NULL = do.NULL, prefix = prefix)
                            },
                            "GRangesList" = {
                              colnames(x, do.NULL = do.NULL, prefix = prefix)
                            },
                            {stop(paste("Unknown class:",class(x)))}))
            }
          })

## R 2.15.0 introduced the GRangesList signature
## TODO not sure if the switch in the method above is still required
setMethod(
          f="colnames",
          signature="GRangesList",
          definition=function(x, do.NULL = TRUE, prefix = "col"){
            if (length(x) == 0){
              return(character(0))
            } else {              
              colnames(elementMetadata(x[[1]]), do.NULL = do.NULL, prefix = prefix)              
            }
          })


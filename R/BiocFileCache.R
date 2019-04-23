##' Manages the data necessary for the examples using BiocFileCache
##'
##' Manages the tutorial and example data using the
##' \code{\linkS4class{BiocFileCache}} package
##'
##' @name BiocFileCache methods
##' @rdname BiocFileCache-methods
##' @aliases .get_cache .get_cache,ANY-method fetchData
##' fetchData,character-method tutorialData tutorialData,ANY-method
##' @docType methods
##' @param fileURL The URL of the file to retrieve. Alternatively, the ID of the
##' file in the BiocFileCache (i.e. the file basename), can be used.
##' @param \dots unused for the time being
##' @section Methods: \describe{
##' \item{.get_cache}{internal function to set up the cache}
##' \item{fetchData}{A function to fetch tutorial data, a file at a time}
##' \item{tutorialData}{the function to retrieve all the tutorial data and
##' cache it, if it is not already available}}
##' @keywords methods
##' @usage fetchData(fileURL)
##' @seealso \code{\linkS4class{BiocFileCache}}
##' @examples
##' tdir <- tutorialData()
##' gAnnot.path <- fetchData("gAnnot.rda")
setMethod(
  f=".get_cache",
  signature="ANY",
  definition=function(...){
      cache <- rappdirs::user_cache_dir(appname="easyRNASeq")
      BiocFileCache::BiocFileCache(cache,ask=FALSE)
  })

setMethod(
    f="fetchData",
    signature="character",
    definition=function(
        fileURL="https://github.com/UPSCb/UPSCb/raw/master/tutorial/easyRNASeq/gAnnot.rda"){

        bfc <- .get_cache()
        filename <- basename(fileURL)
        rid <- bfcquery(bfc,filename,"rname",exact=TRUE)$rid
        if (!length(rid)) {
            rid <- names(bfcadd(bfc, filename, fileURL ))
        }
        if (!isFALSE(bfcneedsupdate(bfc, rid))){bfcdownload(bfc, rid, ask=FALSE)}
        bfcrpath(bfc, rids = rid)
    })

setMethod(
    f="tutorialData",
    signature="ANY",
    definition = function(...){
      tdat <- sapply(TUTORIAL.DATA,fetchData)
      return(dirname(tdat[1]))
    })


#' Count methods for RNAseq object
#'
#' Summarize the read counts per exon, feature, gene, transcript or island.
#' \itemize{
#' \item{\code{exonCounts}: for that summarization, reads are
#' summarized per exons. An "exon" field is necessary in the annotation object
#' for this to work. See
#' \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq annotation
#' methods}} for more details on the annotation object.}
#' \item{\code{featureCounts} is similar to the 'exons' one. This is just a
#' wrapper to summarize count for genomic features that are not exon related.
#' I.e. one could use it to measure eRNAs. Again, a "feature" field is
#' necessary in the annotation object for this to work.}
#' \item{\code{geneCounts} sums the counts per either \code{bestExons} or
#' \code{geneModels}. In either case, the annotation object needs to contain
#' both an "exon" and a "gene" field.}
#' \item{\code{islandCounts} sums the
#' counts per computed islands.}
#' \item{\code{transcriptCounts} sums the counts
#' obtained by exons into their respective transcripts. Note that this often
#' result in counting some reads several times. For this function to work you
#' need both an "exon" and a "transcript" field in your annotation object. To
#' avoid this, one could create transcript specific synthetic exons, i.e.
#' features that would be unique to a transcript. To offer this possibility,
#' transcripts count can be summarized from "features", in which case the
#' annotation object need to have both the "feature" and "transcript" fields
#' defined.  }
#' }
#'
#' \dots{} for \itemize{
#' \item{geneCounts: additional options for the
#' \code{\link[easyRNASeq:easyRNASeq-summarization-internal-methods]{.geneModelSummarization}}}
#' \item{islandCounts: additional options for
#' \code{\link[easyRNASeq:easyRNASeq-island-methods]{findIslands}} }}
#'
#' @aliases exonCounts exonCounts,RNAseq-method
#' featureCounts featureCounts,RNAseq-method
#' geneCounts geneCounts,RNAseq-method
#' islandCounts islandCounts,RNAseq-method
#' transcriptCounts transcriptCounts,RNAseq-method
#' @name easyRNASeq summarization methods
#' @rdname easyRNASeq-summarization-methods
#' @param obj An object derived from class \code{\linkS4class{RNAseq}},can be
#' a \code{matrix} for RPKM, see details
#' @param force For \code{islandCount}, force RNAseq to redo \code{findIsland}
#' @param from either "exons" or "features" can be used to summarize per
#' transcript
#' @param summarization Method use for summarize genes
#' @param ... See details
#' @usage exonCounts(obj)
#' featureCounts(obj)
#' transcriptCounts(obj,from="exons")
#' geneCounts(obj,summarization=c("bestExons","geneModels"),...)
#' islandCounts(obj,force=FALSE,...)
#' @return A numeric vector containing count per exon, feature, gene or
#' transcript.
#' @author Nicolas Delhomme
#' @seealso \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq
#' annotation methods}}
#' \code{\link[easyRNASeq:easyRNASeq-summarization-internal-methods]{.geneModelSummarization}}
#' \code{\link[easyRNASeq:easyRNASeq-island-methods]{findIslands}}
#' @keywords methods
#' @examples
#'
#' library(curl)
#' library(BSgenome.Dmelanogaster.UCSC.dm3)
#'
#'  # get the example data files - we retrieve a set of example bam files
#'  # from GitHub using curl, as well as their index.
#' invisible(sapply(c("ACACTG","ACTAGC"),function(bam){
#'      curl_download(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'                           "master/tutorial/easyRNASeq/",bam,".bam"),paste0(bam,".bam"))
#'      curl_download(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'                           "master/tutorial/easyRNASeq/",bam,".bam.bai"),paste0(bam,".bam.bai"))
#'  }))
#'
#'  # get an example annotation file - we retrieve it from GitHub using curl
#'  invisible(curl_download(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'        "master/tutorial/easyRNASeq/gAnnot.rda"),"gAnnot.rda"))
#'
#'
#' 	# create an RNAseq object
#' 	# summarizing 2 bam files by exons
#' 	rnaSeq <- easyRNASeq(".",
#'                      organism="Dmelanogaster",
#'                      chr.sizes=seqlengths(Dmelanogaster),
#'                      readLength=36L,
#'                      annotationMethod="rda",
#'                      annotationFile="gAnnot.rda",
#'                      format="bam",
#'                      count="exons",
#'                      pattern="[A,C,T,G]{6}\\.bam$",
#'                      outputFormat="RNAseq")
#' 	# summing up the exons by transcript
#' 	rnaSeq <- transcriptCounts(rnaSeq)
#'
## TODO for the exonCounts and other counts, we really need to pay attention to the order of the chromosomes.
## We could have that forced (i.e. having them ordered) at the rnaSeq creation time
## and have a validity action to re-order them any time the field get changed
## TODO check if the above is true or if it is already handled by IRanges and co.
## TODO think of variable width reads

## Count methods
# exons
# TODO maybe for reporting, we should simplify the list
# at least mention that some exons are duplicated
# TODO see if we can extract the aggregate to an internal f and use it for exon and island
# TODO check with alejandro and simon for exons
# TODO think if we can add a function for the user to calculate the coverage

# TODO think if we can have the genomicAnnotation re-setting somewhere else

setMethod(
          f="exonCounts",
          signature="RNAseq",
          definition=function(obj){

            genomicAnnotation(obj) <- switch(class(genomicAnnotation(obj)),
                                             "GRanges"=genomicAnnotation(obj)[seqnames(genomicAnnotation(obj)) %in% names(readCoverage(obj))],
                                             "RangedData"=genomicAnnotation(obj)[names(ranges(obj)) %in% names(readCoverage(obj))])
            if(class(genomicAnnotation(obj))=="GRanges"){
              seqlevels(genomicAnnotation(obj)) <- seqlevels(genomicAnnotation(obj))[seqlevels(genomicAnnotation(obj)) %in% names(readCoverage(obj))]
            }
            exCounts <- .doBasicCount(obj)

            # FIXME - names are not ordered in the right way!!!
            names(exCounts) <- .getName(obj,"exons")
            readCounts(obj)<-.extendCountList(readCounts(obj),exCounts,"exons",filename=fileName(obj))

            # return
            return(obj)
          })

# features
setMethod(
          f="featureCounts",
          signature="RNAseq",
          definition=function(obj){

            genomicAnnotation(obj) <- switch(class(genomicAnnotation(obj)),
                                             "GRanges"=genomicAnnotation(obj)[seqnames(genomicAnnotation(obj)) %in% names(readCoverage(obj))],
                                             "RangedData"=genomicAnnotation(obj)[names(ranges(obj)) %in% names(readCoverage(obj))])
            if(class(genomicAnnotation(obj))=="GRanges"){
              seqlevels(genomicAnnotation(obj)) <- seqlevels(genomicAnnotation(obj))[seqlevels(genomicAnnotation(obj)) %in% names(readCoverage(obj))]
            }
            fCounts <- .doBasicCount(obj)
            names(fCounts) <-  .getName(obj,"features")
            readCounts(obj)<-.extendCountList(readCounts(obj),
                                             fCounts,
                                             "features",
                                             filename=basename(fileName(obj)))

            # return
            return(obj)
          })

# transcripts
setMethod(
          f="transcriptCounts",
          signature="RNAseq",
          definition=function(obj,from="exons"){

            # check if the from is valid
            if(! from %in% c("exons","features")){
              stop("To work, the 'from' arguments should be one of 'exons', 'features'")
            }

            # check if we have the exons summary already
            if(is.null(readCounts(obj,from))){
              obj <- switch(from,
                            "exons" = exonCounts(obj),
                            "features" = featureCounts(obj)
                            )
            }

            # summarize these counts
            tAgg <- aggregate(readCounts(obj,from),list(transcript=.getName(obj,"transcripts")),sum)
            tCounts<-tAgg[,-1,drop=FALSE]
            rownames(tCounts)<-tAgg$transcript
            readCounts(obj)<-.extendCountList(readCounts(obj),tCounts,"transcripts",filename=fileName(obj))

            # return
            return(obj)
          })

# genes
setMethod(
          f="geneCounts",
          signature="RNAseq",
          definition=function(obj,summarization=c("bestExons","geneModels"),...){

            # check the summarization methods
            summarizations <- eval(formals("geneCounts")$summarization)
            if(!summarization %in% summarizations){
              stop(paste(
                         "The given summarization:",
                         summarization,
                         "is not part of the supported summarizations:",
                         paste(summarizations,collapse=", ")))
            }

            # switch
            gCounts <- switch(EXPR=summarization,
                              "bestExons"={.bestExonSummarization(obj)},
                              "geneModels"={

                                # check if we can at all run
                                if(length(readCoverage(obj))==0){
                                  if(length(fileName(obj))>0){
                                    stop("Cannot execute the 'geneCounts' on a multiple sample 'RNAseq' object yet. Re-run the 'easyRNASeq' function with the 'count' argument 'genes' and the 'summarization' argument 'geneModels'.")
                                  } else {
                                    stop("No coverage information available. Use either the 'easyRNASeq' or 'fetchCoverage' methods first.")
                                  }
                                }

                                # TODO do we need that?
                                # it was generated already...
                                # and it does not work for multiple files...
                                # check if we have the gene model already
                                if(nrow(geneModel(obj))==0){
                                  geneModel(obj) <- .geneModelAnnotation(genomicAnnotation(obj),...)
                                }
                                .geneModelSummarization(obj)
                              })
            readCounts(obj)<-.extendCountList(readCounts(obj),gCounts,"genes",subType=summarization,filename=fileName(obj))

            # return
            return(obj)
          })

# islands
setMethod(
          f="islandCounts",
          signature="RNAseq",
          definition=function(obj,force=FALSE,...){

            # check if readIsland exists, otherwise compute it
            if(is.null(readIslands(obj)$names)|force){
              obj <- findIslands(obj,...)
            }

            # summarize using it
            islands=as.integer(unlist(
              aggregate(
                        readCoverage(obj)[match(names(readIslands(obj)),names(readCoverage(obj)))],
                        ranges(readIslands(obj)),
                        sum
                        )))

            iCounts <- ceiling(islands/readLength(obj))
            names(iCounts) <- readIslands(obj)$names
            readCounts(obj)<-.extendCountList(readCounts(obj),iCounts,"islands",filename=fileName(obj))

            # return
            return(obj)
          })

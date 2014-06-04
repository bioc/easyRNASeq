##' RnaSeqParam constructor
##' 
##' TODO ALL WRONG!!!
##' This constructs a \code{\linkS4class{RnaSeqParam}} object.
##' The default parameters are derived from the currently most 
##' common RNA-Seq experimental use-case and are detailed below:
##' \itemize{
##' \item paired is TRUE, \emph{i.e.} paired-end sequencing is expected.
##' \item stranded is FALSE \emph{i.e.} stranded sequencing is not expected.
##' \item yieldSize is set to 100,000. This is the amount of reads iteratively 
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
##' @aliases RnaSeqParam RnaSeqParam,ANY-method
##' @name easyRNASeq RnaSeqParam constructor
##' @rdname easyRNASeq-RnaSeqParam
##' @param annotParam An object derived from class \code{AnnotParam}.
##' @param bamParam An object derived from class \code{BamParam}.
##' @param countBy TODO
##' @param precision A character value, either 'read' or 'bp' that
##' defines the precision at which counting is done, either per read 
##' or per covered bp. 'read' is the default.
##' @examples
## create the necessary AnnotParam
##' annotParam <- AnnotParam(
##'                 datasource=system.file(
##'                                 "extdata",
##'                                 "Dmel-mRNA-exon-r5.52.gff3",
##'                                 package="RnaSeqTutorial"))
##' 
##' ## create the RnaSeqParam
##' rsp <- RnaSeqParam(annotParam=annotParam)
##' 
##' ## change some defaults
##' RnaSeqParam(countBy="features",annotParam=annotParam)
##' RnaSeqParam(bamParam=BamParam(stranded=TRUE,yieldSize=1L),annotParam=annotParam)
##' 
##' 
setMethod(f="RnaSeqParam",
          signature="ANY",
          definition=function(
            annotParam=AnnotParam(),
            bamParam=BamParam(),
            countBy=c("exons","features","genes","transcripts"),
            precision=c("read","bp")){
            

            precision <- match.arg(precision)
            countBy <- match.arg(countBy)
              
            ## what was that for anyway????
#             countBy <- switch(precision,
#                               "read"="transcripts",
#                               "bp"="exons"
#                    )
#             
            new("RnaSeqParam",
                annotParam=annotParam,
                bamParam=bamParam,
                countBy=countBy,
                precision=precision)
          })

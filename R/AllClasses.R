# TODO Remove the RNAseq class at next iteration
# Defunct the IRanges ranges method on it

##==========================
## RNAseq
##==========================
#' Class "RNAseq"
#'
#' A class holding all the necessary information and annotation to summarize
#' couts (number of reads) per features (i.e. exons or transcripts or genes)
#' for RNA-Seq experiments.
#'
#'
#' @name RNAseq class
#' @rdname easyRNASeq-class
#' @aliases RNAseq-class RNAseq
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the
#' form \code{new("RNAseq", ...)}.
#' @author Nicolas Delhomme
#' @seealso \itemize{ \item\code{\linkS4class{GRangesList}}
#' \item\code{\linkS4class{RleList}}
#' \item\code{\link{easyRNASeq-easyRNASeq}
#' \item\code{\link{easyRNASeq-accessors}
#' \item\code{\link{easyRNASeq-annotation-methods}
#' \item\code{\link{easyRNASeq-correction-methods}
#' \item\code{\link{easyRNASeq-coverage-methods}
#' \item\code{\link{easyRNASeq-summarization-methods}
#' \item\code{\link{print-methods}
#' }
#' @keywords classes
#' @examples
#'
#' showClass("RNAseq")
#'
#' @exportClass RNAseq
setClass(
         Class="RNAseq",
         representation=representation(
           chrSize="integer",
           fileName="character",
           geneModel="GRanges",
           genomicAnnotation="GRanges",
           librarySize="numeric",
           organismName="character",
           readCounts="list",
           readCoverage="RleList",
           readIslands="GRanges",
           readLength="integer"
           ),
         prototype=prototype(
           chrSize=integer(0),
           fileName=character(0),
           geneModel=GRanges(),
           genomicAnnotation=GRanges(),
           librarySize=numeric(0),
           organismName=character(1),
           readCounts=list(),
           readCoverage=RleList(),
           readIslands=GRanges(),
           readLength=integer(1)
           )
#         ,
#          validity=function(obj){
#          }
         )

##==========================
## BamParam
##==========================

#' Class "BamParam"
#'
#' A class describing the parameters of a bam file issued
#' from an RNA-Seq experiment.
#'
#' @name BamParam class
#' @rdname easyRNASeq-BamParam-class
#' @aliases BamParam-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the
#' form \code{new("BamParam", ...)} or using the BamParam constructor.
#' @section Slots from the Class: The \code{\linkS4class{BamParam}} class has
#' the following slots:
#' \itemize{
#' \item paired
#' \item stranded
#' \item strandProtocol
#' \item yieldSize
#' }
#' all of which can be accessed using the accordingly names accessor.
#' @author Nicolas Delhomme
#' @seealso \itemize{
#' \item \code{\link[easyRNASeq:easyRNASeq-BamParam-accessors]{BamParam accessors}}
#' \item \code{\linkS4class{RnaSeqParam}}
#' \item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam]{RnaSeqParam constructor}}
#' \item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam accessors}}
#' \item \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}
#' \item \code{\linkS4class{AnnotParam}}
#' \item \code{\link[easyRNASeq:easyRNASeq-AnnotParam]{AnnotParam constructor}}
#' }
#' @keywords classes
#' @examples
#'
#' showClass("BamParam")
#'
#' @exportClass BamParam
setClass("BamParam",
         representation=representation(
           paired="logical",
           stranded="logical",
           strandProtocol="character",
           yieldSize="integer"),
         prototype=prototype(
           paired=TRUE,
           stranded=FALSE,
           yieldSize=1e6L
         ))

##==========================
## AnnotParam
##==========================

#' Class "AnnotParam"
#'
#' A class holding all the necessary parameters to retrieve the necessary
#' annotation for processing an RNA-Seq experiment.
#'
#' @name AnnotParam class
#' @rdname easyRNASeq-AnnotParam-class
#' @aliases AnnotParam-class AnnotParamCharacter-class AnnotParamObject-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the
#' form \code{new("AnnotParamCharacter", ...)} or \code{new("AnnotParamObject", ...)}
#' (both subject to API changes) or using the
#' \code{\link[easyRNASeq:easyRNASeq-AnnotParam]{AnnotParam}} constructor (failsafe, prefered).
#' The class \code{AnnotParam} in itself is virtual and hence cannot be instantiated.
#' @author Nicolas Delhomme
#' @seealso \itemize{
#' \item \code{\linkS4class{RnaSeqParam}}
#' \item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam]{RnaSeqParam constructor}}
#' \item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam accessors}}
#' \item \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}
#' \item \code{\link[easyRNASeq:easyRNASeq-AnnotParam]{AnnotParam constructor}}
#' }
#' @keywords classes
#' @examples
#'
#' showClass("AnnotParam")
#'
#' @exportClass AnnotParam
setClass("AnnotParam")

#' @exportClass AnnotParamCharacter
setClass("AnnotParamCharacter",
         representation=representation(
           datasource="character",
           type="character"
         ),
         prototype=prototype(
           datasource=character(0),
           type="gff3"
         ),
         contains="AnnotParam")

#' @exportClass AnnotParamObject
setClass("AnnotParamObject",
         representation=representation(
           datasource="Vector"
         ),
         prototype=prototype(
           datasource=GRangesList()
         ),
         contains="AnnotParam")

##==========================
## RnaSeqParam
##==========================

#' Class "RnaSeqParam"
#'
#' A class holding all the necessary parameters to process a bam file issued
#' from an RNA-Seq experiment together with the related annotation to compute
#' a count-table using the \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}.
#' The precision slot is used to determine the count unit:
#' \itemize{
#' \item{reads}{default. The standard \code{\link[GenomicAlignments]{summarizeOverlaps-methods}} function is used to extract the read counts}
#' \item{bp}{The \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization functions}} are used to extract the read covered bp counts}
#' }
#' @name RnaSeqParam class
#' @rdname easyRNASeq-RnaSeqParam-class
#' @aliases RnaSeqParam-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the
#' form \code{new("RnaSeqParam", ...)} or using the RnaSeqParam constructor.
#' @author Nicolas Delhomme
#' @seealso \itemize{
#' \item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam]{RnaSeqParam constructor}}
#' \item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam accessors}}
#' \item \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}
#' \item \code{\linkS4class{AnnotParam}}
#' \item \code{\link[easyRNASeq:easyRNASeq-AnnotParam]{AnnotParam constructor}}
#' \item \code{\linkS4class{BamParam}}
#' \item \code{\link[easyRNASeq:easyRNASeq-BamParam]{BamParam constructor}}
#' \item \code{\link[GenomicAlignments]{summarizeOverlaps-methods}}
#' \item \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization functions}}
#' }
#' @keywords classes
#' @examples
#'
#' showClass("RnaSeqParam")
#'
#' @exportClass RnaSeqParam
setClass("RnaSeqParam",
         representation=representation(
           annotParam="AnnotParam",
           bamParam="BamParam",
           countBy="character",
           precision="character"),
         prototype=prototype(
           countBy="transcripts",
           precision="read"
         ))

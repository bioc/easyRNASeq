## Not sure we need a prototype later on

## TODO should we have one readAnnotation slot
## containing all the genomicAnnoation, geneModel, and readIslands info?
## and have it more dynamic?
##' Class "RNAseq"
##' 
##' A class holding all the necessary information and annotation to summarize
##' couts (number of reads) per features (i.e. exons or transcripts or genes)
##' for RNA-Seq experiments.
##' 
##' 
##' @name RNAseq class
##' @rdname easyRNASeq-class
##' @aliases RNAseq-class RNAseq chrSize,RNAseq-method chrSize<-,RNAseq-method
##' easyRNASeq,RNAseq-method exonCounts,RNAseq-method
##' featureCounts,RNAseq-method fetchAnnotation,RNAseq-method
##' fetchCoverage,RNAseq-method fileName,RNAseq-method fileName<-,RNAseq-method
##' findIslands,RNAseq-method geneCounts,RNAseq-method geneModel,RNAseq-method
##' geneModel<-,RNAseq-method genomicAnnotation,RNAseq-method
##' genomicAnnotation<-,RNAseq-method islandCounts,RNAseq-method
##' librarySize,RNAseq-method librarySize<-,RNAseq-method
##' organismName,RNAseq-method organismName<-,RNAseq-method print,RNAseq-method
##' readCounts,RNAseq-method readCounts<-,RNAseq-method
##' readCoverage,RNAseq-method readCoverage<-,RNAseq-method
##' readIslands,RNAseq-method readIslands<-,RNAseq-method
##' readLength,RNAseq-method readLength<-,RNAseq-method RPKM,RNAseq-method
##' show,RNAseq-method transcriptCounts,RNAseq-method
##' @docType class
##' @section Objects from the Class: Objects can be created by calls of the
##' form \code{new("RNAseq", ...)}.
##' @author Nicolas Delhomme
##' @seealso \itemize{ \item\code{\linkS4class{RangedData}}
##' \item\code{\linkS4class{RleList}}
##' \item\code{\link[easyRNASeq:easyRNASeq]{easyRNASeq function}}
##' \item\code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}}
##' \item\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq
##' annotation methods}}
##' \item\code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq
##' correction methods}}
##' \item\code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq
##' coverage methods}}
##' \item\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq
##' summarization methods}} \item\code{\link[easyRNASeq:print-methods]{print}}
##' }
##' @keywords classes
##' @examples
##' 
##' showClass("RNAseq")
##'
##' @exportClass RNAseq
setClass(
         Class="RNAseq",
         representation=representation(
           chrSize="list",
           fileName="character",
           geneModel="RangedData",
           genomicAnnotation="Vector",
           librarySize="numeric",
           organismName="character",
           readCounts="list",
           readCoverage="RleList",
           readIslands="RangedData",
           readLength="integer"
           ),
         prototype=prototype(
           chrSize=list(),
           fileName=character(0),
           geneModel=RangedData(),           
           genomicAnnotation=RangedData(),
           librarySize=numeric(0),
           organismName=character(1),
           readCounts=list(),
           readCoverage=RleList(),           
           readIslands=RangedData(),
           readLength=integer(1)
           )
##         ,
##          validity=function(obj){
##          }
         )


\name{RNAseq-class}
\Rdversion{1.1}
\docType{class}
\alias{RNAseq-class}
\alias{chrSize,RNAseq-method}
\alias{chrSize<-,RNAseq-method}
\alias{easyRNASeq,RNAseq-method}
\alias{exonCounts,RNAseq-method}
\alias{featureCounts,RNAseq-method}
\alias{fetchAnnotation,RNAseq-method}
\alias{fetchCoverage,RNAseq-method}
\alias{fileName,RNAseq-method}
\alias{fileName<-,RNAseq-method}
\alias{findIslands,RNAseq-method}
\alias{geneCounts,RNAseq-method}
\alias{geneModel,RNAseq-method}
\alias{geneModel<-,RNAseq-method}
\alias{genomicAnnotation,RNAseq-method}
\alias{genomicAnnotation<-,RNAseq-method}
\alias{islandCounts,RNAseq-method}
\alias{librarySize,RNAseq-method}
\alias{librarySize<-,RNAseq-method}
\alias{organismName,RNAseq-method}
\alias{organismName<-,RNAseq-method}
\alias{print,RNAseq-method}
\alias{readCounts,RNAseq-method}
\alias{readCounts<-,RNAseq-method}
\alias{readCoverage,RNAseq-method}
\alias{readCoverage<-,RNAseq-method}
\alias{readIslands,RNAseq-method}
\alias{readIslands<-,RNAseq-method}
\alias{readLength,RNAseq-method}
\alias{readLength<-,RNAseq-method}
\alias{RPKM,RNAseq-method}
\alias{show,RNAseq-method}
\alias{transcriptCounts,RNAseq-method}

\title{Class "RNAseq"}
\description{
	A class holding all the necessary information and annotation to summarize couts (number of reads) per features (i.e. exons or transcripts or genes) for RNA-Seq experiments.
}
\section{Objects from the Class}{
	 Objects can be created by calls of the form \code{new("RNAseq", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{chrSize}:}{A named \code{"list"} containing the chromosomes' size}
    \item{\code{fileName}:}{Name of the last processed file}
    \item{\code{geneModel}:}{A \code{\linkS4class{RangedData}} containing gene models}
    \item{\code{genomicAnnotation}:}{Object of class \code{\linkS4class{RangedData}} or \code{\linkS4class{GRangesList}} containing genomic annotations, i.e. exon coordinates, xrefs, \emph{etc.}}
    \item{\code{librarySize}:}{A named numeric list describing the library size of every sample.}
    \item{\code{organismName}:}{A \code{"character"} describing the organism from which the sample was derived.}
    \item{\code{readCounts}:}{A \code{"list"} of \code{"list"} containing summarized counts per summarization features, i.e. per exon, or transcript or gene}
    \item{\code{readCoverage}:}{A \code{\linkS4class{RleList}} containing the provided sample coverage}
    \item{\code{readIslands}:}{A \code{\linkS4class{RangedData}} containing read islands}
    \item{\code{readLength}:}{An \code{"integer"} specifying the read length }
  }
}
\section{Methods}{
  \describe{
    \item{chrSize}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{chrSize<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{easyRNASeq}{\code{signature(normalizeMethod = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq]{easyRNASeq function}} }
    \item{exonCounts}{\code{signature(obj = "RNAseq")}: \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}} }
    \item{fetchAnnotation}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq annotation methods}} }
    \item{fetchCoverage}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq methods}} }
    \item{fileName}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{fileName<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{geneCounts}{\code{signature(obj = "RNAseq")}: \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}} }
    \item{geneModel}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{geneModel<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{genomicAnnotation}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{genomicAnnotation<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{islandCounts}{\code{signature(obj = "RNAseq")}:\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}}}
    \item{organismName}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{organismName<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{print}{\code{signature(x = "RNAseq")}: see \code{\link[easyRNASeq:print-methods]{print}} }
    \item{readCounts}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{readCounts<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{readCoverage}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{readCoverage<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{readLength}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{readLength<-}{\code{signature(obj = "RNAseq")}: see \code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}} }
    \item{RPKM}{\code{signature(obj = "RNAseq")}: \code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq correction methods}} }
    \item{show}{\code{signature(object = "RNAseq")}: see \code{\link[easyRNASeq:show-methods]{show}} }
    \item{transcriptCounts}{\code{signature(obj = "RNAseq")}: \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}} }
    }
}
\author{
	Nicolas Delhomme
}

\seealso{
	\itemize{
	\item{\code{\linkS4class{RangedData}}}
	\item{\code{\linkS4class{RleList}}}
	\item{\code{\link[easyRNASeq:easyRNASeq]{easyRNASeq function}}}
	\item{\code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq accessors}}}
	\item{\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq annotation methods}}}
	\item{\code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq correction methods}}}
	\item{\code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq coverage methods}}}
	\item{\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}}}
	\item{\code{\link[easyRNASeq:print-methods]{print}}}
	}
}
\examples{
showClass("RNAseq")
}
\keyword{classes}

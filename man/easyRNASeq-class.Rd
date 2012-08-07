\docType{class}
\name{RNAseq class}
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
\alias{RNAseq}
\alias{RNAseq class}
\alias{RNAseq-class}
\alias{RPKM,RNAseq-method}
\alias{show,RNAseq-method}
\alias{transcriptCounts,RNAseq-method}
\title{Class "RNAseq"}
\description{
  A class holding all the necessary information and
  annotation to summarize couts (number of reads) per
  features (i.e. exons or transcripts or genes) for RNA-Seq
  experiments.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("RNAseq", ...)}.
}
\examples{
showClass("RNAseq")
}
\author{
  Nicolas Delhomme
}
\seealso{
  \itemize{ \item\code{\linkS4class{RangedData}}
  \item\code{\linkS4class{RleList}}
  \item\code{\link[easyRNASeq:easyRNASeq]{easyRNASeq
  function}}
  \item\code{\link[easyRNASeq:easyRNASeq-accessors]{RNAseq
  accessors}}
  \item\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq
  annotation methods}}
  \item\code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq
  correction methods}}
  \item\code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq
  coverage methods}}
  \item\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq
  summarization methods}}
  \item\code{\link[easyRNASeq:print-methods]{print}} }
}
\keyword{classes}


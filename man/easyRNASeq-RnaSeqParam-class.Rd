% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{RnaSeqParam class}
\alias{RnaSeqParam class}
\alias{RnaSeqParam-class}
\title{Class "RnaSeqParam"}
\description{
A class holding all the necessary parameters to process a bam file issued
from an RNA-Seq experiment together with the related annotation to compute
a count-table using the \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}.
The precision slot is used to determine the count unit:
\itemize{
\item{reads}{default. The standard \code{\link[GenomicAlignments]{summarizeOverlaps-methods}} function is used to extract the read counts}
\item{bp}{The \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization functions}} are used to extract the read covered bp counts}
}
}
\section{Objects from the Class}{
 Objects can be created by calls of the
form \code{new("RnaSeqParam", ...)} or using the RnaSeqParam constructor.
}

\examples{

showClass("RnaSeqParam")

}
\seealso{
\itemize{
\item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam]{RnaSeqParam constructor}}
\item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam accessors}}
\item \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}
\item \code{\linkS4class{AnnotParam}}
\item \code{\link[easyRNASeq:easyRNASeq-AnnotParam]{AnnotParam constructor}}
\item \code{\linkS4class{BamParam}}
\item \code{\link[easyRNASeq:easyRNASeq-BamParam]{BamParam constructor}}
\item \code{\link[GenomicAlignments]{summarizeOverlaps-methods}}
\item \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization functions}}
}
}
\author{
Nicolas Delhomme
}
\keyword{classes}

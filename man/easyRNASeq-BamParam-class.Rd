% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{BamParam class}
\alias{BamParam class}
\alias{BamParam-class}
\title{Class "BamParam"}
\description{
A class describing the parameters of a bam file issued
from an RNA-Seq experiment.
}
\section{Objects from the Class}{
 Objects can be created by calls of the
form \code{new("BamParam", ...)} or using the BamParam constructor.
}

\section{Slots from the Class}{
 The \code{\linkS4class{BamParam}} class has
the following slots:
\itemize{
\item paired
\item stranded
\item strandProtocol
\item yieldSize
}
all of which can be accessed using the accordingly names accessor.
}

\examples{

showClass("BamParam")

}
\seealso{
\itemize{
\item \code{\link[easyRNASeq:easyRNASeq-BamParam-accessors]{BamParam accessors}}
\item \code{\linkS4class{RnaSeqParam}}
\item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam]{RnaSeqParam constructor}}
\item \code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam accessors}}
\item \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq function}}
\item \code{\linkS4class{AnnotParam}}
\item \code{\link[easyRNASeq:easyRNASeq-AnnotParam]{AnnotParam constructor}}
}
}
\author{
Nicolas Delhomme
}
\keyword{classes}

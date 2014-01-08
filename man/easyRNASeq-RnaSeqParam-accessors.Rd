\docType{methods}
\name{easyRNASeq RnaSeqParam accessors}
\alias{RnaSeqParam-accessors}
\alias{annotParam}
\alias{annotParam,RnaSeqParam-method}
\alias{bamParam}
\alias{bamParam,RnaSeqParam-method}
\alias{countBy}
\alias{countBy,RnaSeqParam-method}
\alias{datasource,RnaSeqParam-method}
\alias{easyRNASeq RnaSeqParam accessors}
\alias{paired,RnaSeqParam-method}
\alias{precision}
\alias{precision,RnaSeqParam-method}
\alias{stranded,RnaSeqParam-method}
\alias{yieldSize,RnaSeqParam-method}
\title{Accessors for RnaSeqParam class}
\usage{
\S4method{yieldSize}{RnaSeqParam}(object)
}
\arguments{
  \item{object}{An object derived from class
  \code{RnaSeqParam}.}
}
\value{
The value of the corresponding slot.
}
\description{
These functions and generics define `accessors` (to get and
set values) for \code{\linkS4class{RnaSeqParam}} objects
within the \pkg{easyRNASeq} package. Implemented are:
\itemize{ \item annotParam \item bamParam \item countBy
\item datasource \item paired \item precision \item
stranded \item yieldSize }
}
\examples{
## create the necessary AnnotParam
annotParam <- AnnotParam(
                datasource=system.file(
                                "extdata",
                                "Dmel-mRNA-exon-r5.52.gff3",
                                package="RnaSeqTutorial"))

## create the RnaSeqParam
rsp <- RnaSeqParam(annotParam=annotParam)
## get the yieldSize Parameter
ysize <-yieldSize(rsp)
}
\author{
Nicolas Delhomme
}
\seealso{
\itemize{ \item The
\code{\link[easyRNASeq:easyRNASeq-AnnotParam-class]{AnnotParam}}
class \item The
\code{\link[easyRNASeq:easyRNASeq-BamParam-class]{BamParam}}
class \item The
\code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-class]{RnaSeqParam}}
class } The
\code{\link[easyRNASeq:easyRNASeq-BamParam-accessors]{BamParam
yieldSize}} accessor
}
\keyword{manip}


\docType{methods}
\name{easyRNASeq BamParam accessors}
\alias{BamParam-accessors}
\alias{easyRNASeq BamParam accessors}
\alias{paired}
\alias{paired,BamParam-method}
\alias{stranded}
\alias{stranded,BamParam-method}
\alias{yieldSize,BamParam-method}
\title{Accessors for BamParam class}
\usage{
\S4method{yieldSize}{BamParam}(object)
}
\arguments{
  \item{object}{An object derived from class
  \code{BamParam}.}
}
\value{
The value of the corresponding slot.
}
\description{
These functions and generics define `accessors` (to get and
set values) for \code{\linkS4class{BamParam}} objects
within the \pkg{easyRNASeq} package.
}
\examples{
bp <- BamParam()
## get the yieldSize Parameter
ysize <-yieldSize(bp)
}
\author{
Nicolas Delhomme
}
\seealso{
The
\code{\link[easyRNASeq:easyRNASeq-BamParam-class]{BamParam}}
class The
\code{\link[easyRNASeq:easyRNASeq-RnaSeqParam-accessors]{RnaSeqParam
yieldSize}} accessor
}
\keyword{manip}


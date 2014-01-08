\docType{methods}
\name{simpleRNASeq,character-method}
\alias{simpleRNASeq}
\alias{simpleRNASeq,character-method}
\title{simpleRNASeq method}
\usage{
\S4method{simpleRNASeq}{character}(bamFiles = BamFileList(),
  param = RnaSeqParam(), nnodes = 1, verbose = FALSE)
}
\arguments{
  \item{bamFiles}{a \code{\linkS4class{BamFileList}}
  object}

  \item{nnodes}{The number of CPU cores to use in parallel}

  \item{param}{RnaSeqParam a
  \code{\linkS4class{RnaSeqParam}} object that describres
  the RNA-Seq experimental setup.}

  \item{verbose}{a logical to be report progress or not.}
}
\value{
returns a \code{\linkS4class{SummarizedExperiment}} object.
}
\description{
This function is a wrapper around the more low level
functionalities of the package. It is the simplest way to
get a \code{\linkS4class{SummarizedExperiment}} object from
a set of bam files.
\code{\linkS4class{SummarizedExperiment}} are containers
meant to hold any Next-Generation Sequencing experiment
results and metadata. The simpleRNASeq method replaces the
\code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}}
function to simplify the usability. It does the following:
\itemize{ \item use
\code{\link[GenomicRanges:GRanges-class]{GenomicRanges}}
for loading/pre-processing the data. \item fetch
\code{\link[easyRNASeq:fetchAnnotation]{annotations}}
depending on the selected parameters \item extract the
reads \code{\link[easyRNASeq:fetchCoverage]{coverage}} from
the provided file(s) \item
\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarizes}}
the read counts according to the selected summarization
\item returns a \code{\linkS4class{SummarizedExperiment}}
object. }
}
\examples{
\dontrun{
  TODO
}
}
\author{
Nicolas Delhomme
}
\seealso{
\code{\linkS4class{SummarizedExperiment}}
\code{\linkS4class{RnaSeqParam}}
}
\keyword{methods}


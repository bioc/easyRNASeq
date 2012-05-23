\name{easyRNASeq summarization internal methods}
\alias{.bestExonSummarization}
\alias{.geneModelSummarization}
\alias{easyRNASeq summarization internal methods}
\title{Summarization of reads count per gene}
\usage{
  .bestExonSummarization(obj)
}
\arguments{
  \item{obj}{An object derived from class
  \code{\linkS4class{RNAseq}}}
}
\value{
  A vector containing read counts
}
\description{
  Called by the method \code{geneCounts}. Having performed
  the \code{exonCounts} is a pre-requisite.  \itemize{
  \item.bestExonSummarization Identify the exon showing the
  highest coverage. \item.geneModelSummarization Sum the
  coverage values of the synthetic exons constituting a
  gene model.  }
}
\author{
  Nicolas Delhomme
}
\keyword{internal}


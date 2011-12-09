\name{easyRNASeq summarization internal methods}

\alias{.bestExonSummarization}
\alias{.geneModelSummarization}

\title{Summarization of reads count per gene}

\description{
  Called by the method \code{geneCounts}. Having performed the \code{exonCounts} is a pre-requisite.
  \itemize{  
  \item{.bestExonSummarization} Identify the exon showing the highest coverage.
  \item{.geneModelSummarization} Sum the coverage values of the synthetic exons constituting a gene model.
  }
}

\usage{
  .bestExonSummarization(obj)
  .geneModelSummarization(obj)
}

\arguments{
  \item{obj}{An object derived from class \code{\linkS4class{RNAseq}}}
}

\value{
  A vector containing read counts
}

\author{Nicolas Delhomme}

\keyword{internal}


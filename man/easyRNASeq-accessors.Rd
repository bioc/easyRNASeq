\name{easyRNASeq accessors}
\alias{accessors}
\alias{chrSize}
\alias{chrSize<-}
\alias{easyRNASeq accessors}
\alias{fileName}
\alias{fileName<-}
\alias{geneModel}
\alias{geneModel<-}
\alias{genomicAnnotation}
\alias{genomicAnnotation<-}
\alias{librarySize}
\alias{librarySize<-}
\alias{organismName}
\alias{organismName<-}
\alias{readCounts}
\alias{readCounts<-}
\alias{readCoverage}
\alias{readCoverage<-}
\alias{readIslands}
\alias{readIslands<-}
\alias{readLength}
\alias{readLength<-}
\title{Accessors for RNAseq class}
\arguments{
  \item{obj}{An object derived from class \code{RNAseq}.}

  \item{count}{The type of count you want to access,
  'genes','features','exons','transcripts' or 'islands'}

  \item{summarization}{If count is set to genes, precise
  the type of summarization, 'bestExons' or 'geneModels'}

  \item{unique}{For the 'exons' count only. Should the
  counts returned be unique for their identifier (i.e. the
  matrix row names)?}

  \item{value}{The replacement value.}
}
\value{
  Usually, the value of the corresponding slot, or other
  simple content described on the help page of
  \code{easyRNASeq}.
}
\description{
  These functions and generics define `accessors` (to get
  and set values) for objects in the \pkg{easyRNASeq}
  package.
}
\examples{
rnaSeq<-new("RNAseq")
##set organisme name of an RNAseq object
organismName(rnaSeq) <- "Dmelanogaster"
##get organisme name of an RNAseq object
orgName<-organismName(rnaSeq)
}
\author{
  Nicolas Delhomme
}
\keyword{manip}


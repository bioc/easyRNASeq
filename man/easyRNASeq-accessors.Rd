\name{accessors}

\alias{accessors}
%chrSize
\alias{chrSize<-}
\alias{chrSize}
%fileName
\alias{fileName<-}
\alias{fileName}
%geneModel
\alias{geneModel<-}
\alias{geneModel}
%genomicAnnotation
\alias{genomicAnnotation<-}
\alias{genomicAnnotation}
%librarySize
\alias{librarySize<-}
\alias{librarySize}
%organismName
\alias{organismName<-}
\alias{organismName}
%readCounts
\alias{readCounts<-}
\alias{readCounts}
%readCoverage
\alias{readCoverage<-}
\alias{readCoverage}
%readIslands
\alias{readIslands<-}
\alias{readIslands}
%readLength
\alias{readLength<-}
\alias{readLength}

\title{Accessors for RNAseq class}

\description{

  These functions and generics define `accessors` (to get and set
  values) for objects in the \pkg{easyRNASeq} package.

}

\arguments{

  \item{obj}{An object derived from class \code{RNAseq}.}
  The default is to extract or set the contents of a slot of the corresponding name.
  \item{count}{The type of count you want to access, 'genes','features','exons','transcripts' or 'islands'}
  \item{summarization}{If count is set to genes, precise the type of summarization, 'bestExons' or 'geneModels'}
  \item{unique}{For the 'exons' count only. Should the counts returned be unique for their identifier (i.e. the matrix row names)?}
  \item{value}{The replacement value.}
}

\value{

  Usually, the value of the corresponding slot, or other simple content
  described on the help page of \code{easyRNASeq}.

}

\usage{

##chrSize
chrSize(obj)<-value
chrSize(obj)

##fileName
fileName(obj)<-value
fileName(obj)

##geneModel
geneModel(obj)<-value
geneModel(obj)

##genomicAnnotation
genomicAnnotation(obj)<-value
genomicAnnotation(obj)

##librarySize
librarySize(obj)<-value
librarySize(obj)

##organismName
organismName(obj)<-value
organismName(obj)

##readCounts
readCounts(obj)<-value
readCounts(obj,count=c("exons","features","genes","islands","transcripts"),summarization=c("bestExons","geneModels"),unique=FALSE)

##readCoverage
readCoverage(obj)<-value
readCoverage(obj)

##readIslands
readIslands(obj)<-value
readIslands(obj)

##readLength
readLength(obj)<-value
readLength(obj)
}

\examples{
rnaSeq<-new("RNAseq")
##set organisme name of an RNAseq object
organismName(rnaSeq) <- "Dmelanogaster"
##get organisme name of an RNAseq object
orgName<-organismName(rnaSeq)
}

\author{Nicolas Delhomme}
\keyword{manip}


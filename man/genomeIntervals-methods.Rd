\docType{methods}
\name{genomeIntervals additional methods}
\alias{as}
\alias{as,Genome_intervals,GRangesList-method}
\alias{as,Genome_intervals,RangedData-method}
\alias{as,Genome_intervals-method}
\alias{coerce,Genome_intervals,GRangesList-method}
\alias{coerce,Genome_intervals,RangedData-method}
\alias{genomeIntervals additional methods}
\alias{type}
\alias{type,Genome_intervals-method}
\title{Extension for the genomeIntervals package}
\usage{
\S4method{type}{Genome_intervals}(x)
\S4method{as}{Genome_intervals}(from,to)
}
\arguments{
  \item{from}{An object of class
  \code{\linkS4class{Genome_intervals}}}

  \item{to}{a character string; either RangedData or
  GRangesList}

  \item{x}{An object of class
  \code{\linkS4class{Genome_intervals}}}
}
\value{
\describe{ \item{coerce}{ A \code{\linkS4class{RangedData}}
or \code{\linkS4class{GRangesList}} containing the result
of the coercion.  } \item{type}{ The content of the type
column, usually a factor or a character vector } }
}
\description{
\describe{ \item{coerce}{ This method extends the
genomeIntervals package by offering the functionality to
coerce a
\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals
object}} into a
\code{\link[IRanges:RangedData-class]{RangedData object}}
or \code{\linkS4class{GRangesList}} object.} \item{type}{
Another way to access the content of the gff type column.
} }
}
\examples{
\dontrun{
annot<-readGff3(system.file("extdata","annot.gff",package="RnaSeqTutorial")
gAnnot<-as(from=annot,to="RangedData")
type(annot)
}
}
\author{
Nicolas Delhomme
}
\seealso{
\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals
object}} \code{\link[genomeIntervals:readGff3]{readGff3
function}}
}


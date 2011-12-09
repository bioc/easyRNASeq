\name{genomeIntervals methods}

\alias{coerce,Genome_intervals,RangedData-method}

\title{ Extension for the genomeIntervals package}
\description{
	This method extends the genomeIntervals package by offering the functionality to coerce a 
	\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals object}} into
	a \code{\link[IRanges:RangedData-class]{RangedData object}} one.
}

\arguments{
  \item{from}{An object of class \code{\linkS4class{Genome_intervals}}}
}

\value{
	A \code{\linkS4class{RangedData}} containing the result of the coercion.
}

\examples{
	\dontrun{
	annot<-readGff3(system.file("extdata","annot.gff",package="RnaSeqTutorial")
	gAnnot<-as(annot,"RangedData")
	}
}

\author{Nicolas Delhomme}

\seealso{
	\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals object}}
	\code{\link[genomeIntervals:readGff3]{readGff3 function}}
}

\keyword{methods}


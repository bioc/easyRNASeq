\name{IRanges methods}

\alias{ranges}
\alias{range,RNAseq-method}

\title{Extension of the IRanges package}
\description{
	Return the ranges of the genomic annotation.
}

\arguments{
  \item{obj}{An object of the \code{\linkS4class{RNAseq}} class}
}

\value{
	An \code{\linkS4class{IRanges}} object.
}

\details{
	It retrieves the object stored in the genomicAnnotation slot of the RNAseq object and apply the \code{ranges} function on it. The object retrieved can be of the \code{\linkS4class{RangedData}} or \code{\linkS4class{GRangesList}} class.
}

\examples{
	\dontrun{
	library("RnaSeqTutorial")
	obj <- new('RNAseq',
		organismName="Dmelanogaster",
		readLength=36L,
		chrSize=as.list(seqlengths(Dmelanogaster))
		)

	obj <- fetchAnnotation(obj,
				method="gff",
                                filename=system.file(
						"extdata",
						"annot.gff",
						package="RnaSeqTutorial"))
	ranges(obj)
	}
}

\author{Nicolas Delhomme}

\keyword{methods}


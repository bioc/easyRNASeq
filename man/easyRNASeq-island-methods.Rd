\name{easyRNASeq island methods}
\alias{easyRNASeq island methods}
\alias{findIslands}
\title{Identify expressed regions de-novo}
\arguments{
  \item{obj}{An object of class \code{RNAseq}}

  \item{max.gap}{Maximum gap between two peaks to build an
  island}

  \item{min.cov}{Minimum coverage for an island to be
  returned}

  \item{min.length}{Minimum size of an island to be
  returned}

  \item{plot}{If TRUE, draw plots of coverage distribution.
  Help the user to select an appropriate value for the
  minimum coverage.}

  \item{\dots}{See details}
}
\value{
  An \code{RNAseq} object with the readIsland slot set with
  a RangedData containing the selected islands and the
  readCount slot actualized with a list containing the
  count table per island.
}
\description{
  Process the coverage to locate regions with a minimum
  coverage (min.cov). If regions are separated by a gap
  shorter than a maximum length (max.gap), they are
  unified. Only islands longer than min.length are
  returned.  These functions are now outdated and would
  need to be actualized.
}
\details{
  \dots{} are for providing additional options to the
  \code{\link[graphics:hist]{hist}} plot function.
}
\examples{
\dontrun{
	## NOTE that this function might need to be actualized
	obj <- new('RNAseq',
		organismName="Dmelanogaster",
		readLength=36L,
		chrSize=as.list(seqlengths(Dmelanogaster))
		)

	obj <- fetchCoverage(
			obj,
			format="bam",
                        filename=system.file(
				"extdata",
				"ACACTG.bam",
                            	package="RnaSeqTutorial")
			)

	obj <- findIslands(
			obj,
			max.gap=10L,
			min.cov=10L,
			min.length=200L)
	}
}
\author{
  Nicolas Delhomme
}
\keyword{connection}
\keyword{data}
\keyword{methods}


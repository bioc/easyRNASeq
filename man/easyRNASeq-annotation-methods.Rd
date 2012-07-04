\name{easyRNASeq annotation methods}
\alias{easyRNASeq annotation methods}
\alias{fetchAnnotation}
\alias{knownOrganisms}
\alias{knownOrganisms,missing-method}
\title{Fetch genic annotation from a gff/gtf file or using biomaRt}
\arguments{
  \item{obj}{An object of class \code{RNAseq}}

  \item{method}{one of biomaRt, gff, gtf}

  \item{filename}{If the method is gff or gtf, the actual
  gtf, gff filename}

  \item{ignoreWarnings}{set to TRUE (bad idea! they have a
  good reason to be there) if you do not want warning
  messages.}

  \item{\dots}{See details}
}
\value{
  A \code{\linkS4class{RangedData}} containing the fetched
  annotations.

  A vector containing the known organisms
}
\description{
  The annotation can be retrieved in two ways \itemize{
  \item{biomaRt}{Use biomaRt and Ensembl to get organism
  specific annotation.} \item{gff/gtf}{Use a gff or gtf
  local annotation file.}} When using \pkg{biomaRt}, it is
  important that the \code{organismName} slot of the
  \code{\linkS4class{RNAseq}} object is set the prefix of
  one of the value available using the \pkg{biomaRt}
  \code{\link[biomaRt:listDatasets]{listDatasets}}
  function, e.g. "Dmelanogaster".  When reading from a
  gff/gtf file, a version 3 formatted gff (gtf are modified
  gff3 from Ensembl) is expected. The function
  \pkg{genomeIntervals}
  \code{\link[genomeIntervals:readGff3]{readGff3}} is used
  to read the data in. Another annotation caveat is the
  reference names, \emph{i.e.} the chromosome/scaffold
  names used in the alignment files and those fetched when
  retrieving the genic annotation might differ.
  \pkg{easyRNASeq} tries to be clever in this case and
  guess the correspondance. However, it is not always
  obvious. Organisms were this has been checked can be
  listed with the \emph{knownOrganisms} function.
}
\details{
  \dots{} are for additional arguments, passed to the
  \pkg{biomaRt} \code{\link[biomaRt:getBM]{getBM}} function
  or to the
  \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
  internal function that takes an optional arguments:
  annotation.type that default to "exon". This is used to
  select the proper rows of the gff or gtf file.
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
	}
knownOrganisms()
}
\author{
  Nicolas Delhomme

  Nicolas Delhomme
}
\keyword{connection}
\keyword{data}
\keyword{methods}


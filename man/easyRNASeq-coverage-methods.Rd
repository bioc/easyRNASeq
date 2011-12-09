\name{easyRNASeq coverage methods}
\alias{fetchCoverage}

\title{Compute the coverage from a Short Read Alignment file}

\description{
	\itemize{
	\item{fetchCoverage:}
	Computes the genomic reads' coverage from a read file in bam format or any format supported by \pkg{ShortRead}.
	}
}

\usage{
fetchCoverage(	obj,
		format=c("aln","bam"),
		filename=character(1),
		filter=srFilter(),
		type="SolexaExport",
		chr.sel=c(),
		isUnmappedQuery=FALSE,
		what=c("rname","pos","qwidth"),
		validity.check=TRUE,
		chr.map=data.frame(),
		ignoreWarnings = FALSE, 
		gapped = TRUE,...)
}

\arguments{
	\item{obj}{An \code{\linkS4class{RNAseq}} object}
	\item{chr.map}{A data.frame describing the mapping of original chromosome names towards wished chromosome names. See details.}
	\item{chr.sel}{A vector of chromosome names to subset the final results.}
	\item{filename}{The full path of the file to use}
	\item{filter}{The filter to be applied when loading the data using the "aln" format}
	\item{format}{The format of the reads, one of "aln","bam". If not "bam", all the types supported by the ShortRead package are supported too.}
	\item{gapped}{Is the bam file provided containing gapped alignments?}
	\item{ignoreWarnings}{set to TRUE (bad idea! they have a good reason to be there) if you do not want warning messages.}
	\item{isUnmappedQuery}{additional argument for scanBamFlag \pkg{Rsamtools}}
	\item{type}{The type of data when using the "aln" format. See the \pkg{ShortRead} package.}
	\item{validity.check}{Shall UCSC chromosome name convention be enforced}
	\item{what}{additional argument for ScanBamParam \pkg{Rsamtools}}
	\item{\dots}{additional arguments. See details}
}

\value{
	An \code{\linkS4class{RNAseq}} object. The slot readCoverage contains a SimpleRleList object representing a list of coverage vectors, one per chromosome.
}

\details{
  \dots for fetchCoverage:
  Can be used for readAligned method from package \pkg{ShortRead} or for scanBamFlag method from package \pkg{Rsamtools}.
}

\author{Nicolas Delhomme}

\seealso{
	\code{\linkS4class{Rle}}
	\code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
}

\examples{
	\dontrun{
	library("RnaSeqTutorial")
	library(BSgenome.Dmelanogaster.UCSC.dm3)

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
	}
}

\keyword{methods}


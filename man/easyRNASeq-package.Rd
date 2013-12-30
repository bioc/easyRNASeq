\docType{package}
\name{easyRNASeq package}
\alias{easyRNASeq package}
\alias{easyRNASeq package-package}
\alias{easyRNASeq-package}
\title{Count summarization and normalization pipeline for Next Generation Sequencing data.}
\description{
Offers functionalities to summarize read counts per feature
of interest, e.g. exons, transcripts, genes, etc. Offers
functionalities to normalize the summarized counts using
3rd party packages like
\code{\link[DESeq:newCountDataSet]{DESeq}} or
\code{\link[edgeR:DGEList]{edgeR}}.
}
\details{
\tabular{ll}{ Package: \tab easyRNASeq\cr Type: \tab
Package\cr Version: \tab 1.8.4\cr Date: \tab 2013-12-29\cr
License: \tab Artistic-2.0\cr LazyLoad: \tab yes\cr
Depends: \tab methods, parallel, biomaRt, edgeR, DESeq,
genomeIntervals, LSD, Rsamtools, ShortRead,
RnaSeqTutorial\cr Suggests: \tab
BSgenome.Dmelanogaster.UCSC.dm3 }
}
\section{Methods}{
  The main function
  \code{\link[easyRNASeq:easyRNASeq]{easyRNASeq}} will
  summarize the counts per feature of interest, for as many
  samples as provided and will return a count matrix (N*M)
  where N are the features and M the samples. This data can
  be corrected to \pkg{RPKM} in which case a matrix of
  corrected value is returned instead, with the same
  dimensions. Alternatively a
  \code{\linkS4class{SummarizedExperiment}} can be returned
  and this is expected to be the default in the upcoming
  version of easyRNASeq (as of 1.5.x). If the necessary
  sample information are provided, the data can be
  normalized using either
  \code{\link[DESeq:newCountDataSet]{DESeq}} or
  \code{\link[edgeR:DGEList]{edgeR}} and the corresponding
  package object returned. For more insider details, and
  step by step functions, see: \tabular{ll}{
  \code{\link[easyRNASeq:ShortRead-methods]{ShortRead
  methods}} for pre-processing the data.
  \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq
  annotation methods}} for getting the annotation.
  \code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq
  coverage methods}} for computing the coverage from a
  Short Read Alignment file.
  \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq
  summarization methods}} for summarizing the data.
  \code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq
  correction methods}} for correcting the data (i.e.
  generating RPKM).
  \code{\link[easyRNASeq:edgeR-methods]{edgeR methods}} for
  post-processing the data.
  \code{\link[easyRNASeq:DESeq-methods]{DESeq methods}} for
  post-processing the data.  }
}
\examples{
\dontrun{
	library("RnaSeqTutorial")
	library(BSgenome.Dmelanogaster.UCSC.dm3)

	## creating a count table from 4 bam files
	count.table <- easyRNASeq(filesDirectory=
		    			system.file(
					"extdata",
					package="RnaSeqTutorial"),
					pattern="[A,C,T,G]{6}\\.bam$",
				format="bam",
				readLength=36L,
				organism="Dmelanogaster",
				chr.sizes=as.list(seqlengths(Dmelanogaster)),
				annotationMethod="rda",
				annotationFile=system.file(
				                            "data",
							    "gAnnot.rda",
							    package="RnaSeqTutorial"),
				count="exons")

	}
}
\author{
Nicolas Delhomme, Bastian Schiffthaler, Ismael Padioleau
}
\seealso{
The class RNAseq specification: \code{\linkS4class{RNAseq}}

The default output class specification:
\code{\linkS4class{SummarizedExperiment}}

The imported packages:
\code{\link[biomaRt:useMart]{biomaRt}}
\code{\link[edgeR:DGEList]{edgeR}}
\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals}}
\code{\link[Biostrings:XString-class]{Biostrings}}
\code{\link[BSgenome:BSgenome-class]{BSgenome}}
\code{\link[DESeq:newCountDataSet]{DESeq}}
\code{\link[GenomicRanges:GRanges-class]{GenomicRanges}}
\code{\link[IRanges:IRanges-constructor]{IRanges}}
\code{\link[Rsamtools:scanBam]{Rsamtools}}
\code{\link[ShortRead:readAligned]{ShortRead}}

The suggested packages:
\code{\link[parallel:makeCluster]{parallel}}
\code{\link[GenomicFeatures:TranscriptDb-class]{GenomicFeatures}}
}
\keyword{package}


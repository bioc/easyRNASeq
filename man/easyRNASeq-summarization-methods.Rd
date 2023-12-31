% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyRNASeq-counts.R
\name{easyRNASeq summarization methods}
\alias{easyRNASeq summarization methods}
\alias{exonCounts}
\alias{exonCounts,RNAseq-method}
\alias{featureCounts}
\alias{featureCounts,RNAseq-method}
\alias{geneCounts}
\alias{geneCounts,RNAseq-method}
\alias{islandCounts}
\alias{islandCounts,RNAseq-method}
\alias{transcriptCounts}
\alias{transcriptCounts,RNAseq-method}
\title{Count methods for RNAseq object}
\usage{
exonCounts(obj)
featureCounts(obj)
transcriptCounts(obj,from="exons")
geneCounts(obj,summarization=c("bestExons","geneModels"),...)
islandCounts(obj,force=FALSE,...)
}
\arguments{
\item{obj}{An object derived from class \code{\linkS4class{RNAseq}},can be
a \code{matrix} for RPKM, see details}

\item{force}{For \code{islandCount}, force RNAseq to redo \code{findIsland}}

\item{from}{either "exons" or "features" can be used to summarize per
transcript}

\item{summarization}{Method use for summarize genes}

\item{...}{See details}
}
\value{
A numeric vector containing count per exon, feature, gene or
transcript.
}
\description{
Summarize the read counts per exon, feature, gene, transcript or island.
\itemize{
\item{\code{exonCounts}: for that summarization, reads are
summarized per exons. An "exon" field is necessary in the annotation object
for this to work. See
\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq annotation
methods}} for more details on the annotation object.}
\item{\code{featureCounts} is similar to the 'exons' one. This is just a
wrapper to summarize count for genomic features that are not exon related.
I.e. one could use it to measure eRNAs. Again, a "feature" field is
necessary in the annotation object for this to work.}
\item{\code{geneCounts} sums the counts per either \code{bestExons} or
\code{geneModels}. In either case, the annotation object needs to contain
both an "exon" and a "gene" field.}
\item{\code{islandCounts} sums the
counts per computed islands.}
\item{\code{transcriptCounts} sums the counts
obtained by exons into their respective transcripts. Note that this often
result in counting some reads several times. For this function to work you
need both an "exon" and a "transcript" field in your annotation object. To
avoid this, one could create transcript specific synthetic exons, i.e.
features that would be unique to a transcript. To offer this possibility,
transcripts count can be summarized from "features", in which case the
annotation object need to have both the "feature" and "transcript" fields
defined.  }
}
}
\details{
\dots{} for \itemize{
\item{geneCounts: additional options for the
\code{\link[easyRNASeq:easyRNASeq-summarization-internal-methods]{.geneModelSummarization}}}
\item{islandCounts: additional options for
\code{\link[easyRNASeq:easyRNASeq-island-methods]{findIslands}} }}
}
\examples{
\dontrun{
library(BSgenome.Dmelanogaster.UCSC.dm3)

 # get the example data files
 tdir <- tutorialData()

 # get an example annotation file - we retrieve it from GitHub using curl
 gAnnot.rda <- fetchData("gAnnot.rda")

	# create an RNAseq object
	# summarizing 2 bam files by exons
	rnaSeq <- easyRNASeq(tdir,
                     organism="Dmelanogaster",
                     chr.sizes=seqlengths(Dmelanogaster),
                     readLength=36L,
                     annotationMethod="rda",
                     annotationFile=gAnnot.rda,
                     format="bam",
                     count="exons",
                     pattern="[A,C,T,G]{6}\\\\.bam$",
                     outputFormat="RNAseq")
	# summing up the exons by transcript
	rnaSeq <- transcriptCounts(rnaSeq)
}

}
\seealso{
\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq
annotation methods}}
\code{\link[easyRNASeq:easyRNASeq-summarization-internal-methods]{.geneModelSummarization}}
\code{\link[easyRNASeq:easyRNASeq-island-methods]{findIslands}}
}
\author{
Nicolas Delhomme
}
\keyword{methods}

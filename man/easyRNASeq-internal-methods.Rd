% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyRNASeq-internal-methods.R
\name{easyRNASeq internal methods}
\alias{easyRNASeq internal methods}
\alias{.extractIRangesList}
\alias{.catn}
\alias{.checkArguments}
\alias{.convertToUCSC}
\alias{.getArguments}
\alias{.getName}
\alias{.list.files}
\alias{.list.files-defunct}
\alias{.getWidth}
\alias{.normalizationDispatcher}
\alias{reduce}
\alias{reduce,RNAseq-method}
\alias{strand}
\alias{strand,RNAseq-method}
\alias{strand<-}
\alias{strand<-,RNAseq-method}
\alias{strand<-,RNAseq,ANY-method}
\title{Internal methods of RNAseq objects}
\usage{
.extractIRangesList(obj, chr.sel = c())
}
\arguments{
\item{obj}{An RNAseq object, or for the 'normalizationDispatcher',
depending on the type: a CountDataSet, a DGEList, a matrix, or an RNAseq
object respectively}

\item{chr.sel}{A list of chromosome to restrict the IRanges spaces
returned.}

\item{arg}{The argument name to check for.}

\item{chr.names}{The chromosome names, as a character vector, to be
converted to UCSC ones}

\item{fun}{The name of the function}

\item{organism}{The organism name}

\item{type}{character string specifying the type of object
(normalizationDispatcher)}

\item{value}{the appropriate strand object (strand and strand<-) or the
provided argument value (checkArguments)}

\item{x}{an object of the
\link[GenomicRanges:GRanges-class]{GenomicRanges},
\link[intervals:Intervals_virtual-class]{intervals} or
\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals}
package}

\item{...}{For \code{.getArguments} a list of named parameters to be
matched against a function formal definition. For \code{.catn}, the values
to be printed.}
}
\value{
\item{argString}{a character string representing these arguments
and their value that matched those defined in the formal definition of the
function}
\item{convertedChrNames}{a converted vector of chromosome names}
\item{i.range}{an IRange object} \item{names}{The
annotation names, i.e. a combination of exon, feature, transcript and gene}
\item{normalized.counts}{Depending on the type, a CountDataSet, a DGEList,
a NumericList, or NULL respectively}
}
\description{
These are generic internal methods:
\itemize{
\item Actual
\itemize{
\item .catn Just some pretty printing.
\item .checkArguments check that the provided argument match one
of the formal definition of the function. Stop if not.
\item .extractIRangesList extract an IRanges object from an AlignedRead or a
GAlignments object or a list returned by reading a bam file with
Rsamtools. It returns a list containing the IRangesList and library size.
\item .getArguments For a given function returns the arguments
passed as part of the \dots{} that match that function formals.
\item .getName Get the genomicAnnotation object names. Necessary to deal
with the different possible annotation object:
\code{GRanges} or \code{GRangesList}.
\item .getWidth Get the genomicAnnotation withs. Necessary to deal
with the different possible annotation object:
\code{GRanges} or \code{GRangesList}.
\item .normalizationDispatcher a function to dispatch
the normalization depending on the 'outputFormat' chosen by the user.
\item reduce Allow proper dispatch between the
\link[intervals:Intervals_virtual-class]{intervals} and the
\link[GenomicRanges:GRanges-class]{GenomicRanges} reduce function
\item strand Allow proper dispatch between the
\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} and
the \link[GenomicRanges:GRanges-class]{GenomicRanges} strand function
\item strand<- Allow proper dispatch between the
\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} and
the \link[GenomicRanges:GRanges-class]{GenomicRanges} strand replace
function
}
\item Defunct
\itemize{
\item .convertToUCSC
convert chromosome names to UCSC compliant ones.
\item .list.files check the arguments passed through the \dots to select
only the valid ones (defunct).
}
}
}
\author{
Nicolas Delhomme
}
\keyword{internal}

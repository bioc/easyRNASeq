\name{easyRNASeq internal methods}

\alias{.catn}
\alias{.checkArguments}
\alias{.convertToUCSC}
\alias{.doBasicCount}
\alias{.extendCountList}
\alias{.extractIRangesList}
\alias{.getArguments}
\alias{.getName}
\alias{.normalizationDispatcher}
\alias{reduce}
\alias{strand}
\alias{strand<-}


\title{Internal methods of RNAseq objects}

\description{
	These are generic internal methods:
	\itemize{
		\item{.catn} Just some pretty printing.
		\item{.checkArguments} check that the provided argument match one of the formal definition of the function. Stop if not.
		\item{.convertToUCSC} convert chromosome names to UCSC compliant ones.
		\item{.doBasicCount} A function to calculate the counts for 'exons' or 'features'
		\item{.extendCountList} extend or create the result count list of lists.
		\item{.extractIRangesList} extract an IRanges object from an AlignedRead or a GappedAlignments object or a list returned  by reading a bam file with Rsamtools.
		\item{.getArguments} For a given function returns the arguments passed as part of the \dots that match that function formals.
		\item{.getName} Get the genomicAnnotation object names. Necessary to deal with the different possible annotation object: \code{RangedData} or \code{GRangesList}.
		\item{.normalizationDispatcher} a function to dispatch the normalization depending on the 'outputFormat' chosen by the user.
		\item{reduce} Allow proper dispatch between the \link[intervals:Intervals_virtual-class]{intervals} and the \link[GenomicRanges:GRanges-class]{GenomicRanges} reduce function
		\item{strand} Allow proper dispatch between the \link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} and the \link[GenomicRanges:GRanges-class]{GenomicRanges} strand function
		\item{strand<-} Allow proper dispatch between the \link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} and the \link[GenomicRanges:GRanges-class]{GenomicRanges} strand replace function
	}
}

\usage{
	.catn(...)
	.checkArguments(fun,arg,value)
	.convertedChrNames <- convertToUCSC(chr.names,organism,chr.map)
	.counts <- doBasicCount(obj)
	.countList <- extendCountList(cList, values, type = c("exons", "transcripts", "genes", "islands"), subType = character(1),filename=character(1))
	i.range <- .extractIRangesList(aln,chr.sel=c())
	i.range <- .extractIRangesList(bam,chr.sel=c())
	i.range <- .extractIRangesList(gapped,chr.sel=c())
	argString <- .getArguments(fun,...)
	names <- .getName(obj)
	normalized.counts <- .normalizationDispatcher(obj,type=c("DESeq","edgeR","RPKM","RNAseq"),silent=FALSE,plot=TRUE,...)
}

\arguments{
	\item{aln}{A ShortRead AlignedRead object.}
	\item{arg}{The argument name to check for.}
	\item{bam}{A list as returned by the Rsamtools scanBam function.}
	\item{chr.map}{An optional parameter to provide the actual mappiung between the original and desired chromosome names. Only works with the 'organism' parameter set to the "custom" value. }
	\item{chr.names}{The chromosome names, as a character vector, to be converted to UCSC ones}
	\item{chr.sel}{A list of chromosome to restrict the IRanges spaces returned.}
	\item{cList}{list of lists that contain count results}
	\item{fun}{The name of the function}
	\item{gapped}{An object of the \link[GenomicRanges:GappedAlignments-class]{GappedAlignments} kind.}
	\item{obj}{An RNAseq object, or for the 'normalizationDispatcher', depending on the type:  a CountDataSet, a DGEList, a matrix, or an RNAseq object respectively}
	\item{organism}{The organism name}
	\item{plot}{The organism used to lookup UCSC chromosome names}
	\item{silent}{a boolean to decide whether or not to plot information}
	\item{silent}{a boolean to decide whether or not to print progress information}
	\item{subType}{character string defining a sub type of counts, i.e. for the gene type one of bestExon or geneModel}
	\item{type}{character string specifying the type of count ("exons", "transcripts", "genes" or islands) or the type of object (normalizationDispatcher)}
	\item{value}{the appropriate strand object (strand and strand<-) or the provided argument value (checkArguments)}
	\item{values}{a named vector containing count results}
	\item{x}{an object of the \link[GenomicRanges:GRanges-class]{GenomicRanges}, \link[intervals:Intervals_virtual-class]{intervals} or  \link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} package}
	\item{\dots}{For \code{getArguments} a list of named parameters to be matched against a function formal definition. For \code{catn}, the values to be printed.}
}

\value{
	\item{argString}{a character string representing these arguments and their value that matched those defined in the formal definition of the function}
	\item{convertedChrNames}{a converted vector of chromosome names}
	\item{counts}{an RleList}
	\item{countList}{an updated list of lists that contain count results}
	\item{i.range}{an IRange object}
	\item{names}{The annotation names, i.e. a combination of exon, feature, transcript and gene}
	\item{normalized.counts}{Depending on the type, a CountDataSet, a DGEList, a NumericList, or NULL respectively}
}

\author{Nicolas Delhomme}

\keyword{internal}


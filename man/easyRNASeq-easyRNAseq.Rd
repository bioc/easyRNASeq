\name{easyRNASeq}
\alias{easyRNASeq}
\alias{easyRNASeq,character-method}

\title{easyRNASeq method}
\description{
	This function is a wrapper around the more low level functionalities of the package.
	Is the easiest way to get a count matrix from a set of read files.
	It does the following:
	\itemize{
		\item{\code{\link[easyRNASeq:ShortRead-methods]{use ShortRead/Rsamtools methods}}} for loading/pre-processing the data.
		\item{\code{\link[easyRNASeq:fetchAnnotation]{fetch the annotations}}} depending on the provided arguments
		\item{\code{\link[easyRNASeq:fetchCoverage]{get the reads coverage}}} from the provided file(s)
		\item{\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarize the reads}}} according to the selected summarization features
		\item{\code{\link[easyRNASeq:easyRNASeq-correction-methods]{optionally apply}}} a data correction (i.e. generating RPKM).
		\item{\code{\link[easyRNASeq:edgeR-methods]{use edgeR methods}}} for post-processing the data or
		\item{\code{\link[easyRNASeq:DESeq-methods]{use DESeq methods}}} for post-processing the data (either of them being recommended over RPKM).
	}
}

\usage{
	easyRNASeq(filesDirectory=character(1),
	           organism=character(1),
            	   chr.sizes=c(),
	           readLength=integer(1),
	           annotationMethod=c("biomaRt","env","gff","gtf","rda"),
	           annotationFile=character(1),
		   annotationObject=RangedData(),
	           format=c("aln","bam"),gapped = FALSE,
	           count=c('exons','features','genes','islands','transcripts'),
	           outputFormat=c("DESeq","edgeR","matrix","RNAseq"),
	           pattern=character(1),filenames=character(0),nbCore=1,
	           filter=srFilter(),type="SolexaExport",
	           chr.sel=c(),summarization=c("bestExons","geneModels"),
	           normalize=FALSE,max.gap=integer(1),min.cov=1L,
	           min.length=integer(1),plot=TRUE,
		   conditions=c(),validity.check=TRUE,
		   chr.map=data.frame(),ignoreWarnings=FALSE,
             	   silent=FALSE,...)
}

\arguments{
	\item{annotationFile}{The location (full path) of the annotation file}
	\item{annotationObject}{A \code{\linkS4class{RangedData}} or \code{\linkS4class{GRangesList}} object containing the annotation.}
	\item{annotationMethod}{The method to fetch the annotation, one of "biomaRt","env","gff","gtf" or "rda". All methods but "biomaRt" and "env" require the annotationFile to be set. The "env" method requires the annotationObject to be set.}
	\item{chr.map}{A data.frame describing the mapping of original chromosome names towards wished chromosome names. See details.}
	\item{chr.sel}{A vector of chromosome names to subset the final results.}
	\item{chr.sizes}{A vector or a list containing the chromosomes' size of the selected organism}
	\item{conditions}{A vector of descriptor, each sample must have a descriptor if you use outputFormat DESeq or edgeR. The size of this list must be equal to the number of sample. In addition the vector should be named with the filename of the corresponding samples.}
	\item{count}{The feature used to summarize the reads. One of 'exons','features','genes','islands' or 'transcripts'. See details.}
	\item{filenames}{The name, not the path, of the files to use}
	\item{filesDirectory}{The directory where the files to be used are located}
	\item{filter}{The filter to be applied when loading the data using the "aln" format}
	\item{format}{The format of the reads, one of "aln","bam". If not "bam", all the types supported by the \pkg{ShortRead} package are supported too.}
	\item{gapped}{Is the bam file provided containing gapped alignments?}
	\item{ignoreWarnings}{set to TRUE (bad idea! they have a good reason to be there) if you do not want warning messages.}
	\item{min.cov}{When computing read islands, the minimal coverage to take into account for calling an island}
	\item{min.length}{The minimal size an island should have to be kept}
	\item{max.gap}{When computing read islands, the maximal gap size allowed between two islands to merge them}
	\item{nbCore}{defines how many CPU core to use when computing the geneModels. Use the default parallel library}
	\item{normalize}{A boolean to convert the returned counts in RPKM. Valid when the \code{outputFormat} is left undefined (i.e. when a matrix is returned) and when it is \code{DESeq} or \code{edgeR}. Note that it is not advised to normalize the data prior DESeq or edgeR usage!}
	\item{organism}{A character string describing the organism}
	\item{outputFormat}{By default, easyRNASeq returns a count matrix. If one of \code{DESeq},\code{edgeR},\code{RNAseq} is provided then the respective object will be returned.}
	\item{pattern}{For easyRNASeq, the pattern of file to look for, e.g. "bam$"}
	\item{plot}{Whether or not to plot assessment graphs.}
	\item{readLength}{The read length in bp}
	\item{silent}{set to TRUE if you do not want messages to be printed out.}
	\item{summarization}{A character defining which method to use when summarizing reads by genes. So far, only "geneModels" is available.}
	\item{type}{The type of data when using the "aln" format. See the ShortRead library.}
	\item{validity.check}{Shall UCSC chromosome name convention be enforced}
	\item{\dots}{additional arguments. See details}
}

\value{
	Returns a count table (a matrix of m features x n samples) unless the \code{outputFormat} option has been set, in which case an object of type \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}} or \code{\link[edgeR:DGEList]{edgeR:DGEList}} or \code{\linkS4class{RNAseq}} is returned
}

\details{
  \itemize{
  \item{\dots for the easyRNASeq call}
  Additional arguments, passed to the \pkg{biomaRt} \code{\link[biomaRt:getBM]{getBM}} function or
  to the \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}} internal function that takes an optional
  arguments: annotation.type that default to "exon" (used to select the proper rows of the gff or gtf file) or to the \code{\link[DESeq:estimateDispersions]{DESeq estimateDispersions}} method.
  \item{the annotationObject}
  When the \code{annotationMethods} is set to \code{env} or \code{rda}, a properly formatted \code{RangedData} or \code{GRangesList} object need to be provided. Check the paragraph RangedData in the vignette or the examples at the bottom of this page for examples. The data.frame-like structure of these objects is where \code{easyRNASeq} will look for the exon, feature, transcript, or gene identifier. Depending on the count method selected, it is essential that the akin column name is present in the annotationObject. E.g. when counting "features", the annotationObject has to contain a "feature" field.
  \item{the chr.map}
  The chr.map argument for the easyRNASeq function only works for an "organismName" of value 'custom' with the "validity.check" parameter set to 'TRUE'. 
  This data.frame should contain two columns named 'from' and 'to'. The row should represent the chromosome name in your original data and the wished name in the output of the function.
  \item{count}
  The count can be summarized by exons, features, genes, islands or transcripts. While exons, genes and transcripts are obvious, "features" describes any features provided by the user, e.g. enhancer loci. These are processed as the exons are. For "islands", it is for an under development function that identifies de-novo expression loci and count the number of reads overlapping them.
}
}

\author{Nicolas Delhomme}

\seealso{
	\code{\linkS4class{RNAseq}}
	\code{\link[edgeR:DGEList]{edgeR:DGEList}}
	\code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}
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

	## an example of a chr.map
	chr.map <- data.frame(from=c("2L","2R","MT"),to=c("chr2L","chr2R","chrMT"))

	## an example of a RangedData annotation
	gAnnot <- RangedData(
                     IRanges(
                             start=c(10,30,100),
                             end=c(21,53,123)),
                          space=c("chr01","chr01","chr02"),
                          strand=c("+","+","-"),
                          transcript=c("trA1","trA2","trB"),
                          gene=c("gA","gA","gB"),
                          exon=c("e1","e2","e3"),
                          universe = "Hs19"
                          )

	## an example of a GRangesList annotation
	grngs <- as(gAnnot,"GRanges")
	grngsList<-split(grngs,seqnames(grngs))
}
}

\keyword{methods}

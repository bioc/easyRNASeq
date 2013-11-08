## to re-build the documentation and namespace
## library(roxygen2)
## roxygenize("easyRNASeq/")
##' Count summarization and normalization pipeline for Next Generation Sequencing data.
##'
##' Offers functionalities to summarize read counts per feature of interest, e.g. exons, transcripts, genes, etc.
##' Offers functionalities to normalize the summarized counts using 3rd party packages like \code{\link[DESeq:newCountDataSet]{DESeq}}
##' or \code{\link[edgeR:DGEList]{edgeR}}.
##'
##' \tabular{ll}{
##' Package: \tab easyRNASeq\cr
##' Type: \tab Package\cr
##' Version: \tab 1.5.1\cr
##' Date: \tab 2012-10-15\cr
##' License: \tab Artistic-2.0\cr
##' LazyLoad: \tab yes\cr
##' Depends: \tab methods, parallel, biomaRt, edgeR, DESeq, genomeIntervals, Rsamtools, ShortRead, RnaSeqTutorial\cr
##' Suggests: \tab BSgenome.Dmelanogaster.UCSC.dm3
##' }
##'
##' @section Methods:
##' The main function \code{\link[easyRNASeq:easyRNASeq]{easyRNASeq}} will summarize the counts per
##' feature of interest, for as many samples as provided and will return a
##' count matrix (N*M) where N are the features and M the samples.
##' This data can be corrected to \pkg{RPKM} in which case
##' a matrix of corrected value is returned instead, with the same dimensions.
##' Alternatively a \code{\linkS4class{SummarizedExperiment}} can be returned and this
##' is expected to be the default in the upcoming version of easyRNASeq (as of 1.5.x).
##' If the necessary sample
##' information are provided, the data can be normalized using either \code{\link[DESeq:newCountDataSet]{DESeq}}
##' or \code{\link[edgeR:DGEList]{edgeR}} and the corresponding package object returned.
##' For more insider details, and step by step functions, see:
##' \tabular{ll}{
##' 	\code{\link[easyRNASeq:ShortRead-methods]{ShortRead methods}} for pre-processing the data.
##' 	\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq annotation methods}} for getting the annotation.
##' 	\code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq coverage methods}} for computing the coverage from a Short Read Alignment file.
##' 	\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}} for summarizing the data.
##' 	\code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq correction methods}} for correcting the data (i.e. generating RPKM).
##' 	\code{\link[easyRNASeq:edgeR-methods]{edgeR methods}} for post-processing the data.
##' 	\code{\link[easyRNASeq:DESeq-methods]{DESeq methods}} for post-processing the data.
##' 	}
##' 
##' @name easyRNASeq package
##' @rdname easyRNASeq-package
##' @aliases easyRNASeq-package
##' @docType package
##' @author Nicolas Delhomme
##' @keywords package
##' @seealso
##' 	The class RNAseq specification:
##' 	\code{\linkS4class{RNAseq}}
##'
##'     The default output class specification:
##'     \code{\linkS4class{SummarizedExperiment}}
##' 
##' 	The imported packages:
##' 	\code{\link[biomaRt:useMart]{biomaRt}}
##' 	\code{\link[edgeR:DGEList]{edgeR}}
##' 	\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals}}
##' 	\code{\link[Biostrings:XString-class]{Biostrings}}
##' 	\code{\link[BSgenome:BSgenome-class]{BSgenome}}
##' 	\code{\link[DESeq:newCountDataSet]{DESeq}}
##' 	\code{\link[GenomicRanges:GRanges-class]{GenomicRanges}}	
##' 	\code{\link[IRanges:IRanges-constructor]{IRanges}}
##' 	\code{\link[Rsamtools:scanBam]{Rsamtools}}
##' 	\code{\link[ShortRead:readAligned]{ShortRead}}
##' 
##' 	The suggested packages:
##' 	\code{\link[parallel:makeCluster]{parallel}}
##' 	\code{\link[GenomicFeatures:TranscriptDb-class]{GenomicFeatures}}
##' 
##' @examples
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	library(BSgenome.Dmelanogaster.UCSC.dm3)
##' 
##' 	## creating a count table from 4 bam files
##' 	count.table <- easyRNASeq(filesDirectory=
##' 		    			system.file(
##' 					"extdata",
##' 					package="RnaSeqTutorial"),
##' 					pattern="[A,C,T,G]{6}\.bam$",
##' 				format="bam",
##' 				readLength=36L,
##' 				organism="Dmelanogaster",
##' 				chr.sizes=as.list(seqlengths(Dmelanogaster)),
##' 				annotationMethod="rda",
##' 				annotationFile=system.file(
##' 				                            "data",
##' 							    "gAnnot.rda",
##' 							    package="RnaSeqTutorial"),
##' 				count="exons")
##' 
##' 	}
##' 
NULL

## To define the NAMESPACE
## library(codetoolsBioC)
## library(easyRNASeq)
## writeNamespaceImports("easyRNASeq")
## import classes
##' @importClassesFrom biomaRt Mart
##' @importClassesFrom Biostrings DNAStringSet
##' @importClassesFrom DESeq CountDataSet
##' @importClassesFrom edgeR DGEList
##' @importClassesFrom genomeIntervals Genome_intervals
##' @importClassesFrom GenomicRanges GenomicRanges GRangesList
##' @importClassesFrom IRanges RangedData RleList Vector
##' @importClassesFrom methods ANY character integer
##' list matrix numeric vector
##' @importClassesFrom Rsamtools ScanBamParam
##' @importClassesFrom ShortRead AlignedRead ShortReadQ SRFilter
## import S4 methods
##' @importMethodsFrom Biobase fData
##' @importMethodsFrom BiocGenerics cbind colnames duplicated
##' eval intersect lapply order paste rbind rownames sapply
##' table unique
##' @importMethodsFrom Biostrings initialize
##' @importMethodsFrom DESeq counts estimateDispersions estimateSizeFactors
##' @importMethodsFrom genomeIntervals seq_name
##' @importMethodsFrom GenomicRanges rname strand "strand<-"
##' @importMethodsFrom IRanges "%in%" aggregate as.list as.matrix coerce
##' "colnames<-" coverage elementMetadata
##' "elementMetadata<-" end findOverlaps gsub ifelse
##' match mean median na.omit narrow
##' nchar ncol nrow queryHits ranges reduce rev
##' "rownames<-" runLength runsum runValue show
##' sort split start sub substr tolower
##' "universe<-" unlist values which width
##' @importMethodsFrom Rsamtools scanBam ScanBamParam
##' @importMethodsFrom ShortRead chromosome position readAligned
##' srdistance srFilter
## import methods
##' @importFrom biomaRt getBM listDatasets useDataset useMart
##' @importFrom Biostrings DNAStringSet
##' @importFrom DESeq fitInfo newCountDataSet
##' @importFrom edgeR calcNormFactors DGEList estimateCommonDisp
##' estimateTagwiseDisp maPlot plotMDS.DGEList plotMeanVar
##' @importFrom genomeIntervals getGffAttribute readGff3
##' @importFrom GenomicAlignments readGAlignments
##' @importFrom graphics abline axis boxplot grid hist legend lines
##' mtext par plot rect
##' @importFrom IRanges IRanges IRangesList isSingleString LogicalList
##' RangedData RangesList SplitDataFrameList
##' @importFrom methods as extends is new
##' @importFrom parallel makePSOCKcluster parLapply stopCluster
##' @importFrom Rsamtools scanBamFlag
##' @importFrom ShortRead alignData sread
##' @importFrom utils combn str
NULL




## to re-build the documentation and namespace
## library(roxygen2)
## roxygenize("easyRNASeq/")

## to update the package versions
# pkg <- c("Biobase","BiocGenerics","biomaRt","Biostrings",
#          "DESeq","edgeR","GenomeInfoDb","genomeIntervals",
#          "GenomicAlignments","GenomicRanges",
#          "IRanges","LSD","Rsamtools","S4Vectors","ShortRead",
#          "BiocStyle","BSgenome",
#          "BSgenome.Dmelanogaster.UCSC.dm3","GenomicFeatures",
#          "RnaSeqTutorial","RUnit")
# pkg[! pkg %in% rownames(installed.packages())]
# installed.packages()[pkg,"Version"]

###==========================
### package details
###==========================
##' Count summarization and normalization pipeline for Next Generation Sequencing data.
##'
##' Offers functionalities to summarize read counts per feature of interest, e.g. exons, transcripts, genes, etc.
##' Offers functionalities to normalize the summarized counts using 3rd party packages like \code{\link[DESeq:newCountDataSet]{DESeq}}
##' or \code{\link[edgeR:DGEList]{edgeR}}.
##'
##' \tabular{ll}{
##' Package: \tab easyRNASeq\cr
##' Type: \tab Package\cr
##' Version: \tab 2.3.1\cr
##' Date: \tab 2015-02-11\cr
##' License: \tab Artistic-2.0\cr
##' LazyLoad: \tab yes\cr
##' Depends: \tab methods, parallel, Biobase, BiocGenerics, biomaRt, Biostrings, edgeR, DESeq, genomeIntervals, GenomeInfoDb, GenomicAlignments, GenomicRanges, graphics, IRanges, LSD, Rsamtools, S4Vectors, ShortRead, utils\cr
##' Suggests: \tab BSgenome.Dmelanogaster.UCSC.dm3
##' }
##'
##' @section Methods:
##' The main function \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} will summarize the counts per
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
##' @aliases easyRNASeq-package assay type BamFileList BamFileList-class IRanges
##' RangedData SRFilterResult SummarizedExperiment-class chromosomeFilter
##' compose nFilter RangedData-class
##' @docType package
##' @author Nicolas Delhomme, Bastian Schiffthaler, Ismael Padioleau
##' @keywords package
##' @seealso
##' 	The class RNAseq specification:
##' 	\code{\linkS4class{RNAseq}}
##'
##'   The default output class specification:
##'   \code{\linkS4class{SummarizedExperiment}}
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
##' 	\code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}}
##'
##'   The following classes and functions that are made available from
##'   other packages:
##'   \itemize{
##'     \item{Classes}{
##'       \code{\linkS4class{BamFileList}}
##'       \code{\linkS4class{RangedData}}
##'       \code{\linkS4class{SummarizedExperiment}}
##'     }
##'     \item{Functions/Methods}{
##'       \code{\link[GenomicRanges:SummarizedExperiment-class]{
##'         The SummarizedExperiment assay accessor}
##'       }
##'       \code{\link[Rsamtools:BamFileList]{The BamFileList constructor}}
##'       \code{\link[IRanges:IRanges-constructor]{The IRanges constructor}}
##'       \code{\link[IRanges:RangedData-class]{The RangedData constructor}}
##'       \code{\link[ShortRead:srFilter]{For the SRFilterResult,
##'       chromosomeFilter, compose and nFilter methods}}
##'     }
##'   }
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
##' 					pattern="[A,C,T,G]{6}\\.bam$",
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

###==========================
## To define the NAMESPACE
###==========================
## library(codetoolsBioC)
## library(easyRNASeq)
## writeNamespaceImports("easyRNASeq")
## import classes
##' @importClassesFrom Biostrings DNAStringSet
##' @importClassesFrom DESeq CountDataSet
##' @importClassesFrom edgeR DGEList
##' @importClassesFrom genomeIntervals Genome_intervals
##' @importClassesFrom GenomicAlignments GAlignments GAlignmentPairs
##' @importClassesFrom GenomicRanges GRanges GRangesList SummarizedExperiment
##' @importClassesFrom IRanges RangedData RleList
##' @importClassesFrom S4Vectors Vector DataFrame SimpleList
##' @importClassesFrom methods ANY character "function" integer
##' list matrix missing numeric vector
##' @importClassesFrom Rsamtools BamFile BamFileList
##' @importClassesFrom ShortRead AlignedRead ShortReadQ
## import S4 methods
##' @importMethodsFrom Biobase fData varMetadata
##' @importMethodsFrom BiocGenerics annotation cbind clusterApply
##' colnames counts duplicated estimateDispersions estimateSizeFactors
##' eval get intersect lapply match order paste pmax rbind rownames sapply
##' strand "strand<-" table unique
##' @importMethodsFrom Biostrings type
##' @importMethodsFrom genomeIntervals seq_name
##' @importMethodsFrom GenomeInfoDb seqinfo seqlengths "seqlengths<-"
##' seqlevels "seqlevels<-" seqnames "seqnames<-"
##' @importMethodsFrom GenomicAlignments cigar summarizeOverlaps
##' @importMethodsFrom GenomicRanges assay assays "assays<-" colData "colData<-"
##' "exptData<-" grglist rowData "rowData<-" SummarizedExperiment
##' @importMethodsFrom IRanges aggregate as.list as.matrix as.table
##' "colnames<-" countOverlaps coverage elementLengths end "end<-" findOverlaps
##' gsub mean median narrow nchar queryHits ranges reduce rev "rownames<-" space
##' split start "start<-" sub  tolower "universe<-" unlist values which width
##' @importMethodsFrom S4Vectors "%in%" elementMetadata "elementMetadata<-"
##' ifelse levels mcols Rle runLength runsum runValue substr
##' @importMethodsFrom methods coerce initialize show
##' @importMethodsFrom Rsamtools countBam path scanBam scanBamHeader
##' ScanBamParam yieldSize "yieldSize<-"
##' @importMethodsFrom S4Vectors elementMetadata "elementMetadata<-" mcols
##' @importMethodsFrom ShortRead chromosome id position readAligned
##' srdistance sread srFilter writeFastq
## import methods
##' @importFrom biomaRt getBM listDatasets useDataset useMart
##' @importFrom Biostrings DNAStringSet
##' @importFrom DESeq fitInfo newCountDataSet
##' @importFrom edgeR calcNormFactors DGEList estimateCommonDisp
##' estimateTagwiseDisp maPlot plotMDS.DGEList plotMeanVar plotBCV
##' @importFrom genomeIntervals getGffAttribute parseGffAttributes readGff3
##' @importFrom GenomicAlignments GAlignments readGAlignments readGAlignmentPairs
##' @importFrom GenomicRanges GRanges GRangesList
##' @importFrom graphics abline axis axTicks boxplot grid hist legend lines
##' mtext par plot rect
##' @importFrom IRanges IRanges IRangesList LogicalList
##' RangedData RangesList SplitDataFrameList RleList
##' @importFrom LSD heatscatter
##' @importFrom methods as extends is new
##' @importFrom parallel makePSOCKcluster parLapply stopCluster
##' @importFrom Rsamtools BamFileList bamFlagTest index scanBamFlag
##' @importFrom S4Vectors DataFrame SimpleList
##' @importFrom ShortRead alignData chromosomeFilter compose nFilter
##' SRFilterResult
##' @importFrom utils combn str
## and export!
##' @exportClass BamFileList RangedData SummarizedExperiment
##' @exportMethod assay seqnames split srFilter SummarizedExperiment width writeFastq
##' @export chromosomeFilter compose BamFileList IRanges nFilter RangedData readAligned SRFilterResult
NULL

###==========================
## To detail the deprecation
###==========================
##' The following function have been deprecated:
##' \itemize{
##' \item \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}}
##' \item \code{\link[easyRNASeq:easyRNASeq-coverage-methods]{fetchCoverage}}
##' \item \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{fetchAnnotation}}
##' \item \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{knownOrganisms}}
##' }
##'
##' \itemize{
##' \item The \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} function is superseded by the
##' \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq}} function to consolidate and
##' prune the overall package. The changes are based on user comments and on the
##' general standardization occuring in the field.
##' \item The \code{\link[easyRNASeq:easyRNASeq-coverage-methods]{fetchCoverage}} function only had two
##' parameters deprecated as the consequence of the package consolidation. As the scanBam
##' function is not called directly anymore but through higher level functions (from the
##' GenomicRanges package), the 'what' and 'isUnmappedQuery' parameters were obsolete.
##' }
##' @aliases easyRNASeq easyRNASeq,RNAseq-method
##' fetchCoverage,RNAseq-method fetchCoverage
##' knownOrganisms knownOrganisms,missing-method
##' fetchAnnotation fetchAnnotation,RNAseq-method
##' @name Deprecated functions
##' @rdname easyRNASeq-deprecated
NULL

###==========================
## To detail defunct function
###==========================
##' The following function are defunct:
##' \itemize{
##' \item \code{lotDispersionEstimates,DGEList-method}
##' }
##'
##' \itemize{
##' \item The \code{plotDispersionEstimates,DGEList-method}
##' function is superseded by the \code{\link[edgeR:plotBCV]{plotBCV}} function
##' as the \pkg{edgeR} DGEList object structure changed
##' }
##' @aliases plotDispersionEstimates,DGEList-method
##' @name Defunct functions
##' @rdname easyRNASeq-defunct
NULL

###==========================
## To detail hardcoded parameters
###==========================
##' Objects created when the package is attached.
##'
##' The package creates the following objects when attached
##' \itemize{
##' \item GTF.FIELDS
##' \item ANNOTATION.TYPE
##' }
##'
##' These objects hold the following information
##' \itemize{
##' \item GTF.FIELDS \code{c("gene_id","transcript_id","exon_id","gene_name")}
##' \item ANNOTATION.TYPE \code{c(mRNA="mRNA",exon="exon")}
##' }
##' and are designed as global variables to expose the
##' fact that they are hardcoded. There exist as
##' placeholder in case a user would require different
##' values for these.
##'
##' @aliases easyRNASeq-global-variables .onAttach ANNOTATION.TYPE GTF.FIELDS
##' @name easyRNASeq-global-variables
##' @rdname easyRNASeq-global-variables
##' @param libname a character string giving the library directory where the package defining the namespace was found.
##' @param pkgname a character string giving the name of the package.
##' @seealso \code{\link[base:ns-hooks]{.onAttach}} in the \code{base} package.
##' @keywords internal
".onAttach" <- function(libname,pkgname){
  assign("GTF.FIELDS",c("gene_id","transcript_id","exon_id",
                        "gene_name"),
         envir=as.environment("package:easyRNASeq"))
  assign("ANNOTATION.TYPE",c(mRNA="mRNA",exon="exon"),
         envir=as.environment("package:easyRNASeq"))
}
